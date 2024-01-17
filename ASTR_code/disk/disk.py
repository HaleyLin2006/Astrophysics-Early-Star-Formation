import numpy as np
from tqdm import tqdm
import math
from scipy.interpolate import griddata
import os
import sys

sys.path.insert(1, '/Users/haleylin/Desktop/ASTR Sphinx/ASTR_code/')
sys.path.insert(2, '/Users/haleylin/Desktop/ASTR Sphinx/ASTR_code/star')
sys.path.insert(1, '/Users/haleylin/Desktop/ASTR Sphinx/ASTR_code/envelope')

from parameter import param
from data_handler import data
from star.star_help import star
from envelope.envelope_help import envelope

param = param()
data = data()
star = star()
envelope = envelope()

'''
    QUESTION: sigma_0 -> when the disk grow, r_cent grow, but

                    sigma_0 = M_d / (4π (r_cent-r_i) ) 

            - M_d grows slower than r_cent
            - unit of sigma_0 should be [M]/[L]^2, but here, it is [M]/[L]

    QUESTION: Integration of disk mass should be using arclength to calculate the area

                dA = 2π R sqrt(dx^2 + dy^2) dx dy

                dMdt = density * velocity * dA

    QUESTION: Because sigma_0 is calculated wrongly, the disk spread wrongly and boundary
                is wrong

    QUESTION: Luminosity is fixed for now
'''
class disk():
    def __init__(self):
        data_param = ["CM_Cartesian", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        if os.path.exists(data._get_directory(data_param)):
            print('Loading CM81 Cartesian data...')

            self.radius, self.r_grid, self.cm_cgs_density, self.cm_cgs_v_r, self.cm_cgs_v_theta, cgs_v_phi = data.load(data_param)

            # radius in [cm]
            print("radius",     len(self.radius),          len(self.radius[0]),          len(self.radius[0]))                    # [time][radius]
            print("density",    len(self.cm_cgs_density),  len(self.cm_cgs_density[0]),   len(self.cm_cgs_density[0][0]))        # [time][theta][radius]
            print("v_r",        len(self.cm_cgs_v_r),      len(self.cm_cgs_v_r[0]),       len(self.cm_cgs_v_r[0][0]))            # [time][theta][radius]
            print("v_theta",    len(self.cm_cgs_v_theta),  len(self.cm_cgs_v_theta[0]),   len(self.cm_cgs_v_theta[0][0]))        # [time][theta][radius]

            print(np.array(self.radius[0]))

            # because during cm81 calculation, I forgot to take abs(v_r)...
            # self.cm_cgs_density = np.absolute(np.array(self.cm_cgs_density)).tolist()
            # print(np.array(self.cm_cgs_density).shape)

        else:
            print('CM81 Cartesian Data Does Not exist!')
            self._initialize_CM_cartesian()

        # print("Saving data")
        # data_param = ["CM_Cartesian", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        # data.save(param=data_param, radius=self.radius, r_grid=self.r_grid, density=self.cm_cgs_density, v_r=self.cm_cgs_v_r, v_theta=self.cm_cgs_v_theta, v_phi=cgs_v_phi)
        # print("Data saved")

    def disk_evol(self, change_luminosity=False):
        '''
        1. Check boundary
            where disk density and envelope density match
                disk density = sigma/(h * np.sqrt*(2*π)) * np.exp(-z**2/(2*h**2))
                sigma(t_0) = 1e-5

            if this bondary is greater than 3 grid points, start considering the next step, integrating
                3 grid pts from R_i = sqrt(L_* / (4π sigma T_evap^4))

                (does not depends on centrifugal radius at all??)

            (Fortran) Take the vertical profile of the disk into account for computing
                the mass loading. First, compute the "surface" of the disk: the
                point where disk density and envelope density are equal. If no
                such point exists for a given R, place the surface at 1e90 cm.
            
                If the disk density at the midplane is lower than the envelope
                density, the disk may actually still exist at this radius. So,
                scan upwards to see if the disk density exceeds the envelope
                density at any other z.

        2. Within the boundary, integrate to get M_d
            M_d = 2π S_{boundary} rho(R, z) dV
        
        3. Trace back surface density
            M_d = 2π S_{boundary} R Sigma dA
            With extrapolation trick, Sigma(R = Sigma(R_i)*(R_i/R), given arbitrary R_i
            
            (can we directly write rho(R,z) dV = R Sigma dA ?)
            (in what direction this extrapolation work?)

        4. Interaction betweem disk and star: M_*, L_*
            dM_{d -> *}dt = -2πR_i * u_{R_i} * Sigma_{R_i}
            M_* = S_{0 -> t} dM_{d->*}dt + dM_{env->*} dt
            L_*(M_*)...

        5. Update L_* at the first step, and repeat...
            (when to stop?)
        '''
        self.M_disk = []
        
        # fix the luminosity for now
        r_star = star.get_r_star()
        star_mass = envelope.cgs_env_to_star_mass(r_star)
        self.r_i = self._extrapolate_ri(star.disk_inner_radius(star_mass, const_lum=True))
        self.centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 for t in param.time]
        
        # print(["{:e}".format(i) for i in self.r_i])
        # print(["{:e}".format(i) for i in self.centrifugal_radius])
        # print(["{:e}".format(i) for i in self.radius[0]])zero_height = True
        first_time = True
        for t in tqdm(range(len(param.time)), desc='Disk evloving through time...'):
            # if the disk is smaller than 3 grid points, than don't initiate this progress
            if_start_integration = self.if_start_integration(self.r_i[t], self.centrifugal_radius[t])
            if not if_start_integration:
                print(f"At time {t+1}/100 Disk not yet form")
                self.M_disk.append(param.min_disk_mass)
                continue

            # The sigma constant for the disk to spread out
            
            if first_time:
                h = self._find_closest_indices(self.radius[t], [0.1*r for r in self.radius[t]])
                sigma_0 = self.get_sigma_0(m_d=self.M_disk[-1], r_i=self.r_i[t], r_cent=self.centrifugal_radius[t], h=h, sigma_r=[], first=True, t=t)
            else:
                sigma_0 = self.get_sigma_0(m_d=self.M_disk[-1], r_i=self.r_i[t], r_cent=self.centrifugal_radius[t], h=h, sigma_r=sigma_r, first=False, t=t)
            
            sigma_r = [sigma_0 * (self.r_i[t]/r) for r in self.radius[t]]
            print("sigma_0", sigma_0)

            # determine the boundary of the disk with density condition
            dh = 0
            start = True
            while dh > 3 or start:
                start = False

                past_h = h
                h = self.check_boundary_w_density(sigma_r=sigma_r, h=past_h, t=t)
                if self._check_zero_height(h):
                    h = self._find_closest_indices(self.radius[t], [0.1*r for r in self.radius[t]])
                # change it to fortran's model, when the height is so small, directly use CM81 for disk density

                dh = self._highest_dh(past_h, h)
                print("Run", t, dh)
            
            print(np.array(h))

            dm_i = self.get_next_m_d(h=h, t=t)
            print("prev mass", self.M_disk[-1], "dm", dm_i)
            m_d_i = self.M_disk[-1] + dm_i
            self.M_disk.append(m_d_i)

            first_time = False

        return self.M_disk

    def get_sigma_0(self, m_d, r_i, r_cent, h, sigma_r, first, t):
        '''
            Return sigma_0, the normalization parameter that governs the distribution of surface density (sigma).

            Parameters:
                m_d     (float) : mass of disk for current interation
                r_i     (float) : the inner radius of the disk by evaporation, dependent on luminosity of star
                r_cent  (float) : the centrifugal radius of the disk 
                h       (list)  : index of the height of the disk, initially we have thin disk assumption H = 0.1R
                t       (int)   : the index of param.time specifying when we are dealing the issue
            Returns:
                sigma_0 (float) : normalization parameter that governs the surface density (sigma) distribution
                                    at R direction
                                    Sigma(R) = sigma_0 (R/R_0)^{-1} 
        '''
        if first: 
            return m_d / (4*np.pi * r_i * (r_cent-r_i))

        h = self.radius[t][h]

        # print(m_d, r_i, r_cent, self._find_closest_index(self.radius[0], r_cent))
        r_boundary = 0
        ri = self._find_closest_index(self.radius[t], r_cent)
        for r in range(len(self.radius[t])):
            env_density  = self.cm_cgs_density[t][r][0]
            disk_density = sigma_r[r] / (h[r] * np.sqrt(2*np.pi))
            # print(r, "env", env_density, "disk", disk_density)
            if disk_density < env_density:
                # print(r, disk_density, env_density)
                r_boundary = self.radius[t][r]
                ri = r
                break
        
        # print(np.array([self.cm_cgs_density[t][r][0] for r in range(len(self.radius[t]))]).shape, np.array([sigma_r[r] / (h[r] * np.sqrt(2*np.pi)) for r in range(len(self.radius[t]))]).shape)
        print(ri, self._find_closest_index(self.radius[t], r_cent))

        if r_boundary > r_cent or r_boundary <= r_i:
            print(m_d, r_i, (r_cent-r_i))
            sigma_0 = m_d / (4*np.pi * r_i * (r_cent-r_i))
            print("Check Sigma_0", m_d, sigma_0*(4*np.pi * r_i * (r_cent-r_i)))
        if r_boundary <= r_cent and r_boundary > r_i:
            print(m_d, r_i, (r_boundary-r_i))
            sigma_0 = m_d/ (4*np.pi * r_i * (r_boundary-r_i))
            print("Check Sigma_0", m_d, sigma_0*(4*np.pi * r_i * (r_boundary-r_i)))
        
        return sigma_0

    def check_boundary_w_density(self, sigma_r, h, t):
        '''
            Check the boundary of the disk with density method. This is by identifying at which height does the density structure
            of the disk match with the density of the envelope.

            Parameters:
                sigma_r (list)  : the sigma at the midplane throughout the radius
                h       (list)  : index of the height of the disk, initially we have thin disk assumption H = 0.1R
                t       (int)   : the index of param.time specifying when we are dealing the issue
            Returns:
                height  (list)  : the updated height of the disk according to the density structure
        '''
        h = self.radius[t][h]

        height = []
        for r in range(len(self.radius[t])):
            if self.radius[t][r] >= self.centrifugal_radius[t]:
                height.append(0)
                continue
            # if r > len(h)-1:
            #     break

            rho_d_midplane = sigma_r[r] / (h[r]*np.sqrt(2*np.pi))

            # print(sigma_r[r], (h[r]*np.sqrt(2*np.pi)), self.cm_cgs_density[t][r][0])

            # print(rho_d_midplane, abs(self.cm_cgs_density[t][r][0]))
            if rho_d_midplane < abs(self.cm_cgs_density[t][r][0]):               # here we assume that it is in [time][x][y] format
                # print("continue")
                height.append(0)
                continue

            height_index = 0
            for z in range(len(self.radius[t])):
                # print(rho_d_midplane, (-self.radius[t][z]**2/(2*h[r]**2)))
                rho_d_z = rho_d_midplane * np.exp(-self.radius[t][z]**2/(2*h[r]**2))
                density = self.cm_cgs_density[t][r][z]

                # print(rho_d_z, density)

                if density > rho_d_z:
                    height_index = z
                    # print("Height")
                    break
            
            height.append(height_index)
        
        # print(height)

        # at early time, the result is all 0, but indeed the central plane density is higher than envelope density...
        # print("h", [i/param.AU for i in height])
        return height

    def if_start_integration(self, r_i, r_cent):
        '''
            If the distance between the inner radius of the disk (r_i) and the centrifugal radius is greater
            than 3 grid point, then the disk is too small and we wish skip to the next time instance.

            This function detects such instance where the disk is too small.

            Parameters:
                r_i     (float) : the inner radius of the disk by evaporation, dependent on luminosity of star
                r_cent  (float) : centrifugal radius
            Returns:
                start_int (bool): True, if r_cent-r_i is greater than 3 grid points, False otherwise.
        '''
        r_i_index   = self._find_closest_index(self.radius[0], r_i)
        r_cent_index = self._find_closest_index(self.radius[0], r_cent)

        # print(["{:e}".format(i) for i in self.radius[0]])
        # print("{:e}".format(r_i/param.AU), "{:e}".format(r_cent/param.AU), "{:e}".format(self.radius[0][0]))
        # print(r_i_index, r_cent_index, r_cent_index-r_i_index)

        if r_cent_index-r_i_index > 3:
            return True
        else:
            return False
        
    def get_next_m_d(self, h, t):
        '''
            Integrate mass inflow from envelope to disk.

            Parameters:
                h       (lst): index of height of the disk
            Return:
                dm_d    (float): increment of mass of disk from envelope 
        '''

        h = self.radius[t][h]

        if t==0:
            dt = param.time[t]
        if t!=0:
            dt = param.time[t]-param.time[t-1]

        dmdt = 0
        for r in range(len(self.radius)):
            if self.radius[t][r] > self.centrifugal_radius[t]:
                break
            if r > len(h)-1:
                break

            h_index = self._find_closest_index(self.radius[t], h[r])

            radius      = np.sqrt(self.radius[t][r]**2 + self.radius[t][h_index]**2)
            density     = self.cm_cgs_density[t][r][h_index]
            theta       = np.arctan(self.radius[t][h_index]/self.radius[t][r])
            # print("theta", theta)
            # theta       = 0
            velocity    = abs(np.real(self.cm_cgs_v_r[t][r][h_index])*np.sin(theta)) + abs(np.real(self.cm_cgs_v_theta[t][r][h_index]*np.cos(theta)))

            dr = self.r_grid[t][r+1] - self.r_grid[t][r]
            # dz = self.r_grid[t][h_index+1] - self.r_grid[t][h_index]
            # ds = np.sqrt(dr**2 + dz**2)

            dmdt += 4*np.pi * radius * density * velocity * dr

            # L * M/L^3 * L/T * L 
            
            # if r == 0:
            #     dr = self.radius[t][r]
            # else:
            #     dr = (self.radius[t][r]-self.radius[t][r-1])

            # dmdt += radius * density * velocity * dr

        return dmdt*dt

    def volume(self):
        v = 0
        t = 0
        for r in range(1, len(self.r_grid)):
            for h in range(1, len(self.r_grid)):
                dx = (self.r_grid[t][r]-self.r_grid[t][r-1])
                dy = (self.r_grid[t][h]-self.r_grid[t][h-1])
                da = dx * dy     

                v += self.r_grid[t][r] * da

        v *= 2*np.pi

        print(v)
        r_env  = param.G*param.M/(2*param.cs**2)
        print(self.r_grid[0][-1], r_env)
        print(4/3*np.pi * (r_env**3 - self.r_grid[0][0]**3))
        # print(r_env**2-self.r_grid[0][0]**2)
        print(4/3*np.pi * (r_env**3 - self.r_grid[0][0]**3) / (v))

    def _highest_dh(self, past_h, h):
        max_index = 0
        print(len(h), len(past_h), len(self.radius[0]))
        print(h)
        # print(past_h)
        for i in range(len(self.radius[0])):
            # print(i, h[i], past_h[i])
            if h[i]-past_h[i] > max_index:
                max_index = h[i]-past_h[i]

        return max_index
        
    def _check_zero_height(self, lst):
        for i in lst:
            if i != 0:
                return False
        
        return True

    def _initialize_CM_cartesian(self):
        '''
            Initialize and save CM81 data in cartesian coordinate.

            Parameters:
                (None)
            Returns:
                (None)
        '''
        # import CM81 (envelope data)
        data_param = ["CM", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        if os.path.exists(data._get_directory(data_param)):
            print('=== Loading CM81 data... ===')
            radius, r_grid, _cm_nh_density, _cm_cgs_v_r, _cm_cgs_v_theta, cgs_v_phi = data.load(data_param)
            
            print("radius",     len(radius),           len(radius[0]))                                               # [time][radius]
            print("density",    len(_cm_nh_density),   len(_cm_nh_density[0]),    len(_cm_nh_density[0][0]))         # [time][theta][radius]
            print("v_r",        len(_cm_cgs_v_r),      len(_cm_cgs_v_r[0]),       len(_cm_cgs_v_r[0][0]))            # [time][theta][radius]
            print("v_theta",    len(_cm_cgs_v_theta),  len(_cm_cgs_v_theta[0]),   len(_cm_cgs_v_theta[0][0]))        # [time][theta][radius]

        else:
            print("Data does not exist!")
        
        '''
            NOTE: This way of vectorization will sacrifice part of accuracy, as you can't find the 
                the exact same radius after the polar to cartesian coordinate conversion.
                But in general, it should be in the same order of magintude (after searching couple of 
                cases by hand). I assume the accuracy of linear interpolation will be about the same,
                but it takes much longer time.
        '''
        x = radius[0]
        y = radius[0]

        density = np.zeros((len(param.time), len(x), len(y)))
        v_r     = np.zeros((len(param.time), len(x), len(y)))
        v_theta = np.zeros((len(param.time), len(x), len(y)))
        
        # Vectorized index calculation
        x_i, y_i = np.meshgrid(range(len(x)), range(len(y)))
        r_index, theta_index = self._cartesian_index_to_polar(x_i, y_i, radius)

        for t in tqdm(range(len(param.time)), desc='Changing CM81 from cylindrical coordinate to cartesian coordinate...'):
            density[t, :, :] = _cm_nh_density[t, theta_index, r_index] * param.M_H
            v_r[t, :, :] = _cm_cgs_v_r[t, theta_index, r_index]
            v_theta[t, :, :] = _cm_cgs_v_theta[t, theta_index, r_index]

        data_param = ["CM_Cartesian", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        print("Saving Data...")
        data.save(param=data_param, radius=self.radius, r_grid=r_grid, density=density, v_r=v_r, v_theta=v_theta, v_phi=cgs_v_phi)
        print("Successfully Stored Data")

    def _initialize_TSC_cartesian(self):
        '''
            Initialize and save TSC data in cartesian coordinate.

            Parameters:
                (None)
            Returns:
                (None)
        '''
        # import TSC (envelope data)
        data_param = ["TSC", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        if os.path.exists(data._get_directory(data_param)):
            print('=== Loading TSC84 data... ===')
            radius, r_grid, _tsc_nh_density, _tsc_cgs_v_r, _tsc_cgs_v_theta, cgs_v_phi = data.load(data_param)

            print("radius",     len(radius),            len(radius[0]),             len(self.radius[0]))                    # [time][radius]
            print("density",    len(_tsc_nh_density),   len(_tsc_nh_density[0]),    len(_tsc_nh_density[0][0]))         # [time][theta][radius]
            print("v_r",        len(_tsc_cgs_v_r),      len(_tsc_cgs_v_r[0]),       len(_tsc_cgs_v_r[0][0]))            # [time][theta][radius]
            print("v_theta",    len(_tsc_cgs_v_theta),  len(_tsc_cgs_v_theta[0]),   len(_tsc_cgs_v_theta[0][0]))        # [time][theta][radius]
            print("v_phi",      len(cgs_v_phi),         len(cgs_v_phi[0]),          len(cgs_v_phi[0]))                      # [time][theta] 
        else:
            print("Data does not exist!")

        self.radius = [[value * param.AU for value in row] for row in radius]

        # interpolate TSC density and velocity from cylindrical coordinate to cartesian coordinate
        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        # TSC density profile in cartesian coordinate and cgs unit
        self.tsc_cgs_density = []
        self.tsc_cgs_v_r     = []
        self.tsc_cgs_v_theta = []
        for t in tqdm(range(len(param.time)), desc="Changing TSC from cylindrical coordinate to cartesian coordinate..."):
            density = [] 
            v_r     = []
            v_theta = []
            for theta in range(len(theta_values)):
                density.append([i*param.M_H for i in _tsc_nh_density[t][theta]])
                v_r.append(_tsc_cgs_v_r[t][theta])
                v_theta.append(_tsc_cgs_v_r[t][theta])

            # Convert (r, theta) to (x, y) coordinates
            r_values =  np.array(self.radius[t])
            x_values = []
            y_values = []
            for j, theta in enumerate(theta_values):
                x_values.extend(r_values * np.cos(theta))
                y_values.extend(r_values * np.sin(theta))

            # Create a flattened (x, y) array and corresponding flattened data array
            x_flat = np.array(x_values)
            y_flat = np.array(y_values)

            density_flat    = np.array(density).flatten()
            v_r_flat        = np.array(v_r).flatten()
            v_theta_flat    = np.array(v_theta).flatten()

            # Create a grid in (x, y) coordinates
            x_grid = self.radius[t]
            y_grid = self.radius[t]
            xx, yy = np.meshgrid(x_grid, y_grid)

            # Interpolate data onto the (x, y) grid
            density_interp = griddata((x_flat, y_flat), density_flat, (xx, yy), method='linear')
            v_r_interp = griddata((x_flat, y_flat), v_r_flat, (xx, yy), method='linear')
            v_theta_interp = griddata((x_flat, y_flat), v_theta_flat, (xx, yy), method='linear')

            self.tsc_cgs_density.append(density_interp)
            self.tsc_cgs_v_r.append(v_r_interp)
            self.tsc_cgs_v_theta.append(v_theta_interp)

        data_param = ["TSC_Cartesian", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        print("Saving Data...")
        data.save(param=data_param, radius=self.radius, r_grid=r_grid, density=self.tsc_cgs_density, v_r=self.tsc_cgs_v_r, v_theta=self.tsc_cgs_v_theta, cgs_v_phi=cgs_v_phi)
        print("Successfully Stored Data")
    
    def _cartesian_index_to_polar(self, x, y, radius):
        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        l_x = radius[0][x]
        l_y = radius[0][y]

        r     = np.sqrt(l_x**2 + l_y**2)
        theta = np.arctan(l_y/l_x)

        r_index     = self._find_closest_indices(radius[0], r)
        theta_index = self._find_closest_indices(theta_values, theta)

        return r_index, theta_index

    def _extrapolate_ri(self, r_i):
        first_r_i = 0
        for i in r_i:
            if not math.isnan(i):
                first_r_i = i
                break

        for i in range(len(r_i)):
            if math.isnan(r_i[i]):
                r_i[i] = first_r_i
            else:
                break
        
        return r_i

    def _find_closest_index(self, lst, target_value):
        closest_index = min(range(len(lst)), key=lambda i: abs(lst[i] - target_value))

        return closest_index

    def _find_closest_indices(self, lst, target_values):
        # Ensure lst is sorted for binary search
        sorted_indices = np.argsort(lst)
        sorted_lst = lst[sorted_indices]

        # Use binary search to find the closest indices
        closest_indices = np.searchsorted(sorted_lst, target_values, side='left')
        closest_indices = np.clip(closest_indices, 1, len(lst) - 1)  # Avoids index out of bounds

        # Determine the closer index based on absolute differences
        left_diff = np.abs(lst[sorted_indices[closest_indices - 1]] - target_values)
        right_diff = np.abs(lst[sorted_indices[closest_indices]] - target_values)

        closest_indices = np.where(left_diff < right_diff, closest_indices - 1, closest_indices)
        closest_indices = sorted_indices[closest_indices]

        return closest_indices