import numpy as np
from tqdm import tqdm
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

class disk():
    def __init__(self):
        data_param = ["CM", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        if os.path.exists(data._get_directory(data_param)):
            print('Loading CM81 data...')

            self.radius, self.r_grid, self.cm_cgs_density, self.cm_cgs_v_r, self.cm_cgs_v_theta, cgs_v_phi = data.load(data_param)

            print("radius",     len(self.radius),          len(self.radius[0]),          len(self.radius[0]))                    # [time][radius]
            print("density",    len(self.cm_cgs_density),  len(self.cm_cgs_density[0]),   len(self.cm_cgs_density[0][0]))        # [time][theta][radius]
            print("v_r",        len(self.cm_cgs_v_r),      len(self.cm_cgs_v_r[0]),       len(self.cm_cgs_v_r[0][0]))            # [time][theta][radius]
            print("v_theta",    len(self.cm_cgs_v_theta),  len(self.cm_cgs_v_theta[0]),   len(self.cm_cgs_v_theta[0][0]))        # [time][theta][radius]

        else:
            print('CM81 Data Does Not exist!')
            self._initialize_CM_cartesian()

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)
        r = self.radius[0][10]
        theta = theta_values[10]

        print(self.cm_cgs_density[0][0][0] *param.M_H)

        data_param = ["CM_Cartesian", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        if os.path.exists(data._get_directory(data_param)):
            print('Loading CM81 Cartesian data...')

            self.radius, self.r_grid, self.cm_cgs_density, self.cm_cgs_v_r, self.cm_cgs_v_theta, cgs_v_phi = data.load(data_param)

            print("radius",     len(self.radius),          len(self.radius[0]),          len(self.radius[0]))                    # [time][radius]
            print("density",    len(self.cm_cgs_density),  len(self.cm_cgs_density[0]),   len(self.cm_cgs_density[0][0]))        # [time][theta][radius]
            print("v_r",        len(self.cm_cgs_v_r),      len(self.cm_cgs_v_r[0]),       len(self.cm_cgs_v_r[0][0]))            # [time][theta][radius]
            print("v_theta",    len(self.cm_cgs_v_theta),  len(self.cm_cgs_v_theta[0]),   len(self.cm_cgs_v_theta[0][0]))        # [time][theta][radius]

        else:
            print('CM81 Cartesian Data Does Not exist!')
            self._initialize_CM_cartesian()    

        x_i = self._find_closest_index(self.radius[0], r * np.cos(theta))
        y_i = self._find_closest_index(self.radius[0], r * np.sin(theta))

        # print(x_i, y_i)
        r, theta = self._cartesian_index_to_polar(x_i, y_i, self.radius)

        print(self.cm_cgs_density[0][x_i][y_i])
            
    def disk_evol(self):
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
        self.M_disk = [param.min_disk_mass]

        r_star = star.get_r_star()
        star_mass = envelope.cgs_env_to_star_mass(r_star)
        self.r_i = star.disk_inner_radius(star_mass)
        self.centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 for t in param.time]

        first_time = True
        for t in range(len(param.time)):
            # if the disk is smaller than 3 grid points, than don't initiate this progress
            if_start_integration = self.if_start_integration(self.r_i[t], self.centrifugal_radius[t])
            if not if_start_integration:
                continue

            # if it's the first time, than the height of the disk is 0. Otherwise we have to check height.
            if first_time:
                sigma_0 = self.get_sigma_0(m_d=self.M_disk[-1], r_i=self.r_i[t], r_cent=self.centrifugal_radius[t], h=[0])
            if not first_time:
                pass

            sigma_r = [sigma_0 * (self.r_i[t]/r) for r in self.radius[t]]

            if first_time:
                h = [0.1*r for r in self.radius[t]]
                self.check_boundary_w_density(sigma_r=sigma_r, h=h, t=t)

    def get_sigma_0(self, m_d, r_i, r_cent, h):
        '''
            Return sigma_0, the normalization parameter that governs the distribution of surface density (sigma).

            Parameters:
                m_d     (float) : mass of disk for current interation
                r_i     (float) : the inner radius of the disk by evaporation, dependent on luminosity of star
                r_cent  (float) : the centrifugal radius of the disk 
                h       (list)  : the height of the disk, initially we have thin disk assumption H = 0.1R
            Returns:
                sigma_0 (float) : normalization parameter that governs the surface density (sigma) distribution
                                    at R direction
                                    Sigma(R) = sigma_0 (R/R_0)^{-1} 
        '''
        if len(h)==1 and h[0] == 0:
            sigma_0 = m_d / (4*np.pi * r_i * (r_cent-r_i))
            return sigma_0

    def check_boundary_w_density(self, sigma_r, h, t):
        '''
            Check the boundary of the disk with density method. This is by identifying at which height does the density structure
            of the disk match with the density of the envelope.

            Parameters:
                sigma_r (list)  : the sigma at the midplane throughout the radius
                h       (list)  : the height of the disk, initially we have thin disk assumption H = 0.1R
                t       (int)   : the index of param.time specifying when we are dealing the issue
            Returns:
        '''

        height = []
        found = False
        for r in range(len(self.radius[t])):
            rho_d_midplane = sigma_r[r] / (h[r]*np.sqrt(2*np.pi))
            if rho_d_midplane < self.cm_cgs_density[t][r][0]:               # here we assume that it is in [time][x][y] format
                continue
            for z in range(len(self.radius[t])):
                rho_d_z_prev = rho_d_midplane * np.exp(-(self.radius[z-1])**2/(2*h**2))
                rho_d_z = rho_d_midplane * np.exp(-self.radius[z]**2/(2*h**2))
                if rho_d_z < self.cm_cgs_density[t][r][z] and rho_d_z_prev > self.cm_cgs_density[t][r][z]:
                    height.append(h)
                    found = True
                    break
            if not found:
                height.append(0)

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

        # print(r_i_index, r_cent_index, r_cent_index-r_i_index)

        if r_cent_index-r_i_index > 3:
            return True
        else:
            return False
        
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
        self.radius = [[value * param.AU for value in row] for row in radius]
        data.save(param=data_param, radius=self.radius, r_grid=r_grid, density=density, v_r=v_r, v_theta=v_theta, v_phi=cgs_v_phi)
        print("Successfully Stored Data")

        # # import CM81 (envelope data)
        # data_param = ["CM", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        # if os.path.exists(data._get_directory(data_param)):
        #     print('=== Loading CM81 data... ===')
        #     radius, r_grid, _cm_nh_density, _cm_cgs_v_r, _cm_cgs_v_theta, cgs_v_phi = data.load(data_param)
            
        #     print("radius",     len(radius),           len(radius[0]),            len(radius[0]))                    # [time][radius]
        #     print("density",    len(_cm_nh_density),   len(_cm_nh_density[0]),    len(_cm_nh_density[0][0]))         # [time][theta][radius]
        #     print("v_r",        len(_cm_cgs_v_r),      len(_cm_cgs_v_r[0]),       len(_cm_cgs_v_r[0][0]))            # [time][theta][radius]
        #     print("v_theta",    len(_cm_cgs_v_theta),  len(_cm_cgs_v_theta[0]),   len(_cm_cgs_v_theta[0][0]))        # [time][theta][radius]

        # else:
        #     print("Data does not exist!")

        # self.radius = [[value * param.AU for value in row] for row in radius]

        # # interpolate CM81 density and velocity from cylindrical coordinate to cartesian coordinate
        # theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        # # CM81 density profile in cartesian coordinate and cgs unit
        # self.cm_cgs_density = []
        # self.cm_cgs_v_r     = []
        # self.cm_cgs_v_theta = []
        # for t in tqdm(range(len(param.time)), desc="Changing CM81 from cylindrical coordinate to cartesian coordinate..."):
        #     density = [] 
        #     v_r     = []
        #     v_theta = []
        #     for theta in range(len(theta_values)):
        #         density.append([i*param.M_H for i in _cm_nh_density[t][theta]])
        #         v_r.append(_cm_cgs_v_r[t][theta])
        #         v_theta.append(_cm_cgs_v_r[t][theta])

        #     # Convert (r, theta) to (x, y) coordinates
        #     r_values =  np.array(self.radius[t])
        #     x_values = []
        #     y_values = []
        #     for j, theta in enumerate(theta_values):
        #         x_values.extend(r_values * np.cos(theta))
        #         y_values.extend(r_values * np.sin(theta))

        #     # Create a flattened (x, y) array and corresponding flattened data array
        #     x_flat = np.array(x_values)
        #     y_flat = np.array(y_values)

        #     density_flat    = np.array(density).flatten()
        #     v_r_flat        = np.array(v_r).flatten()
        #     v_theta_flat    = np.array(v_theta).flatten()

        #     # Create a grid in (x, y) coordinates
        #     x_grid = self.radius[t]
        #     y_grid = self.radius[t]
        #     xx, yy = np.meshgrid(x_grid, y_grid)

        #     # Interpolate data onto the (x, y) grid
        #     density_interp = griddata((x_flat, y_flat), density_flat, (xx, yy), method='linear')
        #     v_r_interp     = griddata((x_flat, y_flat), v_r_flat, (xx, yy), method='linear')
        #     v_theta_interp = griddata((x_flat, y_flat), v_theta_flat, (xx, yy), method='linear')

        #     self.cm_cgs_density.append(density_interp)
        #     self.cm_cgs_v_r.append(v_r_interp)
        #     self.cm_cgs_v_theta.append(v_theta_interp)

        # data_param = ["CM_Cartesian", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        # print("Saving Data...")
        # data.save(param=data_param, radius=self.radius, r_grid=r_grid, density=self.cm_cgs_density, v_r=self.cm_cgs_v_r, v_theta=self.cm_cgs_v_theta, v_phi=cgs_v_phi)
        # print("Successfully Stored Data")

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

        # print(r * np.sin(theta))
        # print(l_y)

        r_index     = self._find_closest_indices(radius[0], r)
        theta_index = self._find_closest_indices(theta_values, theta)

        # print(r)
        # print(radius[0][r_index])

        # print(arr)
        # print(radius[0][x])

        return r_index, theta_index


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

