import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata, interp2d, Rbf, interp1d
from scipy.integrate import quad
import numpy as np
from tqdm import tqdm
import os

from envelope.shu77 import shu
from envelope.tsc84 import tsc
from envelope.cm81 import cm81_envelope
from data_handler import data
from parameter import param

param = param()
data = data()

# TODO: make grid first -> change to CM81 and TSC84 grid (z=1/x, r=axt)
# TODO: Disk mass (done) -> Sigma -> disk boundary -> disk mass (...cycle)
# TODO: Change 1d disk mass integration to 2D

#   Envelope
class shu_graph():
    def __init__(self, z_s1, z_g1, r_grid, save=False):
        """
        Initialize Shu (1977), using Fluid model to model the 
        gravitational collapse of isothermal cloud with no rotation.

        Args:
            save (bool): Whether to save the graph. Default is False.

        (Output):
            self.radius     [cm]/AU
            self.nh_density [nh] -> *param.M_H to go back to cgs
            self.v_r        [cgs]
        """
        self.z_s1 = z_s1
        self.z_g1 = z_g1
        self.r_grid = r_grid
        self.save = save

        self.data_param = ["Shu", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]

        if os.path.exists(data._get_directory(self.data_param)):
            print('Loading data...')
            self.radius, self.r_grid, self.nh_density, self.v_r, self.v_theta, self.v_phi = data.load(self.data_param)

        else:
            # Create a Shu instance with specified parameters
            self.shu = shu(z_g1=z_g1, z_s1=z_s1, cs=param.cs, time=param.time)

            # Initialize self-similar variables
            self.a_shu_s1   = self.shu.get_a_s1()
            self.v_shu_s1   = self.shu.get_v_s1()
            self.a_shu_g1   = self.shu.get_a_g1()
            self.v_shu_g1   = self.shu.get_v_g1()
            self.dadz       = self.shu.get_dadz()
            self.dvdz       = self.shu.get_dvdz()

            # Initialize dimension variables
            self.radius_g1      = self.shu.get_r_g1()
            self.radius_s1      = self.shu.get_r_s1()
            self.cgs_density_g1 = self.shu.get_cgs_rho_g1()
            self.cgs_density_s1 = self.shu.get_cgs_rho_s1()
            self.nh_density_g1  = self.shu.get_nh_rho_g1()
            self.nh_density_s1  = self.shu.get_nh_rho_s1()

            self.radius     = []
            self.nh_density = []
            self.v_r        = []
            for i in range(len(param.time)):
                self.radius.append([j/param.AU for j in self.radius_s1[i]][::-1]+[j/param.AU for j in self.radius_g1[i]][::-1])
                self.nh_density.append(self.nh_density_s1[i][::-1]+self.nh_density_g1[i][::-1])
                self.v_r.append([param.cs*j for j in self.v_shu_s1][::-1]+[param.cs*j for j in self.v_shu_g1][::-1])

            if self.save == True:
                data.save(param=self.data_param, radius=self.radius, r_grid=self.r_grid, density=self.nh_density, v_r=self.v_r)
                print('Successfully stored data')


    def dimensionless_a_v_1d(self):
        """
        Plot the dimensionless density and velocity profiles in 1D,
        against the radius (1/z).
        """
        if os.path.exists(data._get_directory(self.data_param)):
            self.shu = shu(z_g1=self.z_g1, z_s1=self.z_s1, cs=param.cs, time=param.time)

        self.a_shu_s1   = self.shu.get_a_s1()
        self.v_shu_s1   = self.shu.get_v_s1()
        self.a_shu_g1   = self.shu.get_a_g1()
        self.v_shu_g1   = self.shu.get_v_g1()

        plt.plot([1/i for i in self.z_g1], self.a_shu_g1)
        plt.plot([1/i for i in self.z_s1], self.a_shu_s1)
        plt.plot([1/i for i in self.z_g1], [-i for i in self.v_shu_g1])
        plt.plot([1/i for i in self.z_s1], [-i for i in self.v_shu_s1])

    def cgs_density_1d(self):
        """
        Plot the density profiles in CGS units in 1D for different time.
        """
        for i in range(len(param.time)):
            plt.plot(self.radius[i], [j*param.M_H for j in self.nh_density[i]], color=param.color_1d[i], label=f't={param.time[i]/1e12}e12')

    def nh_density_1d(self):
        """
        Plot the density profiles in density of hydrogen in 1D for different time.
        """
        for i in range(len(param.time)):
            plt.plot(self.radius[i], self.nh_density[i], color=param.color_1d[i], label=f't={param.time[i]/1e12}e12')

    def cgs_density_2d(self):
        """
        Plot density profiles in CGS units in 2D for different times using 
        a logarithmic colormap with contour value. Note that it shows the entire range we calculated

        LOG SCALE!
        """
        fig, axes = plt.subplots(2, 4)
        axes_flat = axes.flatten()

        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                
                radii = np.linspace(0, 2*param.sample_pts, 2*param.sample_pts)  # Distances from the origin
                molecule_counts = [j*param.M_H for j in self.nh_density[i]]

                # Create a grid
                x = np.linspace(0, 2*param.sample_pts, 2*param.sample_pts)
                y = np.linspace(0, 2*param.sample_pts, 2*param.sample_pts)
                xx, yy = np.meshgrid(x, y)
                r = np.sqrt(xx**2 + yy**2)

                # Interpolate data onto the grid
                interpolated_data = np.interp(r.ravel(), radii, molecule_counts, left=0, right=0)
                interpolated_data = interpolated_data.reshape(xx.shape)

                # Logarithmic colormap
                norm = LogNorm(vmin=param.cgs_vmin, vmax=param.cgs_vmax)

                # Plot the radially symmetrical image
                im = ax.imshow(interpolated_data, cmap='gist_rainbow_r', origin='lower', 
                               extent=[0, 2*param.sample_pts, 0, 2*param.sample_pts], norm=norm)

                # Add contour line at the value 10^(-18)
                contour_value = param.cgs_contour 
                ax.contour(xx, yy, interpolated_data, levels=[contour_value], colors='red', linewidths=2)

                ax.set_title(f't = {param.time[i]/param.t_acc}')

            else:
                plt.axis("off")
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_label('Density')

        fig.suptitle('Shu (1977) 2D Graph (cgs)')

    def nh_density_2d(self):
        """
        Plot density profiles in density of hydrogen in 2D for different times using 
        a logarithmic colormap with contour value. Note that it only show the range in Visser.
        """

        fig, axes = plt.subplots(2, 4)
        axes_flat = axes.flatten()

        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                radii = [j/1e3 for j in self.radius[i]]  # Distances from the origin
                molecule_counts = self.nh_density[i]
                
                # Create a grid
                x = [j/1e3 for j in self.radius[i]] 
                y = [j/1e3 for j in self.radius[i]] 

                xx, yy = np.meshgrid(x, y)
                r = np.sqrt(xx**2 + yy**2)

                # Interpolate data onto the grid
                interpolated_data = np.interp(r.ravel(), radii, molecule_counts, left=0, right=0)
                interpolated_data = interpolated_data.reshape(xx.shape)

                # Logarithmic colormap
                norm = LogNorm(vmin=param.visser_vmin, vmax=param.visser_vmax)

                # Plot the radially symmetrical image
                im = ax.contourf(x, y, interpolated_data, cmap='gist_rainbow_r', levels=np.logspace(4,9,100), norm=norm, aspect='auto')
                
                ax.set_xlim(0, 1.5)
                ax.set_ylim(0, 1.5)

                r_env_edge = patches.Circle((0,0), param.r_env/1e3, fill=False, edgecolor='black', linewidth=2)
                ax.add_patch(r_env_edge)

                # Add contour line at the value
                ax.contour(xx, yy, interpolated_data, levels=[param.visser_contour],         
                           colors='red', linewidths=2)

                ax.set_title(f't = {param.time[i]/param.t_acc}')

            else:
                plt.axis("off")
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_ticks([10**4, 10**5, 10**6, 10**7, 10**8, 10**9])
                cbar.set_label('Density')

        fig.suptitle('Shu (1977) 2D Graph (n_h)')

class tsc_graph():
    def __init__(self, z_s1, z_g1, r_grid, r, save=False):
        """
        Initialize Terebey, Shu & Cassen (1984), using Fluid model to model the 
        gravitational collapse of isothermal cloud with small rotation.

        Args:
            save (bool): Whether to save the graph. Default is False.

        (Output):
            self.radius     [cm]/AU
            self.nh_density [nh] -> *param.M_H to go back to cgs
            self.v_r        [cgs]
        """
        self.z_s1 = z_s1
        self.z_g1 = z_g1
        self.r_grid = r_grid
        self.save = save

        self.tau_2  = [(param.omega*t)**2* param.cs*t /param.AU for t in param.time]
        self.r_c    = [param.cs/param.omega /param.AU for t in param.time]
        self.r_d    = [(param.omega**2*param.G**3*param.M**3) / (16*param.cs**8) / param.AU for t in param.time]
        self.r_env  = [param.G*param.M/(2*param.cs**2) / param.AU for t in param.time]

        self.centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 / param.AU for t in param.time]

        self.data_param = ["TSC", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]

        if os.path.exists(data._get_directory(self.data_param)):
            print('Loading TSC84 data...')
            self.radius, self.r_grid, self.nh_density, self.cgs_v_r, self.cgs_v_theta, self.cgs_v_phi = data.load(self.data_param)

            print("radius",     len(self.radius),       len(self.radius[0]),        len(self.radius[0]))                # [time][radius]
            print("density",    len(self.nh_density),   len(self.nh_density[0]),    len(self.nh_density[0][0]))         # [time][theta][radius]
            print("v_r",        len(self.cgs_v_r),      len(self.cgs_v_r[0]),       len(self.cgs_v_r[0][0]))            # [time][theta][radius]
            print("v_theta",    len(self.cgs_v_theta),  len(self.cgs_v_theta[0]),   len(self.cgs_v_theta[0][0]))        # [time][theta][radius]
            print("v_phi",      len(self.cgs_v_phi),    len(self.cgs_v_phi[0]),     len(self.cgs_v_phi[0]))             # [time][theta] 

        else:
            # Create Shu and TSC instances with specified parameters
            self.shu = shu(z_g1=self.z_g1, z_s1=self.z_s1, cs=param.cs, time=param.time)
            self.tsc = tsc(z_g1=self.z_g1, z_s1=self.z_s1, a_shu_s1=self.shu.get_a_s1(), a_shu_g1=self.shu.get_a_g1(), 
                            v_shu_s1=self.shu.get_v_s1(), v_shu_g1=self.shu.get_v_g1(),
                            store_dadz=self.shu.get_dadz(), store_dvdz=self.shu.get_dvdz(),
                            cs=param.cs, time=param.time)

            # Initialize dimension variables
            print("Initializing...")
            print("Initializing radius")
            self.radius_s1      = self.tsc.get_r_s1()
            self.radius_g1      = self.tsc.get_r_g1()

            print("Interpolating Radius...")


            print("Initializing density")
            # self.cgs_density_s1 = self.tsc.get_cgs_density_s1()
            # self.cgs_density_g1 = self.tsc.get_cgs_density_g1()
            self.nh_density_s1  = self.tsc.get_nh_density_s1().tolist()
            self.nh_density_g1  = self.tsc.get_nh_density_g1().tolist()

            print(len(self.nh_density_s1), len(self.nh_density_s1[0]), len(self.nh_density_s1[0][0]))

            print("Initializing cgs velocity")
            self.cgs_v_r_s1     = self.tsc.get_cgs_v_r_s1().tolist()
            self.cgs_v_r_g1     = self.tsc.get_cgs_v_r_g1().tolist()
            self.cgs_v_theta_s1 = self.tsc.get_cgs_v_theta_s1().tolist()
            self.cgs_v_theta_g1 = self.tsc.get_cgs_v_theta_g1().tolist()
            self.cgs_v_phi      = self.tsc.get_cgs_v_phi()

            # print("Initializing a, v, w")
            # self.quad_a_g1      = self.tsc.get_quad_a_g1()
            # self.quad_a_s1      = self.tsc.get_quad_a_s1()
            # self.quad_v_g1      = self.tsc.get_quad_v_g1()
            # self.quad_v_s1      = self.tsc.get_quad_v_s1()
            # self.quad_w_g1      = self.tsc.get_quad_w_g1()
            # self.quad_w_s1      = self.tsc.get_quad_w_s1()
            # self.mono_a_g1      = self.tsc.get_mono_a_g1()
            # self.mono_a_s1      = self.tsc.get_mono_a_s1()
            # self.mono_v_s1      = self.tsc.get_mono_v_s1()

            # print(param.omega*param.M*param.G/0.975/param.cs**3)
            # print([param.omega*t for t in param.time])

            self.radius         = np.zeros((len(param.time), len(self.radius_s1[0]+self.radius_g1[0])))
            self.nh_density     = np.zeros((len(param.time), param.theta_pts, len(self.nh_density_s1[0][0]+self.nh_density_g1[0][0])))
            self.cgs_v_r        = np.zeros((len(param.time), param.theta_pts, len(self.cgs_v_r_s1[0][0]+self.cgs_v_r_g1[0][0])))
            self.cgs_v_theta    = np.zeros((len(param.time), param.theta_pts, len(self.cgs_v_theta_s1[0][0]+self.cgs_v_theta_s1[0][0])))
            self.cgs_v_phi      = self.cgs_v_phi
            
            for i in tqdm(range(len(param.time)), desc='reshape process'):
                arr_radius = [j/param.AU for j in self.radius_s1[i]][::-1]+[j/param.AU for j in self.radius_g1[i]][::-1]
                for j in range(len(arr_radius)):
                    self.radius[i][j] = arr_radius[j]

                for j in range(param.theta_pts):
                    arr_density = self.nh_density_s1[i][j][::-1]+self.nh_density_g1[i][j][::-1]
                    arr_v_r     = self.cgs_v_r_s1[i][j][::-1]+self.cgs_v_r_g1[i][j][::-1]
                    arr_v_theta = self.cgs_v_theta_s1[i][j][::-1]+self.cgs_v_theta_g1[i][j][::-1]
                    for k in range(len(arr_density)):
                        self.nh_density[i][j][k]    = arr_density[k]
                        self.cgs_v_r[i][j][k]       = arr_v_r[k]
                        self.cgs_v_theta[i][j][k]   = arr_v_theta[k]
            
            print("radius", self.radius.shape)                                                              # [time][radius]
            print("density", self.nh_density.shape)                                                         # [time][theta][radius]
            print("v_r", self.cgs_v_r.shape)                                                                # [time][theta][radius]
            print("v_theta", self.cgs_v_theta.shape)                                                        # [time][theta][radius]
            print(f"v_phi, ({len(self.cgs_v_phi)}, {len(self.cgs_v_phi[0])}, {len(self.cgs_v_phi[0])})")     # [time][theta]

            theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)
            
            self.interp_density = []
            self.interp_v_r     = []
            self.interp_v_theta = []
            for t in tqdm(range(len(param.time)), desc='Interpolation'):
                # original radius 
                r_values =  np.array([i*param.AU for i in self.radius[t]])      
                # interpolated radius is r[t]
                print(r_values[0], r_values[-1])
                print(r[t][0], r[t][-1])

                # 1d extrapolation at the theta direction
                self.interp_density.append([])
                self.interp_v_r.append([])
                self.interp_v_theta.append([])
                for theta in range(len(theta_values)):
                    f_1d_density_interp = interp1d(r_values, self.nh_density[t][theta], fill_value='extrapolate')
                    f_1d_v_r_interp = interp1d(r_values, self.cgs_v_r[t][theta], fill_value='extrapolate')
                    f_1d_v_theta_interp = interp1d(r_values, self.cgs_v_theta[t][theta], fill_value='extrapolate')

                    self.interp_density[-1].append(f_1d_density_interp(r[t]))
                    self.interp_v_r[-1].append(f_1d_v_r_interp(r[t]))
                    self.interp_v_theta[-1].append(f_1d_v_theta_interp(r[t]))

                plt.plot(r_values, self.nh_density[t][100], color='blue')
                print(self.interp_density[t][100])
                plt.plot(r[t], [abs(i) for i in self.interp_density[t][100]], color='red')

                # Reshape 1d interpolation into 2d
                




                # # Original radius
                # r_values =  np.array(self.radius[t])
                
                # x_values = []
                # y_values = []
                # for j, theta in enumerate(theta_values):
                #     x_values.extend(r_values * np.cos(theta))
                #     y_values.extend(r_values * np.sin(theta))
                # # Create a flattened (x, y) array and corresponding flattened data array
                # x_flat = np.array(x_values)
                # y_flat = np.array(y_values)

                # # (x, y) grid we wish to interpolate on
                # x_grid = r[t]
                # y_grid = r[t]
                # xx, yy = np.meshgrid(x_grid, y_grid)

                # print("Org value: x", x_flat[0], x_flat[-1], "y", y_flat[0], y_flat[-1])
                # print("Interp   : x", x_grid[0], x_grid[-1], "y", y_grid[0], y_grid[-1])

                # # data
                # density_list = [] 
                # for theta in range(len(theta_values)):
                #     density_list.append(self.nh_density[t][theta])
                # density_flat = np.array(density_list).flatten()

                # v_r_list = []
                # for theta in range(len(theta_values)):
                #     v_r_list.append(self.cgs_v_r[t][theta])
                # v_r_flat = np.array(v_r_list).flatten()

                # v_theta_list = []
                # for theta in range(len(theta_values)):
                #     v_theta_list.append(self.cgs_v_theta[t][theta])
                # v_theta_flat = np.array(v_theta_list).flatten()

                # rf_density = Rbf(x_flat, y_flat, density_flat, epsilon=2)
                # rf_v_r     = Rbf(x_flat, y_flat, v_r_flat, epsilon=2)
                # rf_v_theta = Rbf(x_flat, y_flat, v_theta_flat, epsilon=2)

                # density_interp = rf_density(xx, yy)
                # v_r_interp     = rf_v_r(xx, yy)
                # v_theta_interp = rf_v_theta(xx, yy)        

                # # print(density_interp.shape, v_r_interp.shape, v_theta_interp.shape)

                # self.interp_density.append(density_interp)
                # self.interp_v_r.append(v_r_interp)
                # self.interp_v_theta.append(v_theta_interp)

                # density_interp = griddata((x_flat, y_flat), density_flat, (xx, yy), method='linear')
                # v_r_interp = griddata((x_flat, y_flat), v_r_flat, (xx, yy), method='linear')
                # v_theta_interp = griddata((x_flat, y_flat), v_theta_flat, (xx, yy), method='linear')

            if self.save == True:
                print("Saving Data...")
                data.save(param=self.data_param, radius=r, r_grid=r_grid, density=self.interp_density, 
                          v_r=self.interp_v_r, v_theta=self.interp_v_theta, v_phi=self.cgs_v_phi)
                print('Successfully stored data')
       
    def dimensionless_a_1d(self):
        """
        Plot the dimensionless density profiles in 1D against radius (1/z).
        """
        if os.path.exists(data._get_directory(self.data_param)):
            self.shu = shu(z_g1=self.z_g1, z_s1=self.z_s1, cs=param.cs, time=param.time)
            self.tsc = tsc(z_g1=self.z_g1, z_s1=self.z_s1, a_shu_s1=self.shu.get_a_s1(), a_shu_g1=self.shu.get_a_g1(), 
                            v_shu_s1=self.shu.get_v_s1(), v_shu_g1=self.shu.get_v_g1(),
                            store_dadz=self.shu.get_dadz(), store_dvdz=self.shu.get_dvdz(),
                            cs=param.cs, time=param.time)
        
        # Plot Shu acceleration profiles
        plt.plot([1/i for i in self.z_g1], [i for i in self.shu.get_a_g1()], color='orange', label='a shu', linewidth=2)
        plt.plot([1/i for i in self.z_s1], [i for i in self.shu.get_a_s1()], color='orange', linewidth=2)
        
        # Plot quadrupolar acceleration profiles
        plt.plot([1/i for i in self.z_g1], [i for i in self.tsc.get_quad_a_g1()], label=r"$-\alpha_Q$", color='red', linewidth=2)
        plt.plot([1/i for i in self.z_s1], [abs(i) for i in self.tsc.get_quad_a_s1()], color='red', linewidth=2)                 
        
        # Plot monopolar acceleration profiles
        plt.plot([1/i for i in self.z_g1], [i for i in self.tsc.get_mono_a_g1()], label=r"$a_m$", color='blue', linewidth=2)
        plt.plot([1/i for i in self.z_s1], [i for i in self.tsc.get_mono_a_s1()], color='blue', linewidth=2)
    
    def dimensionless_quad_v_w_1d(self):
        """
        Plot dimensionless quadrupolar velocity in theta and phi 
        direction profiles in 1D against radius (1/z).
        """
        if os.path.exists(data._get_directory(self.data_param)):
            self.shu = shu(z_g1=self.z_g1, z_s1=self.z_s1, cs=param.cs, time=param.time)
            self.tsc = tsc(z_g1=self.z_g1, z_s1=self.z_s1, a_shu_s1=self.shu.get_a_s1(), a_shu_g1=self.shu.get_a_g1(), 
                            v_shu_s1=self.shu.get_v_s1(), v_shu_g1=self.shu.get_v_g1(),
                            store_dadz=self.shu.get_dadz(), store_dvdz=self.shu.get_dvdz(),
                            cs=param.cs, time=param.time)
        
        # Plot quadrupolar velocity profiles
        plt.plot([1/i for i in self.z_g1], [i for i in self.tsc.get_quad_v_g1()], label=r"$v_Q$", color='blue', linewidth=2)
        plt.plot([1/i for i in self.z_s1], [i for i in  self.tsc.get_quad_v_s1()], color='blue', linewidth=2)       
        
        # Plot quadrupolar angular velocity profiles
        plt.plot([1/i for i in self.z_g1], [i for i in self.tsc.get_quad_w_g1()], label=r"$-w_Q$", color='green', linewidth=2)
        plt.plot([1/i for i in self.z_s1], [-i for i in self.tsc.get_quad_w_s1()], color='green', linewidth=2)

    def dimensionless_mono_v_1d(self):
        """
        Plot dimensionless monopolar velocity profiles in 1D against radius (1/z).
        """
        if os.path.exists(data._get_directory(self.data_param)):
            self.shu = shu(z_g1=self.z_g1, z_s1=self.z_s1, cs=param.cs, time=param.time)
            self.tsc = tsc(z_g1=self.z_g1, z_s1=self.z_s1, a_shu_s1=self.shu.get_a_s1(), a_shu_g1=self.shu.get_a_g1(), 
                            v_shu_s1=self.shu.get_v_s1(), v_shu_g1=self.shu.get_v_g1(),
                            store_dadz=self.shu.get_dadz(), store_dvdz=self.shu.get_dvdz(),
                            cs=param.cs, time=param.time)
        
        # Plot Shu and monopolar velocity profiles
        plt.plot([1/i for i in self.z_s1], [-i for i in self.shu.get_v_s1()], color='orange', label='-v shu', linewidth=2)
        plt.plot([1/i for i in self.z_s1], [i for i in self.tsc.get_mono_v_s1()], label=r"$v_m$", linewidth=2)

    def cgs_density_1d(self):
        """
        Plot density profiles in CGS units in 1D for different times.
        """
        for i in range(len(param.time)):
            print(self.nh_density[i][100])
            plt.plot(self.radius[i], [j*param.M_H for j in self.nh_density[i][0]], label=f't={param.time[i]/1e12}e12')
    
    def nh_density_1d(self):
        """
        Plot density profiles in density of hydrogen in 1D for different times.
        """
        for i in range(len(param.time)):
            plt.plot(self.radius[i], self.nh_density[i][0], label=f't={param.time[i]/1e12}e12')

    def cgs_density_2d(self):
        """
        Plot density profiles in CGS units in 2D for different times using 
        a logarithmic colormap with contour value.
        """
        fig, axes = plt.subplots(2, 2)
        axes_flat = axes.flatten()

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        i = 0
        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                data_values = [] 
                for theta in range(len(theta_values)):
                    data_values.append([j*param.M_H for j in self.nh_density[i][theta]])

                # Example data in (r, theta) format
                r_values = np.linspace(0, 2*param.sample_pts, 2*param.sample_pts)
                # data_values = [i*param.M_H for i in self.nh_density]  # List of lists, each sublist represents data at a theta value

                # Convert (r, theta) to (x, y) coordinates
                x_values = []
                y_values = []
                for j, theta in enumerate(theta_values):
                    x_values.extend(r_values * np.cos(theta))
                    y_values.extend(r_values * np.sin(theta))

                # Create a flattened (x, y) array and corresponding flattened data array
                x_flat = np.array(x_values)
                y_flat = np.array(y_values)
                data_flat = np.array(data_values).flatten()

                # Create a grid in (x, y) coordinates
                x_grid = np.linspace(0, max(x_flat), 2*param.sample_pts)
                y_grid = np.linspace(0, max(y_flat), 2*param.sample_pts)
                xx, yy = np.meshgrid(x_grid, y_grid)

                # Interpolate data onto the (x, y) grid
                data_interp = griddata((x_flat, y_flat), data_flat, (xx, yy), method='linear')

                # Plot the interpolated data using imshow()
                norm = LogNorm(vmin=param.cgs_vmin, vmax=param.cgs_vmax)
                im = ax.imshow(data_interp, extent=[0, max(x_flat), 0, max(y_flat)],
                        cmap='gist_rainbow_r', origin='lower', aspect='auto', norm=norm)
                
                contour_value = param.cgs_contour
                ax.contour(xx, yy, data_interp, levels=[contour_value], colors='red', linewidths=2)

                ax.set_title(f't = {param.time[i]/param.t_acc}')

            else:
                plt.axis('off')
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_label('Density')

        fig.suptitle('TSC (1984) 2D Graph (cgs)')

    def nh_density_2d(self):
        """
        Plot density profiles in density of hydrogen units in 2D for different times using 
        contour plots with a logarithmic colormap.
        """

        fig, axes = plt.subplots(2, 2)
        axes_flat = axes.flatten()

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        i = 0
        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"graphing progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                
                molecule_counts = [] 
                for theta in range(len(theta_values)):
                    molecule_counts.append(self.nh_density[i][theta])

                # Example data in (r, theta) format
                r_values =  np.array([j/1e3 for j in self.radius[i]])
                
                # Convert (r, theta) to (x, y) coordinates
                x_values = []
                y_values = []
                for j, theta in enumerate(theta_values):
                    x_values.extend(r_values * np.cos(theta))
                    y_values.extend(r_values * np.sin(theta))

                # Create a flattened (x, y) array and corresponding flattened data array
                x_flat = np.array(x_values)
                y_flat = np.array(y_values)
                data_flat = np.array(molecule_counts).flatten()

                # Create a grid in (x, y) coordinates
                x_grid = [j/1e3 for j in self.radius[i]]
                y_grid = [j/1e3 for j in self.radius[i]]
                xx, yy = np.meshgrid(x_grid, y_grid)

                print("Org value: x", x_flat[0], x_flat[-1], "y", y_flat[0], y_flat[-1])
                print("Interp   : x", x_grid[0], x_grid[-1], "y", y_grid[0], y_grid[-1])

                # Interpolate data onto the (x, y) grid
                data_interp = griddata((x_flat, y_flat), data_flat, (xx, yy), method='linear')

                # Create a contour plot with a logarithmic colormap
                norm = LogNorm(vmin=param.visser_vmin, vmax=param.visser_vmax)
                im = ax.contourf(xx, yy, data_interp, cmap='gist_rainbow_r', levels=np.logspace(4,9,100), norm=norm)

                ax.set_xlim(0, 10)
                ax.set_ylim(0, 10)

                # Add contour line at the value
                # ax.contour(xx, yy, data_interp, levels=[param.visser_contour],           
                #            colors='red', linewidths=2)
                
                r_env_edge = patches.Circle((0,0), self.r_env[i]/1e3, fill=False, edgecolor='white', linewidth=2)
                ax.add_patch(r_env_edge)
                circ_tau_2 = patches.Circle((0,0), self.tau_2[i]/1e3, fill=False, edgecolor='black', linewidth=2)
                ax.add_patch(circ_tau_2)
                circ_r_c = patches.Circle((0,0), self.r_c[i]/1e3, fill=False, edgecolor='red', linewidth=2)
                ax.add_patch(circ_r_c)
                circ_r_d = patches.Circle((0,0), self.r_d[i]/1e3, fill=False, edgecolor='blue', linewidth=2)
                ax.add_patch(circ_r_d)
                circ_centrifugal_radius = patches.Circle((0,0), self.centrifugal_radius[i]/1e3, fill=False, edgecolor='green', linewidth=2)
                ax.add_patch(circ_centrifugal_radius)
                circ_5_AU = patches.Circle((0,0), 5/1e3, fill=False, edgecolor='yellow', linewidth=2)
                ax.add_patch(circ_5_AU)

                ax.set_title(f't = {param.time[i]/param.t_acc}')
                
            else:
                plt.axis('off')
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_ticks([10**4, 10**5, 10**6, 10**7, 10**8, 10**9])
                cbar.set_label('Density')
    
        fig.suptitle('TSC (1984) 2D Graph (n_h)')
    
    def int_envelope_mass(self):
        print("Integrating...")

        # times AU for all self.radius
        radius_AU = [[value * param.AU for value in row] for row in self.radius]

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)
        
        # 2 * 2*pi S(rho * r^2 sin(theta) )dr  d(theta)
        # d(theta) = pi/2 /theta pts
        self.mass = []
        v = []
        for t in range(len(param.time)):
            print(f"Integrating progress: {t+1}/{len(param.time)}")
            mass_t = 0 
            v_t = 0
            for theta in range(param.theta_pts):
                for i in range(param.sample_pts):
                    # unit test for radius and r_grid
                    if radius_AU[t][i]-self.r_grid[t][i] < 0 or self.r_grid[t][i+1]-radius_AU[t][i]<0:
                        print(self.r_grid[t][i], radius_AU[t][i], self.r_grid[t][i+1])                        
 
                    if radius_AU[t][i+1] > self.tau_2[t]*param.AU and radius_AU[t][i+1] < self.r_env[t]*param.AU:
                        if param.M_H*self.nh_density[t][theta][i]<0:
                            print(t, theta, i, param.M_H*self.nh_density[t][theta][i])
                        
                        if param.M_H*self.nh_density[t][theta][i]>0:
                            mass_t  += param.M_H*self.nh_density[t][theta][i] * radius_AU[t][i]**2 * np.sin(theta_values[theta]) * (self.r_grid[t][i+1]-self.r_grid[t][i])
                            v_t     += radius_AU[t][i]**2 * np.sin(theta_values[theta]) * (self.r_grid[t][i+1]-self.r_grid[t][i])
           
            self.mass.append(4*np.pi*mass_t * np.pi/2 / param.theta_pts)
            v.append(4*np.pi*v_t * np.pi/2 / param.theta_pts)
        
        print("Volume:", [4/3*np.pi * (self.r_env[t]*param.AU)**3 for t in range(len(v))], "Numerical:", v)
        print("Ratio of volume", [v[t] / (4/3*np.pi * (self.r_env[t]*param.AU)**3) for t in range(len(v))])
        print("Mass:", param.M, self.mass)
        print("Adjusted Mass", [self.mass[t] / (v[t] / (4/3*np.pi * (self.r_env[t]*param.AU)**3)) for t in range(len(v)-1)])
    
    def cgs_mass_env_to_star_disk(self, time):
        r_star = []
        for t in [i/31556926 for i in time]:
            if t>20000 and t<20100:
                r_star.append(5 * param.AU * (1-np.sqrt(t-20,000)/100)) + param.r_ps
            else:
                r_star.append(param.r_ps)
        
        dmdt = 0.975 * param.cs**3 / param.G
        
        mass = [0]
        for i in range(len(time)):
            if i==0:
                dt = time[i]
            else:
                dt = time[i]-time[i-1]

            mass.append(mass[-1]+dmdt*dt)
        
        mass.pop(0)

        # print([i/param.M for i in mass])

        return mass

    def _find_closest_index(self, lst, target_value):
            closest_index = min(range(len(lst)), key=lambda i: abs(lst[i] - target_value))

            return closest_index        

class cassen_graph():
    def __init__(self, r, r_grid, save=False):
        self.r = r
        self.tau_2  = [(param.omega*t)**2* param.cs*t /param.AU for t in param.time]
        self.r_c    = [param.cs/param.omega /param.AU for t in param.time]
        self.r_d    = [(param.omega**2*param.G**3*param.M**3) / (16*param.cs**8) / param.AU for t in param.time]
        self.r_env  = [param.G*param.M/(2*param.cs**2) / param.AU for t in param.time]
        self.centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 / param.AU for t in param.time]

        self.data_param = ["CM", param.cs, param.omega, len(param.time), param.inner_sample_pts, param.inner_theta_pts]

        if os.path.exists(data._get_directory(self.data_param)):
            print('Loading CM81 data...')
            self.radius, self.r_grid, self.nh_density, self.cgs_v_r, self.cgs_v_theta, self.cgs_v_phi = data.load(self.data_param)

        else:
            self.r = r
            self.r_grid = r_grid
            self.radius = []
            for i in range(len(param.time)):
                self.radius.append([j/param.AU for j in self.r[i]])
                print(self.tau_2[i]/1e3, self.radius[i][-1]/1e3)

            self.cm_envelope = cm81_envelope(self.r)

            self.nh_density     = self.cm_envelope.get_density()
            self.cgs_v_r        = self.cm_envelope.get_u_r()
            self.cgs_v_theta    = self.cm_envelope.get_u_theta()
            self.cgs_v_phi      = []

            if save == True:
                data.save(param=self.data_param, radius=self.radius, r_grid=self.r_grid, density=self.nh_density, 
                          v_r=self.cgs_v_r, v_theta=self.cgs_v_theta, v_phi=self.cgs_v_phi)
                print('Successfully stored data')
        
    def v_r_envelope(self):
        theta = 0
        for time in range(len(param.time)):
            plt.plot(self.radius[time], [-i for i in self.cgs_v_r[time][theta]])

    def nh_density_envelope_1d(self):
        theta = 0
        for time in range(len(param.time)):
            if time != 1:
                continue
            for theta in range(param.inner_theta_pts):
                # if theta == 1:
                    # print(theta)
                # if theta<param.inner_theta_pts-100:
                #     continue
                plt.plot([i*param.AU for i in self.radius[time]], [abs(i) for i in self.nh_density[time][theta]])
                    # plt.axvline(self.centrifugal_radius[time])

    def density_envelope_2d(self):
        fig, axes = plt.subplots(2, 2)
        axes_flat = axes.flatten()

        theta_values = np.linspace(param.min_theta, np.pi/2, param.inner_theta_pts)

        i = 0
        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"graphing progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                
                molecule_counts = [] 
                for theta in range(len(theta_values)):
                    molecule_counts.append([abs(j) for j in self.nh_density[i][theta]])

                # Example data in (r, theta) format
                r_values =  np.array([j/1e3 for j in self.radius[i]])
                
                # Convert (r, theta) to (x, y) coordinates
                x_values = []
                y_values = []
                for j, theta in enumerate(theta_values):
                    x_values.extend(r_values * np.cos(theta))
                    y_values.extend(r_values * np.sin(theta))

                # Create a flattened (x, y) array and corresponding flattened data array
                x_flat = np.array(x_values)
                y_flat = np.array(y_values)
                data_flat = np.array(molecule_counts).flatten()

                # Create a grid in (x, y) coordinates
                x_grid = [j/1e3 for j in self.radius[i]]
                y_grid = [j/1e3 for j in self.radius[i]]
                xx, yy = np.meshgrid(x_grid, y_grid)

                # Interpolate data onto the (x, y) grid
                data_interp = griddata((x_flat, y_flat), data_flat, (xx, yy), method='linear')

                # Create a contour plot with a logarithmic colormap
                norm = LogNorm(vmin=param.visser_vmin, vmax=param.visser_vmax)
                im = ax.contourf(xx, yy, data_interp, cmap='gist_rainbow_r', levels=np.logspace(4,9,100), norm=norm, vmin=param.visser_vmin, vmax=param.visser_vmax)

                ax.set_xlim(0, 10)
                ax.set_ylim(0, 10)

                r_env_edge = patches.Circle((0,0), self.r_env[i]/1e3, fill=False, edgecolor='white', linewidth=2)
                ax.add_patch(r_env_edge)
                circ_tau_2 = patches.Circle((0,0), self.tau_2[i]/1e3, fill=False, edgecolor='black', linewidth=2)
                ax.add_patch(circ_tau_2)
                circ_r_c = patches.Circle((0,0), self.r_c[i]/1e3, fill=False, edgecolor='red', linewidth=2)
                ax.add_patch(circ_r_c)
                circ_r_d = patches.Circle((0,0), self.r_d[i]/1e3, fill=False, edgecolor='blue', linewidth=2)
                ax.add_patch(circ_r_d)
                circ_centrifugal_radius = patches.Circle((0,0), self.centrifugal_radius[i]/1e3, fill=False, edgecolor='green', linewidth=2)
                ax.add_patch(circ_centrifugal_radius)
                circ_5_AU = patches.Circle((0,0), 5/1e3, fill=False, edgecolor='yellow', linewidth=2)
                ax.add_patch(circ_5_AU)

                ax.set_title(f't = {param.time[i]/param.t_acc}')
                
            else:
                plt.axis('off')
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_ticks([10**4, 10**5, 10**6, 10**7, 10**8, 10**9])
                cbar.set_label('Density')
    
        fig.suptitle('CM (1981) 2D Graph (n_h)')

    def cgs_mass_env_to_disk(self, disk_inner_radius, time=param.time):
        r_star = []
        # print([i/31556926 for i in time])
        for t in [i/31556926 for i in time]:
            if t<20000:
                r_star.append(5*param.AU)
            elif t>20000 and t<20100:
                print("here")
                r_star.append(5*param.AU*(1-np.sqrt((t-20000)/100)) + param.r_ps)
            else:
                r_star.append(param.r_ps)
        
        self.centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 / param.AU for t in param.time]

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        mass = [0]
        for t in range(len(param.time)):
            r_c_index = self._find_closest_index(self.radius[t], self.centrifugal_radius[t])

            if t==0:
                dt = time[t]
            if t!=0:
                dt = time[t]-time[t-1]
            
            dmdt = 0
            theta = -1
            for i in range(r_c_index):
                if param.AU*self.radius[t][i] < disk_inner_radius[t]:
                    continue
                if param.AU*self.radius[t][i] > disk_inner_radius[t]:
                    if i == 0:
                        dr = param.AU*self.radius[t][i]
                    else:
                        dr = param.AU*(self.radius[t][i]-self.radius[t][i-1])

                    radius = param.AU*self.radius[t][i]
                    density = np.real(self.nh_density[t][theta][i]*param.M_H)
                    velocity = (abs(np.real(self.cgs_v_r[t][theta][i]))*np.sin(theta_values[theta]) + abs(np.real(self.cgs_v_theta[t][theta][i]))*np.cos(theta_values[theta]))

                    dmdt += radius * density * velocity * dr

            mass.append(mass[-1]+dmdt*dt)
        
        mass.pop(0)

        # print([i/param.M for i in mass])

        return mass

    def cgs_mass_env_to_star(self, time=param.time):
        r_star = []
        # print([i/31556926 for i in time])
        for t in [i/31556926 for i in time]:
            if t<20000:
                r_star.append(5*param.AU)
            elif t>20000 and t<20100:
                print("here")
                r_star.append(5*param.AU*(1-np.sqrt((t-20000)/100)) + param.r_ps)
            else:
                r_star.append(param.r_ps)

        # print([i/param.AU for i in r_star])
                
        dmdt = 0.975 * param.cs**3 / param.G
        self.centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 for t in time]

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        caveat_angle = [10**(1.5)*(t/param.t_acc)**(0.26) * np.pi/180 for t in param.time]
        for i in range(len(caveat_angle)):
            caveat_angle[i] = self._find_closest_index(theta_values, caveat_angle[i])

        mass = [0]
        for t in range(len(time)):
            r_star_index = self._find_closest_index(self.radius[t], r_star[t])

            if t==0:
                dt = time[t]
            if t!=0:
                dt = time[t]-time[t-1]

            if self.centrifugal_radius[t] < r_star[t]:
                mass.append(mass[-1]+dmdt*dt)
                print(t, "r_star", r_star[t]/param.AU, "mass", mass[-1]/param.M, "dm", dmdt*dt/param.M)

            if self.centrifugal_radius[t] >= r_star[t]:
                dmdt = 0
                for theta in range(len(theta_values)):
                    if theta < caveat_angle[t]:
                        continue

                    density = abs(self.nh_density[t][theta][r_star_index])*param.M_H
                    velocity = abs(self.cgs_v_theta[t][theta][r_star_index])
                    radius = r_star[t]
                    d_theta = np.pi/2 / param.inner_theta_pts

                    dmdt += 4*np.pi * radius * density * velocity * d_theta

                    if t%10==0 and theta%100==0:
                        print(t, theta)
                        print("r_star", r_star[t]/param.AU)
                        print("velocity", velocity)
                        print("density", density)
                        print("d_theta", d_theta)
                        print("dm", dmdt*dt/param.M)
                        print("mass", (mass[-1]+dmdt*dt)/param.M)
                        print("\n")

                mass.append(mass[-1]+dmdt*dt) 

            # print(mass[-1]/param.M)
        
        mass.pop(0)

        # print([i/param.M for i in mass])

        return mass

    def _find_closest_index(self, lst, target_value):
            closest_index = min(range(len(lst)), key=lambda i: abs(lst[i] - target_value))

            return closest_index

class star():
    def __init__(self):
        return
    
    def get_r_star(self, time):
        r_star = []
        for t in [i/31556926 for i in time]:
                if t<20000:
                    r_star.append(5*param.AU)
                elif t>20000 and t<20100:
                    print("here")
                    r_star.append(5*param.AU*(1-np.sqrt((t-20000)/100)) + param.r_ps)
                else:
                    r_star.append(param.r_ps)
    
        # plt.plot([i/param.sun_R for i in r_star])
        return r_star

    def luminosity_accretion(self, time, mass):
        def m_d_ode(u):
            return np.sqrt(1-u) / u**(4/3)
        
        r_star = self.get_r_star(time)
        centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 for t in time]
        
        l_acc = []
        for t in range(len(time)):
            dmdt = 0.975 * param.cs**3 / param.G
            l0 = param.G * mass[t] * dmdt / r_star[t]
            u_star = r_star[t] / centrifugal_radius[t]
            print(u_star)

            eta_star    = 0.5
            eta_d       = 0.75

            print(quad(m_d_ode, u_star, 1))
            m_d = 1/3 * u_star**(1/3) * quad(m_d_ode, u_star, 1)[0]
            # print(quad(m_d_ode, u_star, 1))

            luminosity = 1/(6*u_star) * (6*u_star - 2 + (2-5*u_star)*np.sqrt(1-u_star))
            luminosity += eta_star/2 * (1 - (1-eta_d)*m_d) * (1 - (1-eta_d)*np.sqrt(1-u_star))
            luminosity *= l0

            l_acc.append(luminosity)
            
        # print([i/(3.846*1e33) for i in l_acc])
        
        return l_acc
    
    def luminosity(self, time, mass):
        return self.luminosity_accretion(time, mass)

    def disk_inner_radius(self, time, mass):
        t_evap = 2000
        l_acc = self.luminosity(time, mass)

        r = []
        for t in range(len(time)):
            r.append(np.sqrt(l_acc[t] / (4*np.pi * param.bolzmann * t_evap**4)))
        
        print([i/param.AU for i in r])
        # plt.plot([i/param.AU for i in r])

        return r

import numpy as np
import matplotlib.pyplot as plt

from data_handler import data
from parameter import param
from graphing_help import shu_graph, tsc_graph, cassen_graph, star

param = param()
data = data()

# make grid
r_env  = [param.G*param.M/(2*param.cs**2) for t in param.time]
r_star = []
for t in [i/31556926 for i in param.time]:
            if t<20000:
                r_star.append(5*param.AU)
            elif t>20000 and t<20100:
                r_star.append(5*param.AU*(1-np.sqrt((t-20000)/100)) + param.r_ps)
            else:
                r_star.append(param.r_ps)

sc_inner_r = "{:e}".format(0.1*param.AU)
a, inner_r_exp = sc_inner_r.split('e')
sc_r_env  = "{:e}".format(r_env[0])
a, r_env_exp = sc_r_env.split('e')

r_grid = []
for t in range(len(param.time)):
    r_grid.append(np.logspace(int(inner_r_exp), int(r_env_exp), 2*param.sample_pts+1))

r = []
for t in range(len(param.time)):
    r.append([])
    for i in range(2*param.sample_pts+1):
        if i == 0:
            r_i = 0.5*r_grid[t][i]
            r[-1].append(r_i)
        else:
            r_i = 0.5*(r_grid[t][i-1] + r_grid[t][i])
            r[-1].append(r_i)

r_temp = r
r = r_grid
r_grid = r_temp

for t in range(len(param.time)):
    r[t] = r[t][:-1]

#######################

z_g1_grid = [2/(i+1) for i in np.logspace(param.sample_range, 0, param.sample_pts)]  # also x_s1
z_s1_grid = np.logspace(0, param.sample_range, param.sample_pts)

z_g1 = []
z_s1 = []
for i in range(param.sample_pts):
    if i == 0:
        z_s1_i = z_s1_grid[i]
        z_g1_i = 0.5*z_g1_grid[i]
        z_s1.append(z_s1_i)
        z_g1.append(z_g1_i)
    
    else:
        z_s1_i = 0.5*(z_s1_grid[i-1] + z_s1_grid[i])
        z_g1_i = 0.5*(z_g1_grid[i-1] + z_g1_grid[i])
        z_s1.append(z_s1_i)
        z_g1.append(z_g1_i)
z_g1[-1] =1-param.min_pts

# # print(z_s1_grid[0], z_s1_grid[-1], z_g1_grid[0], z_g1_grid[-1])
# # print(z_s1[0], z_s1[-1], z_g1[0], z_g1[-1])

r_z = []
for t in param.time:
    r_z.append([1/i*param.cs*t for i in z_s1[::-1]]+[1/i*param.cs*t for i in z_g1[::-1]])

for t in range(len(param.time)):
    print("r_z", np.array(r_z[t]))
    print("r", np.array(r[t]))

# z = []
# for t in param.time:
#     z.append(z_s1[::-1]+ z_g1[::-1])

# r_grid = []
# for t in param.time:
#     r_grid.append([1/i*param.cs*t for i in z_s1_grid[::-1]]+[1/i*param.cs*t for i in z_g1_grid[::-1]])

# r_temp = []
# r_grid_temp = []
# for i in range(len(param.time)):
#     r_temp.append([])
#     r_grid_temp.append([])
#     for j in range(len(r[i])):
#         # at the end, add an additional interval 
#         if j == len(r[i])-1:
#             r_temp[-1].append(r[i][j])
#             r_grid_temp[-1].append(r_grid[i][j])
#             r_grid_temp[-1].append(2*r[i][j]-r_grid[i][j])

#         # avoid the inconsistnecy at x=1
#         elif j+1 == param.sample_pts:
#             r_temp[-1].append((r_grid[i][j]+r_grid[i][j+2])/2)
#             r_grid_temp[-1].append(r_grid[i][j+1])
#         elif j == param.sample_pts:
#             continue

#         else:
#             r_temp[-1].append(r[i][j])
#             r_grid_temp[-1].append(r_grid[i][j])

# r = r_temp
# r_grid = r_grid_temp
#################################

# unit test
for i in range(len(param.time)):
    for j in range(len(r[i])):
        if r[i][j]-r_grid[i][j] < 0 or r_grid[i][j+1]-r[i][j]<0:
            print(j, 'ERRORRRRR')
            for k in range(6, 0, -1):
                print(j-k, r_grid[i][j-k])
            print(j, r_grid[i][j], r[i][j], r_grid[i][j+1])
            for k in range(1, 6):
                print(j+k, r_grid[i][j+k])

# fig, axs = plt.subplots(1, 3, figsize=(12, 4))  # 1 row, 3 columns

# Plot the data in each subplot
# for i in range(len(param.time)):
#     print("grid:", r_grid[i][0], r_grid[i][-1])
#     print("radius:", r[i][0], r[i][-1])
#     axs[i].plot(r_grid[i], '--bo', label='List 1', color='blue')
#     axs[i].plot(r[i], '--bo', label='List 2', color='green')

#     axs[i].set_yscale('log')
#     axs[i].set_xscale('log')

# plt.plot(r_grid[0],'--bo', color='red')
# plt.plot(r[0], '--bo', color='blue')

# CM81 Radius
# tau_2 = [(param.omega*t)**2* param.cs*t /param.AU for t in param.time]
# [(param.omega*t)**2* param.cs*t /param.AU for t in param.time]
# # print(tau_2)

# cm_r_grid = []
# for time in range(len(param.time)):
#     cm_r_grid.append(np.linspace(param.min_pts, tau_2[time]*param.AU, param.inner_sample_pts))

# cm_r = []
# for time in range(len(param.time)):
#     cm_r.append([])
#     for i in range(param.inner_sample_pts):
#         if i == 0:
#             continue
        
#         else:
#             cm_r_i = 0.5*(cm_r_grid[time][i-1]+cm_r_grid[time][i])
#             cm_r[-1].append(cm_r_i)

# print("here", [cm_r_grid[i][-1]/param.AU for i in range(len(param.time))])
# print("here", [cm_r[i][-1]/param.AU for i in range(len(param.time))])

# plt.plot(cm_r_grid[0], '--bo', color='red')
# plt.plot(cm_r[0], '--bo', color='blue')

# integrate volume 
# r_env = [param.G*param.M/(2*param.cs**2) / param.AU for t in param.time]

# theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

# v = []
# for t in range(len(param.time)):
#     v_t = 0
#     for theta in range(param.theta_pts):
#         if np.sin(theta_values[theta]) < 0:
#             print(np.sin(theta_values[theta]))
#             print(theta_values[theta])
                
#         for i in range(param.sample_pts):
            
#             if r[t][i+1] < r_env[t]*param.AU:    
#                 v_t += r[t][i+1]**2 * (r_grid[t][i+1]-r_grid[t][i]) * np.sin(theta_values[theta])
    
#     v.append(4 * np.pi * v_t * np.pi/2 / param.theta_pts)

# print("numerical volume:", v)
# print("volume:", [4/3*np.pi * (r_env[t]*param.AU-r_grid[t][0])**3 for t in range(len(param.time))])
# print("ratio:", [v[t] / (4/3*np.pi * (r_env[t]*param.AU-r_grid[t][0])**3) for t in range(len(param.time))])

########## Shu (1977) Graph ##########
# shu = shu_graph(z_s1=z_s1, z_g1=z_g1, r_grid=r_grid, save=True)
# shu.dimensionless_a_v_1d()
# shu.cgs_density_1d()
# shu.nh_density_1d()
# shu.cgs_density_2d()
# shu.nh_density_2d()



########## TSC (1984) Graph ##########
# print(len(r), len(r[0]))
t = tsc_graph(z_s1=z_s1, z_g1=z_g1, r_grid=r_grid, r=r, save=False)
# t.dimensionless_a_1d()
# t.dimensionless_quad_v_w_1d()
# t.dimensionless_mono_v_1d()
# t.cgs_density_1d()
# t.nh_density_1d()
# t.cgs_density_2d()
# t.nh_density_2d()
# t.int_envelope_mass()
# t.cgs_mass_envelope()

# env_to_star_disk_mass = t.cgs_mass_env_to_star_disk(time=param.time)
# plt.plot([i/param.t_acc for i in param.time], [1 - i/param.M for i in env_to_star_disk_mass])

# mass_disk = []
# for i in range(len(param.time)):
#     index = int(param.time[i]/param.t_acc * 100)
#     mass_disk.append(mass[index-1]-env_to_star_mass[i])
# plt.plot([i/param.t_acc for i in param.time], [i/param.M for i in mass_disk])
# t.cgs_sigma_disk()

################ Star ################
# star = star()
# star.get_r_star(time)
# l_acc = star.luminosity_accretion(param.time, env_to_star_disk_mass)
# plt.plot([i/param.t_acc for i in param.time], [i/(3.846*1e33) for i in l_acc])
# disk_inner_radius = star.disk_inner_radius(param.time, env_to_star_disk_mass)

########## CM (1981) Graph ###########
# cm = cassen_graph(r=r, r_grid=r_grid, save=True)
# cm.v_r_envelope()
# cm.nh_density_envelope_1d()
# cm.density_envelope_2d()

# env_to_disk_mass = cm.cgs_mass_env_to_disk(time=param.time, disk_inner_radius=disk_inner_radius)
# plt.plot([i/param.t_acc for i in param.time], [i/param.M for i in env_to_disk_mass])
# env_to_star_mass = cm.cgs_mass_env_to_star()
# plt.plot([i/param.t_acc for i in param.time], [i/param.M for i in env_to_star_mass])

plt.xscale('log')
plt.yscale('log')

plt.show()