import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
import numpy as np
import os

from envelope.shu77 import shu
from envelope.tsc84 import tsc
from data_handler import data
from parameter import param

# TODO: save data using npz
#       -> create folder for each different param + a csv record param
# TODO: integrate and check total mass

param = param()
data = data()

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
    def __init__(self, z_s1, z_g1, r_grid, save=False):
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

        self.inner_edge     = [(param.omega*t)**2* param.cs*t/param.AU for t in param.time]
        self.outer_edge     = [1/(param.omega*t) * param.cs*t/param.AU for t in param.time]

        self.data_param = ["TSC", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]

        if os.path.exists(data._get_directory(self.data_param)):
            print('Loading data...')
            self.radius, self.r_grid, self.nh_density, self.cgs_v_r, self.cgs_v_theta, self.cgs_v_phi = data.load(self.data_param)

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

            print("Initializing density")
            self.cgs_density_s1 = self.tsc.get_cgs_density_s1()
            self.cgs_density_g1 = self.tsc.get_cgs_density_g1()
            self.nh_density_s1  = self.tsc.get_nh_density_s1()
            self.nh_density_g1  = self.tsc.get_nh_density_g1()

            print("Initializing cgs velocity")
            self.cgs_v_r_s1     = self.tsc.get_cgs_v_r_s1()
            self.cgs_v_r_g1     = self.tsc.get_cgs_v_r_g1()
            self.cgs_v_theta_s1 = self.tsc.get_cgs_v_theta_s1()
            self.cgs_v_theta_g1 = self.tsc.get_cgs_v_theta_g1()
            self.cgs_v_phi      = self.tsc.get_cgs_v_phi()

            print("Initializing a, v, w")
            self.quad_a_g1      = self.tsc.get_quad_a_g1()
            self.quad_a_s1      = self.tsc.get_quad_a_s1()
            self.quad_v_g1      = self.tsc.get_quad_v_g1()
            self.quad_v_s1      = self.tsc.get_quad_v_s1()
            self.quad_w_g1      = self.tsc.get_quad_w_g1()
            self.quad_w_s1      = self.tsc.get_quad_w_s1()
            self.mono_a_g1      = self.tsc.get_mono_a_g1()
            self.mono_a_s1      = self.tsc.get_mono_a_s1()
            self.mono_v_s1      = self.tsc.get_mono_v_s1()

            # print(param.omega*param.M*param.G/0.975/param.cs**3)
            # print([param.omega*t for t in param.time])

            self.radius         = []
            self.nh_density     = []
            self.cgs_v_r        = []
            self.cgs_v_theta    = []
            self.cgs_v_phi      = self.cgs_v_phi

            for i in range(len(param.time)):
                print(f"reshape process, {i+1}/{len(param.time)}" )
                self.radius.append([j/param.AU for j in self.radius_s1[i]][::-1]+[j/param.AU for j in self.radius_g1[i]][::-1])

                self.nh_density.append([])
                self.cgs_v_r.append([])
                self.cgs_v_theta.append([])
                for j in range(param.theta_pts):
                    self.nh_density[-1].append(self.nh_density_s1[i][j][::-1]+self.nh_density_g1[i][j][::-1])
                    self.cgs_v_r[-1].append(self.cgs_v_r_s1[i][j][::-1]+self.cgs_v_r_g1[i][j][::-1])
                    self.cgs_v_theta[-1].append(self.cgs_v_theta_s1[i][j][::-1]+self.cgs_v_theta_g1[i][j][::-1])
                
            print(len(r_grid[0]), len(self.radius[0]))

            print("radius",     len(self.radius),       len(self.radius[0]))                # [time][radius]
            print("density",    len(self.nh_density),   len(self.nh_density[0][0]))         # [time][theta][radius]
            print("v_r",        len(self.cgs_v_r),      len(self.cgs_v_r[0][0]))            # [time][theta][radius]
            print("v_theta",    len(self.cgs_v_theta),  len(self.cgs_v_theta[0][0]))        # [time][theta][radius]
            print("v_phi",      len(self.cgs_v_phi),    len(self.cgs_v_phi[0]))             # [time][theta]

            if self.save == True:
                data.save(param=self.data_param, radius=self.radius, r_grid=self.r_grid, density=self.nh_density, 
                          v_r=self.cgs_v_r, v_theta=self.cgs_v_theta, v_phi=self.cgs_v_phi)
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
            plt.plot(self.radius[i], [j*param.M_H for j in self.nh_density[i][0]], color=param.color_1d[i], label=f't={param.time[i]/1e12}e12')
    
    def nh_density_1d(self):
        """
        Plot density profiles in density of hydrogen in 1D for different times.
        """
        for i in range(len(param.time)):
            plt.plot(self.radius[i], self.nh_density[i][0], color=param.color_1d[i], label=f't={param.time[i]/1e12}e12')

    def cgs_density_2d(self):
        """
        Plot density profiles in CGS units in 2D for different times using 
        a logarithmic colormap with contour value.
        """
        fig, axes = plt.subplots(2, 2)
        axes_flat = axes.flatten()

        theta_values = np.linspace(param.min_theta, 2 * np.pi, param.theta_pts)

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

        theta_values = np.linspace(param.min_theta, 2*np.pi, param.theta_pts)

        i = 0
        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                
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

                # Interpolate data onto the (x, y) grid
                data_interp = griddata((x_flat, y_flat), data_flat, (xx, yy), method='linear')

                # Create a contour plot with a logarithmic colormap
                norm = LogNorm(vmin=param.visser_vmin, vmax=param.visser_vmax)
                im = ax.contourf(xx, yy, data_interp, cmap='gist_rainbow_r', levels=np.logspace(4,9,100), norm=norm)

                self.inner_edge     = [(param.omega*t)**2* param.cs*t/param.AU for t in param.time]
                self.outer_edge     = [1/(param.omega*t) * param.cs*t/param.AU for t in param.time]

                ax.set_xlim(0, 10)
                ax.set_ylim(0, 10)

                # Add contour line at the value
                ax.contour(xx, yy, data_interp, levels=[param.visser_contour],           
                           colors='red', linewidths=2)
                
                r_env_edge = patches.Circle((0,0), param.r_env/1e3, fill=False, edgecolor='white', linewidth=2)
                ax.add_patch(r_env_edge)
                circ_inner = patches.Circle((0,0), self.inner_edge[i]/1e3, fill=False, edgecolor='black', linewidth=2)
                ax.add_patch(circ_inner)
                circ_outer = patches.Circle((0,0), self.outer_edge[i], fill=False, edgecolor='black', linewidth=2)
                ax.add_patch(circ_outer)

                ax.set_title(f't = {param.time[i]/param.t_acc}')
                
            else:
                plt.axis('off')
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_ticks([10**4, 10**5, 10**6, 10**7, 10**8, 10**9])
                cbar.set_label('Density')
    
        fig.suptitle('TSC (1984) 2D Graph (n_h)')
    
    def int_mass(self):
        print("Integrating...")

        # times AU for all self.radius
        radius_AU = [[value * param.AU for value in row] for row in self.radius]
        outer_edge = [i*param.AU for i in self.outer_edge]
        inner_edge = [i*param.AU for i in self.inner_edge]
        
        # 2 * 2*pi S(rho * r^2)dr * (pi/2 /theta pts)
        self.mass = [] 
        for t in range(len(param.time)):
            print(f"Integrating process: {t+1}/{len(param.time)}")
            mass_t = 0 

            for theta in range(param.theta_pts):
                for i in range(param.sample_pts):
                    # unit test
                    if radius_AU[t][i]-self.r_grid[t][i] < 0 or self.r_grid[t][i+1]-radius_AU[t][i]<0:
                        print(self.r_grid[t][i], radius_AU[t][i], self.r_grid[t][i+1])
                    
                    # print(radius_AU[t][i], inner_edge[t], outer_edge[t])

                    if radius_AU[t][i] > inner_edge[t] and radius_AU[t][i] < outer_edge[t]:
                        if param.M_H*self.nh_density[t][theta][i]<0:
                            print(t, theta, i, param.M_H*self.nh_density[t][theta][i])
                        
                        if param.M_H*self.nh_density[t][theta][i]>0:
                            mass_t += param.M_H*self.nh_density[t][theta][i] * radius_AU[t][i]**2 * (self.r_grid[t][i+1]-self.r_grid[t][i])
           
            self.mass.append(4*np.pi*mass_t * np.pi/2 / param.theta_pts)

        print(param.M, self.mass)
        print(self.inner_edge, self.outer_edge)