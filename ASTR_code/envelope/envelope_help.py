import numpy as np
import os
import sys

sys.path.insert(1, '/Users/haleylin/Desktop/ASTR Sphinx/ASTR_code/')
from parameter import param
from data_handler import data

param = param()
data = data()

class envelope():
    def __init__(self):
        # import TSC (envelope data)
        data_param = ["TSC", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]

        if os.path.exists(data._get_directory(data_param)):
            print('Loading TSC84 data...')
            self.radius, r_grid, self.tsc_nh_density, self.tsc_cgs_v_r, self.tsc_cgs_v_theta, cgs_v_phi = data.load(data_param)
        else:
            print("Data does not exist!")

        # import CM81 (envelope data)
        data_param = ["CM", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]
        if os.path.exists(data._get_directory(data_param)):
            print('Loading CM81 data...')
            radius, r_grid, self.cm_nh_density, self.cm_cgs_v_r, self.cm_cgs_v_theta, cgs_v_phi = data.load(data_param)
        else:
            print("Data does not exist!")
    
    def cgs_env_mass(self):
        # the mass calculate here is the mass infalling within centrifugal radius
        # but then, we subtract it by initial mass, therefore having the evolution of envelope mass    
        dmdt = 0.975 * param.cs**3 / param.G
        
        mass = [0]
        for i in range(len(param.time)):
            if i==0:
                dt = param.time[i]
            else:
                dt = param.time[i]-param.time[i-1]

            mass.append(mass[-1]+dmdt*dt)
        
        mass.pop(0)

        # print([i/param.M for i in mass])

        return [1 - i/param.M for i in mass]

    def cgs_env_to_disk_mass(self, disk_inner_radius): 
        self.centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 / param.AU for t in param.time]

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        mass = [0]
        for t in range(len(param.time)):
            r_c_index = self._find_closest_index(self.radius[t], self.centrifugal_radius[t])

            if t==0:
                dt = param.time[t]
            if t!=0:
                dt = param.time[t]-param.time[t-1]
            
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
                    density = np.real(self.cm_nh_density[t][theta][i]*param.M_H)
                    velocity = (abs(np.real(self.cm_cgs_v_r[t][theta][i]))*np.sin(theta_values[theta]) + abs(np.real(self.cm_cgs_v_theta[t][theta][i]))*np.cos(theta_values[theta]))

                    dmdt += radius * density * velocity * dr

            mass.append(mass[-1]+dmdt*dt)
        
        mass.pop(0)

        # print([i/param.M for i in mass])

        return mass

    def cgs_env_to_star_mass(self, r_star):                
        dmdt = 0.975 * param.cs**3 / param.G
        self.centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 for t in param.time]

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        caveat_angle = [10**(1.5)*(t/param.t_acc)**(0.26) * np.pi/180 for t in param.time]
        for i in range(len(caveat_angle)):
            caveat_angle[i] = self._find_closest_index(theta_values, caveat_angle[i])

        mass = [0]
        for t in range(len(param.time)):
            r_star_index = self._find_closest_index([i*param.AU for i in self.radius[t]], r_star[t])

            if t==0:
                dt = param.time[t]
            if t!=0:
                dt = param.time[t]-param.time[t-1]

            if self.centrifugal_radius[t] < r_star[t]:
                mass.append(mass[-1]+dmdt*dt)
                # print(t, "r_star", r_star[t]/param.AU, "mass", mass[-1]/param.M, "dm", dmdt*dt/param.M)

            if self.centrifugal_radius[t] >= r_star[t]:
                dmdt = 0
                for theta in range(len(theta_values)):
                    if theta < caveat_angle[t]:
                        continue
                    
                    # print(r_star_index)
                    density = abs(self.cm_nh_density[t][theta][r_star_index])*param.M_H
                    velocity = abs(self.cm_cgs_v_theta[t][theta][r_star_index])
                    radius = r_star[t]
                    d_theta = np.pi/2 / param.theta_pts

                    dmdt += 4*np.pi * radius * density * velocity * d_theta

                    # if t%10==0 and theta%100==0:
                    #     print(t, theta)
                    #     print("r_star", r_star[t]/param.AU)
                    #     print("velocity", velocity)
                    #     print("density", density)
                    #     print("d_theta", d_theta)
                    #     print("dm", dmdt*dt/param.M)
                    #     print("mass", (mass[-1]+dmdt*dt)/param.M)
                    #     print("\n")

                mass.append(mass[-1]+dmdt*dt) 

            # print(mass[-1]/param.M)
        
        mass.pop(0)

        # print([i/param.M for i in mass])

        return mass
    
    def _find_closest_index(self, lst, target_value):
            closest_index = min(range(len(lst)), key=lambda i: abs(lst[i] - target_value))

            return closest_index
