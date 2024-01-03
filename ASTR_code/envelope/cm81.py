import numpy as np
import matplotlib.pyplot as plt
from scipy.special import legendre
from tqdm import tqdm
import sys

sys.path.insert(1, '/Users/haleylin/Desktop/ASTR Sphinx/ASTR_code/')
from parameter import param

param = param()

class cm81_envelope():
    def __init__(self, r):
        self.r = r
        self.r_c = self.get_centrifugal_radius()

    def get_centrifugal_radius(self):
        return [1/16*param.cs*(0.975)**3 * t**3 * param.omega**2 for t in param.time]
    
    def get_u_r(self):
        r_c = self.get_centrifugal_radius()
        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        # self.u_r = []
        self.u_r = np.zeros((len(param.time), param.theta_pts, len(self.r[0])))
        for time in tqdm(range(len(param.time)), desc='get_u_r'):
            # self.u_r.append([])
            for theta in range(param.theta_pts):
                # self.u_r[-1].append([])
                for i in range(len(self.r[time])):
                    # print([r_c[time]/self.r[time][i], 1-(r_c[time]/self.r[time][i]), 0, -np.cos(theta_values[theta])])
                    theta_0 = np.arccos(np.roots([r_c[time]/self.r[time][i], 1-(r_c[time]/self.r[time][i]), 0, -np.cos(theta_values[theta])])[2])
                    u_r_i = -np.sqrt(param.G*param.M /self.r[time][i]) * np.sqrt(1+np.cos(theta_values[theta])/theta_0)

                    self.u_r[time][theta][i] = u_r_i
                    # self.u_r[-1][-1].append(u_r_i)
        return self.u_r
    
    def get_u_theta(self):
        r_c = self.get_centrifugal_radius()
        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)

        # self.u_r = []
        self.u_theta = np.zeros((len(param.time), param.theta_pts, len(self.r[0])))
        for time in tqdm(range(len(param.time)), desc='get_u_theta'):
            # self.u_r.append([])
            for theta in range(param.theta_pts):
                # self.u_r[-1].append([])
                for i in range(len(self.r[time])):
                    # print([r_c[time]/self.r[time][i], 1-(r_c[time]/self.r[time][i]), 0, -np.cos(theta_values[theta])])
                    theta_0 = np.arccos(np.roots([r_c[time]/self.r[time][i], 1-(r_c[time]/self.r[time][i]), 0, -np.cos(theta_values[theta])])[2])
                    u_theta_i = -np.sqrt(param.G*param.M /self.r[time][i]) * np.sqrt(1+np.cos(theta_values[theta])/theta_0) * (np.cos(theta_0)-np.cos(theta))/np.sin(theta)

                    self.u_theta[time][theta][i] = u_theta_i
                    # self.u_r[-1][-1].append(u_r_i)
        return self.u_theta
    
    def get_density(self):
        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)
        leg = legendre(2)

        dmdt = (0.975)*param.cs**3/param.G
        r_c = self.get_centrifugal_radius()
        u_r  = self.get_u_r()
        
        self.density = []
        for time in tqdm(range(len(param.time)), desc='get_density'):
            self.density.append([])
            for theta in range(param.theta_pts):
                self.density[-1].append([])
                for i in range(len(self.r[time])):
                    theta_0 = np.arccos(np.roots([r_c[time]/self.r[time][i], 1-(r_c[time]/self.r[time][i]), 0, -np.cos(theta_values[theta])])[2])
                    density_i = dmdt/(4*np.pi*self.r[time][i]**2 * u_r[time][theta][i]) / (1+2*r_c[time]/self.r[time][i]*leg(np.cos(theta_0)))

                    self.density[-1][-1].append(density_i/param.M_H)

        return self.density
    
    def dmdt_envelope(self, u_r, density):
        '''
        Calculate dmdt of envelope accrete onto star by 

            2 * 2*pi S_[0,pi/2] (r^2 * sin(theta) * density * u_r) d(theta) 

        Input:
            u_r [cm/s]
            density [nh/cm^3]

        Ooutput:
            dmdt[time]
        '''

        theta_values = np.linspace(param.min_theta, np.pi/2, param.theta_pts)
        
        dmdt = []
        for t in range(len(param.time)):
            closest_r = self._find_closest_index(self.r[t], self.r_c[t])

            dmdt_i = 0
            for theta in range(param.theta_pts):   
                dmdt_i += self.r[t][closest_r]**2 * np.sin(theta_values[theta]) * density[t][theta][closest_r]*param.M_H * (-u_r[t][theta][closest_r])
                
                # if t==0 and theta%100 == 0:
                #     print(dmdt_i, u_r[t][theta][closest_r])
                    
            dmdt.append(4*np.pi * dmdt_i * np.pi/2/param.theta_pts)
            
        return dmdt
    
    def _find_closest_index(self, lst, target_value):
        closest_index = min(range(len(lst)), key=lambda i: abs(lst[i] - target_value))

        return closest_index


# tau_2 = [(param.omega*t)**2* param.cs*t /param.AU for t in param.time]
# print(tau_2)

# r_grid = []
# for time in range(len(param.time)):
#     r_grid.append(np.linspace(param.min_pts, tau_2[time]*param.AU, param.inner_sample_pts))

# r = []
# for time in range(len(param.time)):
#     r.append([])
#     for i in range(param.inner_sample_pts):
#         if i == 0:
#             continue
        
#         else:
#             r_i = 0.5*(r_grid[time][i-1]+r_grid[time][i])
#             r[-1].append(r_i)

# plt.plot(r_grid[0], '--bo', color='red')
# plt.plot(r[0], '--bo', color='blue')

# unit test
# for i in range(len(param.time)):
#     for j in range(len(r[i])):
#         if r[i][j]-r_grid[i][j] < 0 or r_grid[i][j+1]-r[i][j]<0:
#             print(j, 'ERRORRRRR')
#             for k in range(6, 0, -1):
#                 print(j-k, r_grid[i][j-k])
#             print(j, r_grid[i][j], r[i][j], r_grid[i][j+1])
#             for k in range(1, 6):
#                 print(j+k, r_grid[i][j+k])


# for i in range(len(param.time)):
#     print(r_grid[i][-1]/param.AU, r[i][-1]/param.AU, tau_2[i])

# c = cm81_envelope(r=r)
# c.draw()

# plt.yscale('log')
# plt.xscale('log')

# plt.show()