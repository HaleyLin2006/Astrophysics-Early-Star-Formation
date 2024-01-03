import numpy as np
import matplotlib.pyplot as plt

from data_handler import data
from parameter import param
from graphing_help import shu_graph, tsc_graph, cassen_graph

from envelope.envelope_help import envelope
from disk.disk import disk
from star.star_help import star

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

# for t in range(len(param.time)):
#     print("r_z", np.array(r_z[t]))
#     print("r", np.array(r[t]))

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
# t = tsc_graph(z_s1=z_s1, z_g1=z_g1, r_grid=r_grid, r=r, save=True)
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

########## CM (1981) Graph ###########
# cm = cassen_graph(r=r, r_grid=r_grid, save=True)
# cm.v_r_envelope()
# cm.nh_density_envelope_1d()
# cm.density_envelope_2d()

# star = star()
# envelope = envelope()
disk = disk()

################ Star ################
# r_star = star.get_r_star()
# star_mass = envelope.cgs_env_to_star_mass(r_star=r_star)
# l_acc = star.luminosity_accretion(param.time, star_mass)
# plt.plot([i/param.t_acc for i in param.time], [i/(3.846*1e33) for i in l_acc])
# disk_inner_radius = star.disk_inner_radius(star_mass)

############### Envelope ################
# env_mass = envelope.cgs_env_mass()
# plt.plot([i/param.t_acc for i in param.time], [i for i in env_mass])
# env_to_disk_mass = envelope.cgs_env_to_disk_mass(disk_inner_radius=disk_inner_radius)
# plt.plot([i/param.t_acc for i in param.time], [i/param.M for i in env_to_disk_mass])
# env_to_star_mass = envelope.cgs_env_to_star_mass(r_star)
# plt.plot([i/param.t_acc for i in param.time], [i/param.M for i in env_to_star_mass])

############### Disk ################
# disk.disk_evol()

# plt.xscale('log')
# plt.yscale('log')

plt.show()