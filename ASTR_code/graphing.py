import numpy as np
import matplotlib.pyplot as plt

from data_handler import data
from parameter import param
from graphing_help import shu_graph, tsc_graph

param = param()
data = data()

# make grid
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
z_g1[-1] =1-1e-10

print(z_s1_grid[0], z_s1_grid[-1], z_g1_grid[0], z_g1_grid[-1])
print(z_s1[0], z_s1[-1], z_g1[0], z_g1[-1])

r = []
for t in param.time:
    r.append([1/i*param.cs*t for i in z_s1[::-1]]+[1/i*param.cs*t for i in z_g1[::-1]])

z = []
for t in param.time:
    z.append(z_s1[::-1]+ z_g1[::-1])

r_grid = []
for t in param.time:
    r_grid.append([1/i*param.cs*t for i in z_s1_grid[::-1]]+[1/i*param.cs*t for i in z_g1_grid[::-1]])

r_temp = []
r_grid_temp = []
for i in range(len(param.time)):
    r_temp.append([])
    r_grid_temp.append([])
    for j in range(len(r[i])):
        # at the end, add an additional interval 
        if j == len(r[i])-1:
            r_temp[-1].append(r[i][j])
            r_grid_temp[-1].append(r_grid[i][j])
            r_grid_temp[-1].append(2*r[i][j]-r_grid[i][j])

        # avoid the inconsistnecy at x=1
        elif j+1 == param.sample_pts:
            r_temp[-1].append((r_grid[i][j]+r_grid[i][j+2])/2)
            r_grid_temp[-1].append(r_grid[i][j+1])
        elif j == param.sample_pts:
            continue

        else:
            r_temp[-1].append(r[i][j])
            r_grid_temp[-1].append(r_grid[i][j])

r = r_temp
r_grid = r_grid_temp

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

plt.plot(r_grid[0],'--bo', color='red')
plt.plot(r[0], '--bo')

########## Shu (1977) Graph ##########
# shu = shu_graph(z_s1=z_s1, z_g1=z_g1, r_grid=r_grid, save=False)
# shu.dimensionless_a_v_1d()
# shu.cgs_density_1d()
# shu.nh_density_1d()
# shu.cgs_density_2d()
# shu.nh_density_2d()

########## TSC (1984) Graph ##########
# t = tsc_graph(z_s1=z_s1, z_g1=z_g1, r_grid=r_grid, save=False)
# t.dimensionless_a_1d()
# t.dimensionless_quad_v_w_1d()
# t.dimensionless_mono_v_1d()
# t.cgs_density_1d()
# t.nh_density_1d()
# t.cgs_density_2d()
# t.nh_density_2d()
# t.int_mass()

# plt.xscale('log')
plt.yscale('log')

plt.show()