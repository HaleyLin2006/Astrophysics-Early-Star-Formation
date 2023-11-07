# from parameter import param
# import matplotlib.pyplot as plt

# param = param()

# param.sample_pts = 10
# z_g1_grid = [2/(i+1) for i in np.logspace(param.sample_range, 0, param.sample_pts)]  # also x_s1
# z_s1_grid = np.logspace(0, param.sample_range, param.sample_pts)

# z_g1 = []
# z_s1 = []
# for i in range(param.sample_pts):
#     if i == 0:
#         z_s1_i = 0.99*z_s1_grid[i]
#         z_g1_i = 0.99*z_g1_grid[i]
#         z_s1.append(z_s1_i)
#         z_g1.append(z_g1_i)
#     else:
#         z_s1_i = 0.5*(z_s1_grid[i-1] + z_s1_grid[i])
#         z_g1_i = 0.5*(z_g1_grid[i-1] + z_g1_grid[i])
#         z_s1.append(z_sx1_i)
#         z_g1.append(z_g1_i)

# r_s1 = []
# r_g1 = []
# r_s1_grid = []
# r_g1_grid = []
# for t in param.time:
#     r_s1.append([1/i*param.cs*t for i in z_s1])
#     r_g1.append([1/i*param.cs*t for i in z_g1])
#     r_s1_grid.append([1/i*param.cs*t for i in z_s1_grid])
#     r_g1_grid.append([1/i*param.cs*t for i in z_g1_grid])

# r = []
# r_grid = []
# for i in range(len(param.time)):
#     r.append([j for j in r_s1[i]][::-1]+[j for j in r_g1[i]][::-1])
#     r_grid.append([j for j in r_s1_grid[i]][::-1]+[j for j in r_g1_grid[i]][::-1])

# plt.plot(range(param.sample_pts), [i for i in z_s1], '--bo')
# plt.plot(range(param.sample_pts), [i for i in z_s1_grid], '--bo', color='red')

# print(z_s1)
# print(z_s1_grid)

# # plt.xscale('log')
# plt.yscale('log')
# plt.show()

import numpy as np
import matplotlib.pyplot as plt
from parameter import param

param = param()

param.sample_pts = 1000

z_g1_grid = [2/(i+1) for i in np.logspace(param.sample_range, 0, param.sample_pts)]  # also x_s1
z_s1_grid = np.logspace(0, param.sample_range, param.sample_pts).tolist()

z_g1 = []
z_s1 = []
for i in range(param.sample_pts):
    if i == 0:
        z_s1_i = 0.5*z_s1_grid[i+1]
        z_g1_i = 0.5*z_g1_grid[i+1]
        z_s1.append(z_s1_i)
        z_g1.append(z_g1_i)

        # print(i, z_s1_grid[i], z_s1_i, z_s1_grid[i+1])
        # print(z_g1_grid[i], z_g1_i, z_g1_grid[i+1])
    else:
        z_s1_i = 0.5*(z_s1_grid[i-1] + z_s1_grid[i])
        z_g1_i = 0.5*(z_g1_grid[i-1] + z_g1_grid[i])
        z_s1.append(z_s1_i)
        z_g1.append(z_g1_i)

    # if i%50==0 and i != param.sample_pts-1:
        # print(i, z_s1_grid[i-1], z_s1_i, z_s1_grid[i])
        # print(z_g1_grid[i-1], z_g1_i, z_g1_grid[i])

# for i in range(param.sample_pts):
#     if i%50==0:
#         if i == 0:
#             print(i, 0, z_s1[i], z_s1_grid[0])
#         else:
#             print(i, z_s1_grid[i-1], z_s1[i], z_s1_grid[i])


r = []
for t in param.time:
    r.append([1/i*param.cs*t for i in z_s1[::-1]]+[1/i*param.cs*t for i in z_g1[::-1]])

r_grid = []
for t in param.time:
    r_grid.append([1/i*param.cs*t for i in z_s1_grid[::-1]]+[1/i*param.cs*t for i in z_g1_grid[::-1]])

radius_AU = r

print(len(r[0]))
print(len(radius_AU[2]))

# print(r_grid[0][0])
# print(radius_AU[0][0])

radius_AU_temp = []
r_grid_temp = []

for i in range(len(param.time)):
    radius_AU_temp.append([])
    r_grid_temp.append([])
    for j in range(len(radius_AU[i])):
        # at the end, add an additional interval 
        if j == len(radius_AU[i])-1:
            radius_AU_temp[-1].append(radius_AU[i][j])
            r_grid_temp[-1].append(r_grid[i][j])
            r_grid_temp[-1].append(2*radius_AU[i][j]-r_grid[i][j])

            # print(j, r_grid[i][j], radius_AU[i][j], 2*radius_AU[i][j]-r_grid[i][j])

            # unit test
            # if radius_AU[i][j]-r_grid[i][j] < 0 or 2*radius_AU[i][j]-radius_AU[i][j]<0:
            #     print(j, 'ERRORRRRR')
            #     print(j, r_grid[i][j], radius_AU[i][j], 2*radius_AU[i][j]-r_grid[i][j])

        # avoid the inconsistnecy at x=1
        elif j+1 == param.sample_pts:
            radius_AU_temp[-1].append((r_grid[i][j]+r_grid[i][j+2])/2)
            # print(j, r_grid[i][j+2])
            r_grid_temp[-1].append(r_grid[i][j+1])
            # print(j, r_grid[i][j], (r_grid[i][j]+r_grid[i][j+2])/2, r_grid[i][j+2])

            # unit test
            # if radius_AU[i][j]-r_grid[i][j] < 0 or r_grid[i][j+2]-(r_grid[i][j]+r_grid[i][j+2])/2<0:
            #     print(j, 'ERRORRRRR')

        elif j == param.sample_pts:
            continue

        else:
            radius_AU_temp[-1].append(radius_AU[i][j])
            r_grid_temp[-1].append(r_grid[i][j])

            # print(j, r_grid[i][j], radius_AU[i][j], r_grid[i][j+1])

            # unit test
            # if radius_AU[i][j]-r_grid[i][j] < 0 or r_grid[i][j+1]-radius_AU[i][j]<0:
            #     print(j, 'ERRORRRRR')
            #     print(j, r_grid[i][j], radius_AU[i][j], r_grid[i][j+1])
            #     for k in range(1, 11):
            #         print(j+k, r_grid[i][j+k])
    
radius_AU = radius_AU_temp
r_grid = r_grid_temp

print(len(radius_AU[1]))
print(len(r_grid[1]))


# unit test
for i in range(len(param.time)):
    for j in range(len(radius_AU[i])):
        if radius_AU[i][j]-r_grid[i][j] < 0 or r_grid[i][j+1]-radius_AU[i][j]<0:
            print(j, 'ERRORRRRR')
            for k in range(6, 0, -1):
                print(j-k, r_grid[i][j-k])
            print(j, r_grid[i][j], radius_AU[i][j], r_grid[i][j+1])
            for k in range(1, 6):
                print(j+k, r_grid[i][j+k])


def find_closest_index(arr, target):
    closest_index = 0
    closest_difference = abs(arr[0] - target)

    for i in range(1, len(arr)):
        difference = abs(arr[i] - target)
        if difference < closest_difference:
            closest_index = i
            closest_difference = difference

    return closest_index

# my_array = [1, 4, 7, 10, 13, 16]
# target_value = 9

# closest_index = find_closest_index(my_array, target_value)
# print(f"The index of the closest value to {target_value} is {closest_index}.")

theta_values = np.linspace(param.min_theta, 2*np.pi, param.theta_pts)
theta_values 


# for j in range(len(param.time)):
# plt.plot([i for i in r_grid[0]], '--bo', color='red')
# plt.plot([i for i in r[0]], '--bo', color='blue')

plt.yscale('log')


# plt.show()
