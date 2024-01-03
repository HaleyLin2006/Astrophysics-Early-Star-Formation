import numpy as np
import os
from data_handler import data
from parameter import param
import matplotlib.pyplot as plt

param = param()
data = data()

# r_c    = [param.cs/param.omega /param.AU for t in param.time]
# r_d    = [(param.omega**2*param.G**3*param.M**3) / (16*param.cs**8) / param.AU for t in param.time]

# t_pts = 10
# t_range = -3
# time = [i*param.t_acc for i in np.logspace(t_range, 0, t_pts)]
# print(time)

param.time = [i*param.t_acc for i in np.linspace(0,1)]
param.time.pop(0)
print(np.linspace(0,1))

centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 / param.AU for t in param.time]
print(centrifugal_radius)
tau_2 = [(param.omega*t)**2* param.cs*t /param.AU for t in param.time]
print(tau_2)

plt.plot(centrifugal_radius, color='blue', label='centrifugal radius')
plt.plot(tau_2, color='red', label='tau^2')

plt.ylabel('AU')
plt.xlabel('time')
plt.yscale('log')

plt.legend()
plt.show()
# data_param = ["Shu", param.cs, param.omega, len(param.time), param.sample_pts, param.theta_pts]

# data.delete(data_param)

# data.save(param=data_param, radius=[[i for i in range(10)] for j in range(3)], density=[[[i for i in range(10)] for j in range(3)] for k in range(5)], v_r=[], v_theta=[], v_phi=[])

# radius, r_grid, density, v_r, v_theta, v_phi = data.load(param=data_param)

# print('raidus', radius)
# print('density', density)
# print('v_r', v_r)
# print('v_theta', v_theta)
# print('v_phi', v_phi)