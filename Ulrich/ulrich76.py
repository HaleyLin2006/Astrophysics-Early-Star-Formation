import matplotlib.pyplot as plt
import numpy as np

#constants
G = 6.67e-8
M = 1.989e+33
AU = 1.496e+13
a = 0.2*1e5
t = [1e12, 2e12, 4e12, 8e12]
x = [i/20 for i in range(1, 20)]
r = [[i*a*t[0] for i in x], [i*a*t[1] for i in x], [i*a*t[2] for i in x], [i*a*t[3] for i in x]]
angle = [np.pi/6, np.pi/4, np.pi/3, np.pi/2]

################### Ulrich ###################
##### v_r #####
v_r = [[[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))]]
#structure of v_r
#[[[value at different r]=angle]=time]; [time][angle]=[value at different r]
for time in range(len(t)):
    # len(t) = len(r)
    for j in range(len(angle)):
        theta = angle[j]
        for i in r[time]:
            theta_0 = np.arccos(np.roots([50*AU/i, 0, 1-(50*AU/i), -np.cos(theta)])[2])
            v_r[time][j].append(np.sqrt(M*G/i)*np.sqrt(1+np.cos(theta)/np.cos(theta_0)))

marker = ["o", "^", "s", "*"]
for i in range(len(v_r)):
    for time in range(len(t)):
        plt.plot([i/AU for i in r[time]], v_r[time][i], marker[i])
        print(time, angle[i])

plt.yscale('log')
plt.xscale('log')
plt.xlabel("AU", fontsize=20)
plt.ylabel("velocity", fontsize=20)
plt.show()