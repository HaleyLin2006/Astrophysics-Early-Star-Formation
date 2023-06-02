import matplotlib.pyplot as plt
from scipy.integrate import odeint
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


'''
##### density #####
density = [[[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))]]

#theta; i
for time in range(len(t)):
    for j in range(len(angle)):
        theta = angle[j]
        for i in r[time]:
            theta_0 = np.arccos(np.roots([50*AU/i, 0, 1-(50*AU/i), -np.cos(theta)])[2])
            # dm for ulrich
            #dm = 1e-7*M/3.154e+7*np.exp(-np.cos(theta_0)/2) 

            #dm for shu
            dm = 2e-6*M/3.154e7 

            #print("dm=", dm, "theta=", theta, "theta_0=", theta_0)
            d = dm/(8*np.pi*np.sqrt(G*M*i**3))
            d *= (1+np.cos(theta)/np.cos(theta_0))**(-1/2)
            d *= (np.cos(theta)/(2*np.cos(theta_0))+(50*AU/i)*np.cos(theta_0)**2)**(-1)

            #d = (dm/(8*np.pi*np.sqrt(G*M*i**3)))*np.sqrt(np.cos(theta_0)/(np.cos(theta)+np.cos(theta_0)))/(np.cos(theta)/(2*np.cos(theta_0))+50*AU/i*np.cos(theta_0)**2)
            density[time][j].append(d)

marker = ["o", "^", "s", "*"]
#for i in range(len(density)):
i = 1
for time in range(len(t)):
    angle_name = [30, 45, 60, 90]
    color = ['red', 'blue', 'green', 'purple']
    plt.plot([i/AU for i in r[time]], density[time][i], marker[time], label=f'{t[time]/1e12}e12 {angle_name[i]}', color=color[time])


plt.yscale('log')
plt.xscale('log')
plt.xlabel("AU", fontsize=20)
plt.ylabel("density", fontsize=20)
plt.legend()
plt.show()
'''
'''
################### Shu ###################
def ode(p, x):
    
    v = p[0]
    a = p[1]

    dvdx = (a*(x-v)-(2/x))*(x-v)/((x-v)**2-1)
    dadx = (a*(a-(2/x)*(x-v))*(x-v))/((x-v)**2-1)

    return [dvdx, dadx]
#init condition
p0 = [-5.44, 71.5]

a = 0.2*1e5
t = [1e12, 2e12, 4e12, 8e12]

x = [i/20 for i in range(1, 20)]
sol = odeint(ode, p0, x)

r = [[i*a*t[0] for i in x], [i*a*t[1] for i in x], [i*a*t[2] for i in x], [i*a*t[3] for i in x]]
v = sol[:, 0]
u = [-(a*i) for i in v]
alpha = sol[:, 1]
rho = [[(i/(4*np.pi*6.67e-8*t[0]**2)) for i in alpha], [(i/(4*np.pi*6.67e-8*t[1]**2)) for i in alpha], [(i/(4*np.pi*6.67e-8*t[2]**2)) for i in alpha], [(i/(4*np.pi*6.67e-8*t[3]**2)) for i in alpha]]

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

plt.xlabel(r"$\log(r)$ AU", fontsize=15)
plt.ylabel(r"$\log(\rho)$", fontsize=15)
#plt.ylabel(r"$\log(\mid u \mid)$", fontsize=15)
for time in range(len(t)):
    plt.plot([i/AU for i in r[time]], zero_to_nan(rho[time]), label=f"{t[time]/10**12}e12")
    #plt.plot([i/AU for i in r[time]], u, label=f"{t[time]/10**12}e12")

plt.yscale('log')
plt.xscale('log')
plt.xlabel("AU", fontsize=20)
plt.ylabel("density", fontsize=20)
plt.legend()
plt.show()
'''