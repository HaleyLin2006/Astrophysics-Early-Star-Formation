from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

AU = 1.496e+13

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


print("a", alpha)
print("v", v)

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
plt.legend()
plt.show()
