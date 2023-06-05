from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# equation of motion given by equation 11 and 12
def ode(param, x):
    a = param[0]
    v = param[1]

    dvdx = 1/((x-v)**2 -1) * (x-v) * (a*(x-v) - 2/x)
    dadx = a/((x-v)**2 -1) * (x-v) * (a - 2/x*(x-v))

    return [dadx, dvdx]

# number of data points sampled
fine = 10000
# integrate backward from x=1 (divide it back at the end, so that x=1)
x = [i for i in range(fine, 1, -1)]

# initial condition, given by table 2, solution for x=1
p0 = [2, 0]

sol = odeint(ode, p0, x)

a = [i for i in sol[:, 0]]
v = [-i for i in sol[:, 1]]

#'''
AU = 1.496e+13                      # physical constant, Astronomical unit in (cm)
cs = 0.2*1e5                        # condition given by paper (cm/s)
time = [1e12, 2e12, 4e12, 8e12]     # time tested given by the paper (s)

# change the unit from dimensionless solution (a, v) back to normal variable
for t in time:
    r = [i*cs*t/fine for i in x]
    u = [(cs*i/fine) for i in v]
    rho = [(i/(4*np.pi*6.67e-8*t**2)) for i in a]

    plt.plot(r, u, label=f'u, {t/1e12} e12')
    #plt.plot(r, rho, label=f'rho, {t/1e12}e12')
#'''

# apply if want to graph dimensionless (a, v)
'''
plt.plot([i/fine for i in x], a, color='red', label='a')
plt.plot([i/fine for i in x], v, color='blue', label='v')
#'''

#print(a)
plt.legend()
plt.yscale('log')
plt.xscale('log')

plt.xlabel(r'$\log r, (cm)$')
plt.ylabel(r'$\log |u|, (cgs)$')
#plt.ylabel(r'$\log \rho, (cgs)$')

plt.show()