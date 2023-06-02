from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

def ode(param, x):
    a = param[0]
    v = param[1]

    dvdx = 1/((x-v)**2 -1) * (x-v) * (a*(x-v) - 2/x)
    dadx = a/((x-v)**2 -1) * (x-v) * (a - 2/x*(x-v))

    
    return [dadx, dvdx]

fine = 10 
x = [i for i in range(fine, 1, -1)]
p0 = [2, 0]

sol = odeint(ode, p0, x)

a = [i for i in sol[:, 0]]
v = [-i for i in sol[:, 1]]

'''
AU = 1.496e+13
cs = 0.2*1e5
time = [1e12, 2e12, 4e12, 8e12]
count = 1
for t in time:
    r = [i*cs*t for i in z]
    u = [(cs*i/1) for i in v]
    rho = [(i/(4*np.pi*6.67e-8*t**2)) for i in a]

    #plt.plot(r, u, label=f'v{count}')
    plt.plot(r, rho, label=f'rho{count}')
    count += 1
'''

#plt.plot(r, rho, label='rho')

#print(v)

plt.plot([i/fine for i in x], a, color='red', label='a')
plt.plot([i/fine for i in x], v, color='blue', label='v')

#print(a)
plt.legend()
plt.yscale('log')
plt.xscale('log')
plt.show()