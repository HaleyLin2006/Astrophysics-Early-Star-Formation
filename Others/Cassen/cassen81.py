from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np

time = [0.1, 0.5, 1]
M = 1.989e+33   # 1 solar mass

def s(x):
    return ((1-x**(1/2)) * (3-2*x)) / (4*x*(1-x)**(3/2)) * x**(3/4)

u = np.logspace(-3, 0, 1000, endpoint=False)

sigma = []
for t in time:
    sigma.append([])
    print(t)
    for x in u:
        sigma_i = M**2/(4 * t**7 * x**(7/4))
        s_int = integrate.quad(s, x, 1)
        normalized = 16*x*t**7 / M**2
        sigma[-1].append(normalized * sigma_i * s_int[0])

color = ['blue', 'red','green']
for i in range(len(sigma)):
    plt.plot(u, sigma[i], color=color[i])


plt.yscale('log')
plt.show()