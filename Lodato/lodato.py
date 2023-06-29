import matplotlib.pyplot as plt
import numpy as np
from scipy.special import iv

num_data_pt = 100

############ Figure 4 ############

R0 = 1
r = np.linspace(0, 2., num_data_pt)
tau = [0.01, 0.05, 0.1, 0.15]
sigma = []
m = 1
for t in tau:
    sigma.append([])
    for i in range(len(r)):
        x = r[i]/R0
        sigma_i = m/(np.pi*R0**2) * (x**(-1/4)/t) * np.exp(-(1+x**2)/t) * iv(1/4, 2*x/t)
        sigma[-1].append(sigma_i)

#for i in range(len(sigma)):
#    plt.plot(r, sigma[i], label=tau[i])

############ Figure 5 ############

R1 = 1
r = np.linspace(0, 100, num_data_pt)
time = [1, 2, 4, 7]
sigma = []
C = 10       # normalization constant
# assume viscosity v=r
for t in time:
    sigma.append([])
    for i in range(len(r)):
        sigma_i = C*t**(-3/2)/(3*np.pi*r[i]) * np.exp(-(r[i]/R1)/t)
        sigma[-1].append(sigma_i)

for i in range(len(sigma)):
    plt.plot(r, sigma[i], label=f't={time[i]}')

plt.xscale('log')
plt.yscale('log')

plt.legend()
plt.show()
