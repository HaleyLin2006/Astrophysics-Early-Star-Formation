import matplotlib.pyplot as plt
import numpy as np
from scipy.special import iv

num_data_pt = 100

# change time, need to manually integrate i.e. R1 is vary and note that the scaling factor need to accoundt for different density

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

R1 = 1000
r = np.linspace(0, 1000, num_data_pt)
time = [1, 2, 4, 7]
sigma = []
C = 10       # normalization constant
# assume viscosity v=r
for t in time:
    sigma.append([])
    for i in range(len(r)):
        sigma_i = C*t**(-3/2)/(3*np.pi*r[i]) * np.exp(-(r[i]/R1)/t)
        sigma[-1].append(sigma_i)

#'''
for i in range(len(sigma)):
    plt.plot(r, sigma[i], label=f't={time[i]}')
#'''

############ Figure 5 ############

# need to check time
# change C

'''
init_r = [1, 10, 100, 1000]
color = ['red', 'blue', 'green', 'purple']
r = np.linspace(0, 1000, num_data_pt)
time = [1, 2, 4, 7]
C = 1000       # normalization constant
# assume viscosity v=r
sigma_sum = r
count = 0
for r1 in init_r:
    for t in time:
        sigma = []
        for i in range(len(r)):
            sigma_i = C*t**(-3/2)/(3*np.pi*r[i]) * np.exp(-(r[i]/r1)/t)
            sigma.append(sigma_i)
        plt.plot(r, sigma, color=color[count])
    count += 1
#'''

# Second Graph
#for i in range(len(sigma)):
#    plt.plot(r, sigma[i], label=f't={time[i]}')

plt.xlabel(r'$R$')
plt.ylabel(r'$\Sigma$')

plt.xscale('log')
plt.yscale('log')

plt.legend()
plt.show()
