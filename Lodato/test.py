import matplotlib.pyplot as plt
import numpy as np
from scipy.special import iv

num_data_pt = 100
r = np.linspace(0, 2., num_data_pt)
R0 = 1

tau = [0.01, 0.05, 0.1, 0.15]
sigma = []
#'''
for t in range(len(tau)):
    sigma.append([])
    for i in range(len(r)):
        x = r[i]/R0
        sigma_i = 2/(np.pi*R0**2) * (x**(-1/4)/t) * np.exp(-(1+x)**2) * iv(1/4, 2*x) * (2*x/t)
        sigma[-1].append(sigma_i) 
#'''
for i in range(len(sigma)):
    plt.plot(r, sigma[i], label=tau[i])

print(sigma[0])

plt.legend()
plt.show()
