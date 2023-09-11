import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from tabulate import tabulate
import numpy as np

#r and rho dependent on beta
r = []
rho = []
beta = []
for b in range(0, 801):
    b /= 100
    beta.append(b)
    rho_i = 3/(3+10*b+7*(b**2))
    r_i = np.sqrt(b*((1+b)**(2/3)))

    r.append(r_i)
    rho.append(rho_i)

    b *= 100

#v dependent on r
v = []
for i in r:
    v.append(i**(1/7))

show = [["r", "v", "rho", "beta"]]
for i in range(len(r)):
    #print("r", r[i], "rho", rho[i])
    show.append([r[i], v[i], rho[i], beta[i]])

print(tabulate(show))

#Graphing
fig = plt.figure()
density = fig.add_subplot(2,1,1)
velocity = fig.add_subplot(2,1,2)

#plt.yscale('log')
#plt.xscale('log')
plt.xlabel(r"$\frac{r}{r_0}$", fontsize=25)

#density graph
density.plot(r, rho)
density.set_ylabel(r"$\frac{\rho}{\rho_0}$", fontsize=25)
for tick in density.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in density.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)

minor_locator = AutoMinorLocator(4)
density.xaxis.set_minor_locator(minor_locator)
density.yaxis.set_minor_locator(minor_locator)

#velocity graph
velocity.plot(r, v)
velocity.set_ylabel(r"$\frac{v}{v_0}$", fontsize=25)
for tick in velocity.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
for tick in velocity.yaxis.get_major_ticks():
    tick.label.set_fontsize(14)

minor_locator = AutoMinorLocator(4)
velocity.xaxis.set_minor_locator(minor_locator)
velocity.yaxis.set_minor_locator(minor_locator)

plt.show()