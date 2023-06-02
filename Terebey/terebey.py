# correct density but wrong w

import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

def quadrupolar_greater1_ode(param, z):
    a = param[0]
    v = param[1]
    w = param[2]
    q = param[3]
    p = param[4]

    dadz = 1/(z**2-1) * ((-12*z**2*w) + 2*z**2*(-a/z + 2*v -3*q*z**4 + 2*p/z))
    dvdz = 1/(z**2-1) * (1/z*(-a/z + 2*v - 3*q*z**4 + 2*p/z) - 6*z*w)
    dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
    dqdz = -5*a/z**6
    dpdz = a/5/z

    if (dadz <0):
        dadz = 0

    return [dadz, dvdz, dwdz, dqdz, dpdz]

def quadrupolar_smaller1_ode(param, z):
    a = param[0]
    v = param[1]
    w = param[2]
    q = param[3]
    p = param[4]

    # this is prob wrong, 
    a0 = 71.5
    v0 = 5.44

    dadz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*z*(a*v0+v*a0) + a-2*v - 6*z*a0*w) + 2*a/a0 + 3*a0*v - 3*q*a0*z**4 + 2*p*a0/z)
    dvdz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*a/a0**2 + 3*v - z**2 - 3*q*z**4 + 2*p/z) + 2*a*v0*z + a/a0 +2*v*z - 2*v/a0 - 6*w*z)
    dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
    dqdz = -a/5/z**6
    dpdz = -a/5/z

    if (dadz <0):
        dadz = 0

    return [dadz, dvdz, dwdz, dqdz, dpdz]

def monopolar_smaller1_ode(param, z):
    a = param[0]
    v = param[1]
    m = param[2]

    a0 = 71.5
    v0 = 5.44
    m0 = 1/z**2*a0*(1/z-v0)

    dadz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(a*v0*z + a + 2*v*a0*z - v) + 2*a/a0 +3*a0*v + a0*(1/(6*z) + z**2*m - 2/3*(m0/2)**4*z**3))
    dvdz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*a/a0**2 + 3*v + 1/(6*z) + z**2*m -2/3*(m0/2)**4*z**3) + 1/a0*(a*v0*z + a + 2*v*a0*z -2*v))
    dmdz = 1/z**2*(a-1/2)

    if (dadz < 0):
        dadz = 0

    return [dadz, dvdz, dmdz]

fine = 1000

z_g1 = [i/fine  for i in range(1, fine)]
z_s1 = [i for i in range(1, fine)]
z_mono = [i for i in range(1, fine)]


K = 1.5*10**(-3)
param0_g1 = [0, 0, 0, K, 0]

sol_g1 = odeint(quadrupolar_greater1_ode, param0_g1, z_g1)

a_g = sol_g1[:, 0]
v_g = sol_g1[:, 1]
w_g = sol_g1[:, 2] #continuous through x=1

param0_s1 = [a_g[-1]+2/3*v_g[-1], 2/3*v_g[-1], w_g[-1], K*(1+1/35), -2*K/245]

sol_s1 = odeint(quadrupolar_smaller1_ode, param0_s1, z_s1)

a_s = sol_s1[:, 0]
v_s = sol_s1[:, 1]

param0_mono = [1/2, 0, 0]
sol_mono = odeint(monopolar_smaller1_ode, param0_mono, z_mono)

a_m = sol_mono[:, 0]
v_m = sol_mono[:, 1]

plt.plot([1/i for i in z_s1], sol_s1[:, 0], label=r"$\alpha_Q$", color='red')    #x coordinate
#plt.plot([1/i for i in z_s1], [i for i in sol_s1[:, 1]], label=r"$-v_Q$", color='blue')  
#plt.plot([1/i for i in z_s1], sol_s1[:, 2], label=r"$-w_Q$", color='green')

#plt.plot([1/i for i in z_mono], sol_mono[:, 0], label = r"$a_m$")
#plt.plot([1/i for i in z_mono], sol_mono[:, 1], label = r"$-v_m$")
#print([1/i for i in z_s1])

plt.plot([1/i for i in z_g1], sol_g1[:, 0], label=r"$\alpha_Q$", color='red')    #x coordinate
#plt.plot([1/i for i in z_g1], sol_g1[:, 1], label=r"$-v_Q$", color = 'blue')  
#plt.plot([1/i for i in z_g1], sol_g1[:, 2], label=r"$-w_Q$", color = 'green')        

plt.xscale('log')
plt.yscale('log')

plt.xlabel('x')
plt.legend()

plt.show()