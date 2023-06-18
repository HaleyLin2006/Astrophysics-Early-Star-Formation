# sample more for x = 1+
# adopt da0dz, dv0dz
# v_shu_s1 is wrong (approx 1e3 magnitude), causing v_q to be not smooth

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode, odeint

########################## Equation of Motions ##########################

#Shu 1977 EOM for a0, v0 used in quadrupolar x<1 and monopolar x<1 solution
store_dadz = []
store_dvdz = []
def shu_ode(z, param):
    a = param[0]
    v = param[1]

    dadz = -a/((1-z*v)**2 - z**2) * (1/z-v) * (a - 2*z*(1/z-v))
    dvdz = -1/((1-z*v)**2 - z**2) * (1/z-v) * (a*(1/z-v) - 2*z)

    if dadz<0:
        dadz = 0

    global store_dadz
    store_dadz.append(dadz)
    global store_dvdz
    store_dvdz.append(dvdz)

    # progress recorder
    print("shu, z:", z)

    return [dadz, dvdz]

# quadrupolar EOM for x>1 
def quadrupolar_greater1_ode(param, z):
    a = param[0]
    v = param[1]
    w = param[2]
    q = param[3]
    p = param[4]

    # equation 69(a)-(e), changing to z=1/x parameter as our boundary condition works for $x\to\infty$
    dadz = 1/(z**2-1) * ((-12*z**2*w) + 2*z**2*(-a/z + 2*v -3*q*z**4 + 2*p/z))
    dvdz = 1/(z**2-1) * (1/z*(-a/z + 2*v - 3*q*z**4 + 2*p/z) - 6*z*w)
    dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
    dqdz = -5*a/z**6
    dpdz = a/5/z

    # physical constraint, dadz can never be negative, i.e. you can never take away density
    if (dadz <0):
        dadz = 0

    return [dadz, dvdz, dwdz, dqdz, dpdz]

# quadrupolar EOM for x<1
def quadrupolar_smaller1_ode(z, param, a0, v0):
    a = param[0]
    v = param[1]
    w = param[2]
    q = param[3]
    p = param[4]

    # starting from equation 69, 1. change to z=1/x parameter, 2. adapt boundary conditions from equation 26(b), 67(b), and 68
    dadz = -1/((1-z*v0)**2 - z**2) * ((1/z-v0)*(2*z*(a*v0+v*a0) + a-2*v - 6*z*a0*w) + 2*a/a0 + 3*a0*v - 3*q*a0*z**4 + 2*p*a0/z)

    dvdz = -1/((1-z*v0)**2 - z**2) * ((1/z-v0)*(2*a/a0**2 + 3*v - z**2 - 3*q*z**4 + 2*p/z) + 2*a*v0*z + a/a0 +2*v*z - 2*v/a0 - 6*w*z)
    dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
    dqdz = -a/5/z**6
    dpdz = -a/5/z
    #'''

    # physical constraint, dadz can never be negative, i.e. you can never take away density
    if (dadz < 0):
        dadz = 0

    # progress recorder
    print("quadruploar x<1, z:", z)

    return [dadz, dvdz, dwdz, dqdz, dpdz]

############# Numerical Solutions to the Equation of Motions #############
sample_pts = 1000
sample_range = 3        # in interval [1e(-sample_range), 1e(sample_range)]

# z=1/x
z_g1 = np.logspace(-sample_range, 0, sample_pts)
z_s1 = np.logspace(0, sample_range, sample_pts)

#'''
# initial condition for Shu, z -> 0 (x -> infty)
A = 2.0001
param0_g1 = [A*z_g1[0]**2, -(A-2)*z_g1[0]]

# Solution for Shu, z<1, (x>1)
sol_shu_g1 = np.zeros((len(z_g1), len(param0_g1)))
sol_shu_g1[0, :] = param0_g1
f_shu_g1 = ode(shu_ode)
f_shu_g1.set_integrator("vode", method = "bdf", atol=1e-20, nsteps=1000)
f_shu_g1.set_initial_value(param0_g1, z_g1[0])

for i in range(1, len(z_g1)):
    sol_shu_g1[i, :] = f_shu_g1.integrate(z_g1[i])
    f_shu_g1.set_initial_value(sol_shu_g1[i, :], z_g1[i])

a_shu_g1 = [i for i in sol_shu_g1[:, 0]]
v_shu_g1 = [i for i in sol_shu_g1[:, 1]]

# initial condition for Shu, z>1, (x<1)
param0_s1 = [a_shu_g1[-1], v_shu_g1[-1]]

# Solution for Shu, z>1, (x<1)
sol_shu_s1 = np.zeros((len(z_s1), len(param0_s1)))
sol_shu_s1[0, :] = param0_s1
f_shu_s1 = ode(shu_ode)
f_shu_s1.set_integrator("vode", method="bdf", atol=1e-20, nsteps=1000)
f_shu_s1.set_initial_value(param0_s1, z_s1[0])

for i in range(1, len(z_s1)):
    sol_shu_s1[i, :] = f_shu_s1.integrate(z_s1[i])
    f_shu_s1.set_initial_value(sol_shu_s1[i, :], z_s1[i])

a_shu_s1 = [i for i in sol_shu_s1[:, 0]]
v_shu_s1 = [i for i in sol_shu_s1[:, 1]]

# boundary condition for quandrupolar, x>1 model
# adopting equation 65 and 70, K given in III-c of the paper
K = 1.5*10**(-3)
# a, v, w, q, p
param0_g1 = [0, 0, 0, K, 0]

# solution for quadrupolar, x>1 model
sol_g1 = odeint(quadrupolar_greater1_ode, param0_g1, z_g1)
a_g = sol_g1[:, 0]
v_g = sol_g1[:, 1]
w_g = sol_g1[:, 2] #continuous through x=1

# boaundary condition for quadrupolar, x<1 model
# adopting equation 67 and 68 to go over the x=1 discontinuity. a_shu and v_shu from Shu(1976)
# here, we manually take limit by sampling points close to 1, such that a(1+e) = a(1) = a(1-e), where e is an arbitrary small number
# a, v, w, q, p
param0_s1 = [a_g[-1]+2/3*v_g[-1], 2/3*v_g[-1], w_g[-1], K*(1+1/35), -2*K/245]

# solution for quadrupolar, x<1 model
sol_s1 = np.zeros((len(z_s1), len(param0_s1)))
sol_s1[0, :] = param0_s1
f_s1 = ode(quadrupolar_smaller1_ode)
f_s1.set_initial_value(param0_s1, z_s1[0])
f_s1.set_integrator("vode", method = "bdf", atol=1e-20, nsteps=1000)
f_s1.set_f_params(a_shu_s1[0], v_shu_s1[0])

for i in range(1, len(z_s1)):
    sol_s1[i, :] = f_s1.integrate(z_s1[i])
    f_s1.set_initial_value(sol_s1[i, :], z_s1[i])
    #print(i)
    f_s1.set_f_params(a_shu_s1[i], v_shu_s1[i])         # need to plug in dadz, dvdz

######################## Plotting the Solutions ########################
# Shu
plt.plot([1/i for i in z_g1], [i for i in a_shu_g1], color='red', label='a')
plt.plot([1/i for i in z_g1], [-i for i in v_shu_g1], color='orange', label='-v shu')

plt.plot([1/i for i in z_s1], [i for i in sol_shu_s1[:, 0]], color='red', label='a')
plt.plot([1/i for i in z_s1], [-i for i in sol_shu_s1[:, 1]], color='orange', label='-v shu')      #def cheating

# quadrupolar, x>1 model
# 1/i to change z=1/x back to x coordinate
#plt.plot([1/i for i in z_g1], sol_g1[:, 0], label=r"$\alpha_Q$", color='red')                  # alpha_Q 
plt.plot([1/i for i in z_g1], sol_g1[:, 1], label=r"$v_Q$", color= 'blue')                     # -v_Q
plt.plot([1/i for i in z_g1], sol_g1[:, 2], label=r"$w_Q$", color='green')                     # w_Q

# quandrupolar, x<1 model
#plt.plot([1/i for i in z_s1], sol_s1[:, 0], label=r"$\alpha_Q$", color='red')                 # alpha_Q
plt.plot([1/i for i in z_s1], [i for i in sol_s1[:, 1]], label=r"$-v_Q$", color='blue')       # -v_Q
plt.plot([1/i for i in z_s1], [i for i in sol_s1[:, 2]], label=r"$w_Q$", color='green')         # w_Q

#monopolar is incorrect
# monopolar, x<1 model
#plt.plot([1/i for i in z_mono], sol_m[:, 0], label = r"$a_m$")
#plt.plot([1/i for i in z_mono], sol_m[:, 1], label = r"$v_m$")

plt.legend()
plt.yscale('log')
plt.xscale('log')

plt.show()