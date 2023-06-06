import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint, ode

########################## Equation of Motions ##########################

#Shu 1977 EOM for a0, v0 used in quadrupolar x<1 and monopolar solution
store_dadx = []
store_dvdx = []
def shu_ode(param, x):
    a = param[0]
    v = param[1]

    dvdx = 1/((x-v)**2 -1) * (x-v) * (a*(x-v) - 2/x)
    dadx = a/((x-v)**2 -1) * (x-v) * (a - 2/x*(x-v))

    # extract dadx and dvdx for future usage (x<1)
    global store_dadx
    store_dadx.append(dadx)
    global store_dvdx
    store_dvdx.append(dvdx)
    
    return [dadx, dvdx] 

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
# adapting Shu (1976) solution for v_0, a_0 needed for the solution
def quadrupolar_smaller1_ode(z, param, a0, v0):
    a = param[0]
    v = param[1]
    w = param[2]
    q = param[3]
    p = param[4]

    '''
    A_Q = a*z**2*(2/z - dv0dz) + v*z**2*(2/z - da0dz) - 6*z*a0*w
    B_Q = a/(a0**2)*(z**2)*da0dz + (2+dv0dz)*v - 3*q*z**4 + 2*p/z

    dadz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*A_Q + a0*B_Q)
    
    v0 = -v0

    dvdz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*B_Q + 1/a0*A_Q)
    dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
    dqdz = -a/5/z**6
    dpdz = -a/5/z

    '''
    # starting from equation 69, 1. change to z=1/x parameter, 2. adapt boundary conditions from equation 26(b), 67(b), and 68
    dadz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*z*(a*v0+v*a0) + a-2*v - 6*z*a0*w) + 2*a/a0 + 3*a0*v - 3*q*a0*z**4 + 2*p*a0/z)
    
    # cheating?
    v0 = -v0

    dvdz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*a/a0**2 + 3*v - z**2 - 3*q*z**4 + 2*p/z) + 2*a*v0*z + a/a0 +2*v*z - 2*v/a0 - 6*w*z)
    dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
    dqdz = -a/5/z**6
    dpdz = -a/5/z
    #'''

    # physical constraint, dadz can never be negative, i.e. you can never take away density
    if (dadz < 0):
        dadz = 0

    return [dadz, dvdz, dwdz, dqdz, dpdz]

def monopolar_smaller1_ode(z, param, a0, v0):
    a = param[0]
    v = param[1]
    m = param[2]

    m0 = 1/z**2*a0*(1/z-v0)
    dpsidz = -m - 1/(6*z**3)

    # 1. change to z=1/x parameter, 2. following equation 73 and 74, following constraints of equation 78, 79
    dadz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*a*v0*z + a + 2*v*a0*z - 2*v) + 2*a/a0 +3*a0*v + a0*(-z**2*dpsidz - 2/3*(m0/2)**4*z**3))
    dvdz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*a/a0**2 + 3*v -z**2*dpsidz - 2/3*(m0/2)**4*z**3) + 1/a0*(2*a*v0*z + a + 2*v*a0*z - 2*v))
    dmdz = 1/z**2*(a-1/2)

    # cheating?    
    if (dadz < 0):
        dadz = -dadz
    
    return [dadz, dvdz, dmdz]

############# Numerical Solutions to the Equation of Motions #############

# number of data points sampled
fine = 1000

# the "scale" for the following calculation, in z=1/x coordinate
z_shu = [i for i in range(fine, 1, -1)]
z_g1 = [i/fine  for i in range(1, fine)]
z_s1 = [i for i in range(1, fine)]
z_mono = [i for i in range(1, fine)]

# boundary conditions for Shu 1977
p0 = [2, 0]

#solution for Shu 1977
sol_shu = odeint(shu_ode, p0, z_shu)
a_shu = [i for i in sol_shu[:, 0]]
v_shu = [-i for i in sol_shu[:, 1]]

print(store_dvdx[0])

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
f_s1.set_f_params(a_shu[0], v_shu[0])

for i in range(1, len(z_s1)):
    sol_s1[i, :] = f_s1.integrate(z_s1[i])
    f_s1.set_initial_value(sol_s1[i, :], z_s1[i])
    f_s1.set_f_params(a_shu[i], v_shu[i])


# boundary condition for monopolar, x<1 model
# adopting equation 78 and 79. a_shu, v_shu from Shu(1976)
# a, v, m
param0_mono = [1/2, 0, 0]

#'''
# solution for monopolar, x<1 model
sol_m = np.zeros((len(z_mono), len(param0_mono)))
sol_m[0, :] = param0_mono
f_m = ode(monopolar_smaller1_ode)
f_m.set_initial_value(param0_mono, z_mono[0])
f_m.set_f_params(a_shu[0], v_shu[0])

for i in range(1, len(z_mono)):
    sol_m[i, :] = f_m.integrate(z_mono[i])
    f_m.set_initial_value(sol_m[i, :], z_mono[i])
    f_m.set_f_params(a_shu[i], v_shu[i])
#'''

######################## Plotting the Solutions ########################

# Shu 1977
plt.plot([i/fine for i in z_shu], sol_shu[:, 0], label=r"$\alpha_0$")                          # alpha_0
#plt.plot([i/fine for i in z_shu], [-i for i in sol_shu[:, 1]], label=r"$-v_0$")                # v_0

# quadrupolar, x>1 model
# 1/i to change z=1/x back to x coordinate
plt.plot([1/i for i in z_g1], sol_g1[:, 0], label=r"$\alpha_Q$", color='red')                  # alpha_Q 
#plt.plot([1/i for i in z_g1], sol_g1[:, 1], label=r"$v_Q$", color= 'blue')                     # -v_Q
#plt.plot([1/i for i in z_g1], sol_g1[:, 2], label=r"$w_Q$", color='green')                     # w_Q


# quandrupolar, x<1 model
plt.plot([1/i for i in z_s1], sol_s1[:, 0], label=r"$\alpha_Q$", color='red')                 # alpha_Q
#plt.plot([1/i for i in z_s1], [-i for i in sol_s1[:, 1]], label=r"$-v_Q$", color='blue')       # -v_Q
#plt.plot([1/i for i in z_s1], sol_s1[:, 2], label=r"$w_Q$", color='green')                     # w_Q

#monopolar is incorrect
# monopolar, x<1 model
#plt.plot([1/i for i in z_mono], sol_m[:, 0], label = r"$a_m$")
#plt.plot([1/i for i in z_mono], sol_m[:, 1], label = r"$v_m$")


plt.xscale('log')
plt.yscale('log')

plt.xlabel('x')
plt.legend()

plt.show()


# problems:
# need to adopt dv0dx, da0dx x-dependent --> maybe change shu to z coordinate