# correct w but wrong density

import matplotlib.pyplot as plt
from scipy.integrate import odeint

########################## Equation of Motions ##########################

# quadrupolar solution for x>1 
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

# quadrupolar solution for x<1
def quadrupolar_smaller1_ode(param, z, test):
    a = param[0]
    v = param[1]
    w = param[2]
    q = param[3]
    p = param[4]

    # adapting Shu (1976) solution for v_0, a_0 needed for the solution
    a_shu = param[5]
    v_shu = param[6]

    #'''
    a0 = a_shu
    v0 = v_shu
    #'''
    
    #a0 = 71.5
    #v0 = 5.44

    # starting from equation 69, 1. change to z=1/x parameter, 2. adapt boundary conditions from equation 26(b), 67(b), and 68
    dadz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*z*(a*v0+v*a0) + a-2*v - 6*z*a0*w) + 2*a/a0 + 3*a0*v - 3*q*a0*z**4 + 2*p*a0/z)
    dvdz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*a/a0**2 + 3*v - z**2 - 3*q*z**4 + 2*p/z) + 2*a*v0*z + a/a0 +2*v*z - 2*v/a0 - 6*w*z)
    dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
    dqdz = -a/5/z**6
    dpdz = -a/5/z

    #da_shudz = a_shu/((z-v_shu)**2 -1) * (z-v_shu) * (a_shu - 2/z*(z-v_shu))
    #dv_shudz = 1/((z-v_shu)**2 -1) * (z-v_shu) * (a_shu*(z-v_shu) - 2/z)

    da_shudz = a0/((z*v0-1)**2 +z**2) * (1/z-v0) * (a0 - 2*z*(1/z-v0))
    dv_shudz = 1/((z*v0-1)**2 +z**2) * (1/z-v0) * (a0*(1/z-v0) - 2*z)

    # physical constraint, dadz can never be negative, i.e. you can never take away density
    if (dadz <0):
        dadz = 0

    return [dadz, dvdz, dwdz, dqdz, dpdz, da_shudz, dv_shudz]

def monopolar_smaller1_ode(param, z):
    a = param[0]
    v = param[1]
    m = param[2]

    # adapting Shu (1976) solution for v_0, a_0 needed for the solution
    a_shu = param[3]
    v_shu = param[4]

    a0 = a_shu
    v0 = v_shu
    m0 = 1/z**2*a0*(1/z-v0)

    # 1. change to z=1/x parameter, 2. following equation 73 and 74, following constraints of equation 78, 79
    dadz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(a*v0*z + a + 2*v*a0*z - 2*v) + 2*a/a0 +3*a0*v + a0*(1/(6*z) + z**2*m - 2/3*(m0/2)**4*z**3))
    dvdz = 1/((z*v0-1)**2 + z**2) * ((1/z-v0)*(2*a/a0**2 + 3*v + 1/(6*z) + z**2*m -2/3*(m0/2)**4*z**3) + 1/a0*(a*v0*z + a + 2*v*a0*z -2*v))
    dmdz = 1/z**2*(a-1/2)

    da_shudz = -a_shu/((z*v_shu-1)**2 + z**2) * ((1/z-v_shu)*(a_shu - 2*z*(1/z - v_shu)))
    dv_shudz = 1/((z*v_shu-1)**2 + z**2) * ((1/z-v_shu)*(a_shu*(1/z - v_shu) - 2*z))
    
    return [dadz, dvdz, dmdz, da_shudz, dv_shudz]

############# Numerical Solutions to the Equation of Motions #############

#numbers of data points sampled
fine = 1000

# the "scale" for the following calculation, in z=1/x coordinate
z_g1 = [i/fine  for i in range(1, fine)]
z_s1 = [i for i in range(1, fine)]
z_mono = [i for i in range(1, fine)]

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
# a, v, w, q, p, a_shu, v_shu
param0_s1 = [a_g[-1]+2/3*v_g[-1], 2/3*v_g[-1], w_g[-1], K*(1+1/35), -2*K/245, 2, 0]

# solution for quadrupolar, x<1 model
sol_s1 = odeint(quadrupolar_smaller1_ode, param0_s1, z_s1, ("cool", ))

a_s = sol_s1[:, 0]
v_s = sol_s1[:, 1]

# boundary condition for monopolar, x<1 model
# adopting equation 78 and 79. a_shu, v_shu from Shu(1976)
# a, v, m, a_shu, v_shu
param0_mono = [1/2, 0, 0, 2, 0]

# solution for monopolar, x<1 model
sol_mono = odeint(monopolar_smaller1_ode, param0_mono, z_mono)

a_m = sol_mono[:, 0]
v_m = sol_mono[:, 1]

######################## Plotting the Solutions ########################

# quadrupolar, x>1 model
# 1/i to change z=1/x back to x coordinate
#plt.plot([1/i for i in z_g1], sol_g1[:, 0], label=r"$\alpha_Q$", color='red')
plt.plot([1/i for i in z_g1], sol_g1[:, 1], label=r"$v_Q$", color= 'blue')  
plt.plot([1/i for i in z_g1], sol_g1[:, 2], label=r"$w_Q$", color='green')


# quandrupolar, x<1 model
#plt.plot([1/i for i in z_s1], sol_s1[:, 0], label=r"$\alpha_Q$", color='red')    
plt.plot([1/i for i in z_s1], [-i for i in sol_s1[:, 1]], label=r"$-v_Q$", color='blue')  
plt.plot([1/i for i in z_s1], sol_s1[:, 2], label=r"$w_Q$", color='green')
#plt.plot([1/i for i in z_s1], sol_s1[:, 5], label=r"$a_0$", color='orange')
#plt.plot([1/i for i in z_s1], [i for i in sol_s1[:, 6]], label=r"$v_0$", color='yellow')

# monopolar, x<1 model
#plt.plot([1/i for i in z_mono], sol_mono[:, 0], label = r"$a_m$")
#plt.plot([1/i for i in z_mono], sol_mono[:, 1], label = r"$v_m$")
#plt.plot([1/i for i in z_mono], [-i for i in sol_mono[:, 4]], label = r"$-v_0$") # v shu
#plt.plot([1/i for i in z_mono], sol_mono[:, 3], label = r"$a_0$") # a shu
 

plt.xscale('log')
plt.yscale('log')

plt.xlabel('x')
plt.legend()

plt.show()