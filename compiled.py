import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.integrate import odeint, ode
from scipy.special import legendre
from scipy import interpolate
import numpy as np

#universal constant in cgs
AU = 1.496e+13
G = 6.67e-8
M = 1.989e+33   # 1 solar mass
cs = 0.2*1e5    # sound speed in Shu (1977) [cm/s]
time = [1e12, 2e12, 4e12, 8e12] # time in Shu (1977) [s]
color = ['blue', 'orange', 'green', 'red']

class penston():
    def __init__(self):
        self.r = []
        self.rho = []
        self.v = []

        self.density()
        self.v_r()

    def density(self):
        beta = []
        for b in range(0, 801):
            b /= 100

            beta.append(b)
            rho_i = 3/(3+10*b+7*(b**2))
            r_i = np.sqrt(b*((1+b)**(2/3)))
            self.r.append(r_i)
            self.rho.append(rho_i)

            b *= 100
    
    def v_r(self):
        for i in self.r:
            self.v.append(i**(1/7))

    def get_r(self):
        return self.r
    
    def get_density(self):
        return self.rho
    
    def get_v_r(self):
        return self.v

class shu():
    # solve the EOM with z=1/x parameter provided by Shu (1977, S77)
    #   z_g1  (-)   : a list of radius the EOM will be solve on for z>1
    #   z_s1  (-)   : a list of radius the EOM will be solve on for z<1
    #   cs ([cm/s]) : sound speed given by S77 to bring dimensionless solution back with dimension
    #   time ([s])  : time is provided to bring dimensionless solution back with dimension
    #   num_data_pts: number of data points used
    def __init__(self, z_g1, z_s1, cs=0.2*1e5, time=[1e12, 2e12, 4e12, 8e12]):
        self.z_g1 = z_g1
        self.z_s1 = z_s1
        self.cs = cs
        self.t = time

        self.store_dadz = []
        self.store_dvdz = []
        def shu_ode(z, param):
            a = param[0]
            v = param[1]

            dadz = -a/((1-z*v)**2 - z**2) * ((1/z-v) * (a - 2*z*(1/z-v)))
            dvdz = -1/((1-z*v)**2 - z**2) * ((1/z-v) * (a*(1/z-v) - 2*z))

            if dadz<0:
                dadz = 0

            #global store_dadz
            self.store_dadz.append(dadz)
            #global store_dvdz
            self.store_dvdz.append(dvdz)

            # progress recorder
            print("shu, z:", z)

            return [dadz, dvdz]

        # initial condition for Shu, z -> 0 (x -> infty)
        A = 2.0001
        param0_g1 = [A*(self.z_g1[0]**2), -(A-2)*self.z_g1[0]]

        # Solution for Shu, z<1, (x>1)
        sol_shu_g1 = np.zeros((len(self.z_g1), len(param0_g1)))
        sol_shu_g1[0, :] = param0_g1
        f_shu_g1 = ode(shu_ode)
        f_shu_g1.set_integrator("vode", method = "bdf", atol=1e-20, nsteps=1000)
        f_shu_g1.set_initial_value(param0_g1, self.z_g1[0])

        for i in range(1, len(self.z_g1)):
            sol_shu_g1[i, :] = f_shu_g1.integrate(self.z_g1[i])
            f_shu_g1.set_initial_value(sol_shu_g1[i, :], self.z_g1[i])

        self.a_shu_g1 = [i for i in sol_shu_g1[:, 0]]
        self.v_shu_g1 = [i for i in sol_shu_g1[:, 1]]

        self.store_dadz = []
        self.store_dvdz = []

        # initial condition for Shu, z>1, (x<1)
        param0_s1 = [self.a_shu_g1[-1], self.v_shu_g1[-1]]

        # Solution for Shu, z>1, (x<1)
        sol_shu_s1 = np.zeros((len(self.z_s1), len(param0_s1)))
        sol_shu_s1[0, :] = param0_s1
        f_shu_s1 = ode(shu_ode)
        f_shu_s1.set_integrator("vode", method="bdf", atol=1e-20, nsteps=1000)
        f_shu_s1.set_initial_value(param0_s1, self.z_s1[0])

        for i in range(1, len(self.z_s1)):
            sol_shu_s1[i, :] = f_shu_s1.integrate(self.z_s1[i])
            f_shu_s1.set_initial_value(sol_shu_s1[i, :], self.z_s1[i])

        self.a_shu_s1 = [i for i in sol_shu_s1[:, 0]]
        self.v_shu_s1 = [i for i in sol_shu_s1[:, 1]]

    # return dadz x<1
    def get_dadz(self):
        return self.store_dadz
    
    # return dvdz x<1
    def get_dvdz(self):
        return self.store_dvdz

    ###### Radius #####
    # return dimensionless radius x>1
    def get_z_g1(self):
        return self.z_g1
    # return radius with dimension [cm] for x>1
    def get_r_g1(self):
        r = []
        for t in self.t:
            r.append([1/i*self.cs*t for i in self.z_g1])
        
        return r
    
    # return dimensionless radius x<1
    def get_z_s1(self):
        return self.z_s1
    # return radius with dimension [cm] for x<1
    def get_r_s1(self):
        r = []
        for t in self.t:
            r.append([1/i*self.cs*t for i in self.z_s1])
        
        return r

    ###### Density #####
    # return dimensionless parametric solution: density x>1
    def get_a_g1(self):
        return self.a_shu_g1
    # return density with dimension [g/cm^3] x>1
    def get_rho_g1(self):
        rho = []
        for t in self.t:
            rho.append([(i/(4*np.pi*G*t**2)) for i in self.a_shu_g1])
        
        return rho

    # return dimensionless parametric solution: density x<1
    def get_a_s1(self):
        return self.a_shu_s1
    # return density with dimension [g/cm^3] x<1
    def get_rho_s1(self):
        rho = []
        for t in self.t:
            rho.append([(i/(4*np.pi*G*t**2)) for i in self.a_shu_s1])
        
        return rho
    
    ###### Radial Velocity #####
    # return dimensionless pramatric solution: velocity x>1
    def get_v_g1(self):
        return self.v_shu_g1
    # return velocity with dimension [cm/s] x>1
    def get_u_g1(self):
        u = []
        for t in self.t:
            u.append([-(self.cs*i) for i in self.v_shu_g1])

        return u

    # return dimensionless pramatric solution: velocity x<1
    def get_v_s1(self):
        return self.v_shu_s1
    # return velocity with dimension [cm/s] x<1
    def get_u_s1(self):
        u = []
        for t in self.t:
            u.append([-(self.cs*i) for i in self.v_shu_s1])

        return u

class ulrich():
    # EOM of Ulrich 1976
    #   x           : radius that will be calculated upon\
    #   cs  ([cm/s]): sound speed 
    #   time        : time that will be calulated with
    #   m           : mass of star (set to one solar mass)
    #   angle       : angle that will be considered
    #   shu         : set if compare with Shu, if so, then shu=True
    #
    # note that the general struture for the out put is [time][angle][value at each radius]
    # however, Ulrich solution generally do not depend on time
    def __init__(self, x, cs=0.2*1e5, m=M, time=[1e12, 2e12, 4e12, 8e12], angle=[np.pi/6, np.pi/4, np.pi/3, np.pi/2], shu=False):
        self.angle = angle
        self.m = m
        self.time = time
        self.x = x
        self.cs = cs

        self.rho = self.density(shu)
        self.velo_r = self.vel_r()
        self.velo_theta = self.vel_theta()
        self.velo_phi = self.vel_phi()
    
    def density(self, shu=False):
        a = []
        for t in range(len(self.time)):
            a.append([])
            for j in range(len(self.angle)):
                theta = self.angle[j]
                a[-1].append([])
                radius = [i*cs*time[t] for i in self.x]         # /AU ?
                for r in radius:
                    theta_0 = np.arccos(np.roots([50*AU/r, 0, 1-(50*AU/r), -np.cos(theta)])[2])
                    # dm for ulrich
                    if shu:
                        dm = 2e-6*self.m/3.154e7
                    if not shu:
                        dm = 1e-7*self.m/3.154e+7*np.exp(-np.cos(theta_0)/2)

                    d = dm/(8*np.pi*np.sqrt(G*self.m*r**3))
                    d *= (1+abs(np.cos(theta)/np.cos(theta_0)))**(-1/2)
                    #print(np.cos(theta_0))
                    d *= (np.cos(theta)/(2*np.cos(theta_0))+(50*AU/r)*np.cos(theta_0)**2)**(-1)

                    a[-1][-1].append(d)
    
        return a
    
    def vel_r(self):
        v_r = []
        for t in range(len(self.time)):
            v_r.append([])
            for j in range(len(self.angle)):
                theta = self.angle[j]
                v_r[-1].append([])
                radius = [i*cs*time[t] for i in self.x]
                for r in radius:
                    theta_0 = np.arccos(np.roots([50*AU/r, 0, 1-(50*AU/r), -np.cos(theta)])[2])
                    v_r[-1][-1].append(np.sqrt(self.m*G/r)*np.sqrt(1+np.cos(theta)/np.cos(theta_0)))

        return v_r
    
    def vel_theta(self):
        v_th = []
        for t in range(len(self.time)):
            v_th.append([])
            for j in range(len(self.angle)):
                theta = self.angle[j]
                v_th[-1].append([])
                radius = [i*cs*time[t] for i in self.x]
                for r in radius:
                    theta_0 = np.arccos(np.roots([50*AU/r, 0, 1-(50*AU/r), -np.cos(theta)])[2])
                    v_th[-1][-1].append(np.sqrt(M*G/r)*(np.cos(theta_0)-np.cos(theta))*np.sqrt((np.cos(theta_0)+np.cos(theta))/(np.cos(theta_0)*np.sin(theta))))

        return v_th

    def vel_phi(self):
        v_ph = []
        for t in range(len(self.time)):
            v_ph.append([])
            for j in range(len(self.angle)):
                theta = self.angle[j]
                v_ph[-1].append([])
                radius = [i*cs*time[t] for i in self.x]
                for r in radius:
                    theta_0 = np.arccos(np.roots([50*AU/r, 0, 1-(50*AU/r), -np.cos(theta)])[2])
                    v_ph[-1][-1].append(np.sqrt(M*G/r)*(np.sin(theta_0)/np.sin(theta_0))*np.sqrt(1-np.cos(theta)/np.cos(theta_0)))
        
        return v_ph

    def get_r(self):
        r = []
        for t in range(len(self.time)):
            r.append([i*self.cs*time[t]/AU for i in self.x])

        return r
    
    def get_a(self):
        return self.rho

    def get_v_r(self):
        return self.velo_r
    
    def get_v_theta(self):
        return self.velo_theta

    def get_v_phi(self):
        return self.velo_phi

class terebey():
    # EOM for Terebey, Shu, Cassen (1984)
    #   z_g1, z_s1              : the radius coordinate to calculate upon, in z=1/x
    #   a_shu_s1                : the Shu (1977) solution for density a0 (array)
    #   v_shu_s1                : the Shu (1977) solution for velocity v0 (array)
    #   store_dadz, store_dvdz  : dadz and dvdz from Shu (1977) calculation
    #   cs ([cm/s]) : sound speed given by S77 to bring dimensionless solution back with dimension
    #   time ([s])  : time is provided to bring dimensionless solution back with dimension
    def __init__(self, z_g1, z_s1, a_shu_s1, a_shu_g1, v_shu_s1, store_dadz, store_dvdz, cs=0.2*1e5, time=[1e12, 2e12, 4e12, 8e12]):
        # quadrupolar EOM for x>1 
        self.z_g1 = z_g1
        self.z_s1 = z_s1
        self.a_shu_s1 = a_shu_s1
        self.a_shu_g1 = a_shu_g1
        self.cs = cs
        self.time = time
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
            if (dadz < 0):
                dadz = -dadz
            
            return [dadz, dvdz, dwdz, dqdz, dpdz]

        # quadrupolar EOM for x<1
        def quadrupolar_smaller1_ode(z, param, a0, v0, da0dz, dv0dz):
            a = param[0]
            v = param[1]
            w = param[2]
            q = param[3]
            p = param[4]

            m0 = 1/z**2 * a0 * (1/z-v0)
            psi = -1/(3*z**2) - q*z**3 - p/z**2
            dpsidz = 2/(3*z) - 3*q*z**2 + 2*p/z

            A = a*z**2*(2*v0/z - dv0dz) + v*z**2*(2*a0/z - da0dz) - 6*z*a0*w
            B = -a/(a0**2)*da0dz*z**2 + (2-z**2*dv0dz)*v - z**2*dpsidz - 2/3*(m0/2)**4*z**3

            # starting from equation 69, 1. change to z=1/x parameter, 2. adapt boundary conditions from equation 26(b), 67(b), and 68\
            dadz = -1/((1-z*v0)**2 - z**2) * ((1/z-v0)*A + a0*B)
            dvdz = -1/((1-z*v0)**2 - z**2) * ((1/z-v0)*B + A/a0)
            dwdz = 1/(z*v0-1) * ((2/z+v0)*w + a/a0 + psi + (m0/2)**4/3*z**2)
            dqdz = -a/5/z**6
            dpdz = -a/5/z

            # progress recorder
            print("quadruploar x<1, z:", z)

            return [dadz, dvdz, dwdz, dqdz, dpdz]

        # monopolar EOM for x>1
        def monopolar_greater1_ode(param, z): 
            a = param[0]
            v = param[1]
            m = param[2]

            dpsidz = -m - 1/(6*z**3)

            B = a/z + 2*v -z**2*dpsidz - 2/(3*z)

            dadz = 1/(z**2-1) * (a*B)
            dvdz = 1/(z**2-1) * (B/z)
            dmdz = -1/z**4 * (a-1/2)

            return [dadz, dvdz, dmdz]

        # monopolar EOM for x<1
        def monopolar_smaller1_ode(z, param, a0, v0, da0dz, dv0dz):
            a = param[0]
            v = param[1]
            m = param[2]

            m0 = 1/z**2 * a0 * (1/z-v0)
            dpsidz = -m - 1/(6*z**3)

            # 1. change to z=1/x parameter, 2. following equation 73 and 74, following constraints of equation 78, 79
            A = a*z**2*(2*v0/z - dv0dz) + v*z**2*(2*a0/z - da0dz)
            B = a/(a0**2)*z**2*da0dz + (2-z**2*dv0dz)*v - z**2*dpsidz - 2/3*(m0/2)**4*z**3

            dadz = 1/((1-z*v0)**2 - z**2) * ((1/z-v0)*A + a0*B)
            dvdz = -1/((1-z*v0)**2 - z**2) * ((1/z-v0)*B + A/a0)
            dmdz = -1/z**4 * (a-1/2)

            # cheating?    
            #if (dadz < 0):
            #    dvdz = -dadz

            print("monopolar x<1, z:", z)
            
            return [dadz, dvdz, dmdz]

        # boundary condition for quandrupolar, x>1 model
        # adopting equation 65 and 70, K given in III-c of the paper
        K = 1.5*10**(-3)
        # a, v, w, q, p
        param0_g1 = [0, 0, 0, K, 0]

        # solution for quadrupolar, x>1 model
        sol_g1 = odeint(quadrupolar_greater1_ode, param0_g1, z_g1)
        self.a_g = sol_g1[:, 0]
        self.v_g = sol_g1[:, 1]
        self.w_g = sol_g1[:, 2] #continuous through x=1

        # boaundary condition for quadrupolar, x<1 model
        # adopting equation 67 and 68 to go over the x=1 discontinuity. a_shu and v_shu from Shu(1976)
        # here, we manually take limit by sampling points close to 1, such that a(1+e) = a(1) = a(1-e), where e is an arbitrary small number
        # a, v, w, q, p
        param0_s1 = [self.a_g[-1]+2/3*self.v_g[-1], 2/3*self.v_g[-1], self.w_g[-1], K*(1+1/35), -2*K/245]

        # solution for quadrupolar, x<1 model
        sol_s1 = np.zeros((len(z_s1), len(param0_s1)))
        sol_s1[0, :] = param0_s1
        f_s1 = ode(quadrupolar_smaller1_ode)
        f_s1.set_initial_value(param0_s1, z_s1[0])
        f_s1.set_integrator("vode", method = "bdf", atol=1e-20, nsteps=1000)
        f_s1.set_f_params(a_shu_s1[0], v_shu_s1[0], store_dadz[0], store_dvdz[0])

        for i in range(1, len(z_s1)):
            sol_s1[i, :] = f_s1.integrate(z_s1[i])
            f_s1.set_initial_value(sol_s1[i, :], z_s1[i])
            f_s1.set_f_params(a_shu_s1[i], v_shu_s1[i], store_dadz[i], store_dvdz[i])

        self.a_s = sol_s1[:, 0]
        self.v_s = sol_s1[:, 1]
        self.w_s = sol_s1[:, 2]

        param0_mono_g1 = [1/2, 0, 0]
        # solution for monopolar, x>1 model
        sol_mono_g1 = odeint(monopolar_greater1_ode, param0_mono_g1, z_g1)
        self.a_mono_g = sol_mono_g1[:, 0]
        self.v_mono_g = sol_mono_g1[:, 1]
        self.m_mono_g = sol_mono_g1[:, 2] #continuous through x=1

        # initial condition for monopolar x<1 model
        param0_mono_s1 = [self.a_mono_g[-1], self.v_mono_g[-1], self.m_mono_g[-1]]

        # solution for monopolar, x<1 model
        sol_mono_s1 = np.zeros((len(z_s1), len(param0_mono_s1)))
        sol_mono_s1[0, :] = param0_mono_s1
        f_mono_s1 = ode(monopolar_smaller1_ode)
        f_mono_s1.set_integrator("vode", method = "bdf", atol=1e-20, nsteps=1000)
        f_mono_s1.set_initial_value(param0_mono_s1, z_s1[0])
        f_mono_s1.set_f_params(a_shu_s1[0], v_shu_s1[0], store_dadz[0], store_dvdz[0])

        for i in range(1, len(z_s1)):
            sol_mono_s1[i, :] = f_mono_s1.integrate(z_s1[i])
            f_mono_s1.set_initial_value(sol_mono_s1[i, :], z_s1[i])
            f_mono_s1.set_f_params(a_shu_s1[i], v_shu_s1[i], store_dadz[i], store_dvdz[i])

        self.a_mono_s1 = [i for i in sol_mono_s1[:, 0]]
        self.v_mono_s1 = [i for i in sol_mono_s1[:, 1]]

    def get_r_g1(self):
        r = []
        for t in self.time:
            r.append([1/i*self.cs*t for i in self.z_g1])
        
        return r

    def get_r_s1(self):
        r = []
        for t in self.time:
            r.append([1/i*self.cs*t for i in self.z_s1])
        
        return r

    ##### Density #####
    def get_quad_a_g1(self):
        return self.a_g
    
    def get_quad_a_s1(self):
        return self.a_s
    
    def get_mono_a_g1(self):
        return self.a_mono_g
    
    def get_mono_a_s1(self):
        return self.a_mono_s1
    
    ##### Velocity #####
    def get_quad_v_g1(self):
        return self.v_g
    
    def get_quad_v_s1(self):
        return self.v_s
    
    def get_quad_w_g1(self):
        return self.w_g
    
    def get_quad_w_s1(self):
        return self.w_s
    
    def get_mono_v_g1(self):
        return self.v_mono_g
    
    def get_mono_v_s1(self):
        return self.v_mono_s1
    
    def get_density_g1(self, theta, tau=1e-3):
        # for omega=1e-14, tau^2=1e-3
        a_2 = []
        leg = legendre(2)

        for i in range(len(self.a_g)):
            a_2_i = self.a_mono_g[i] + self.a_g[i]*leg(np.cos(theta))
            a_2.append(a_2_i)
        
        rho = []
        for i in range(len(a_2)):
            rho_i = self.a_shu_g1[i] + tau*a_2[i]
            rho.append(rho_i)

        density = []
        for t in self.time:
            density.append([(i/(4*np.pi*G*t**2)) for i in rho])

        return density
    
    def get_density_s1(self, theta, tau=1e-3):
        # for omega=1e-14, tau=1e-3
        a_2 = []
        leg = legendre(2)

        for i in range(len(self.a_s)):
            a_2_i = self.a_mono_s1[i] + self.a_s[i]*leg(np.cos(theta))
            a_2.append(a_2_i)
        
        rho = []
        for i in range(len(a_2)):
            rho_i = self.a_shu_s1[i] + tau*a_2[i]
            rho.append(rho_i)

        density = []
        for t in self.time:
            density.append([(i/(4*np.pi*G*t**2)) for i in rho])

        return density

class lodato():
    # surface density of Lodato (2007)
    #   x       : coordinate to be examine on
    #   cs      : sound speed
    #   time    : time
    #   C       : normalization constant
    def __init__(self, x, cs=0.2*1e5, time=[1e12, 2e12, 4e12, 8e12], C=1):
        self.r = x
        self.time = time
        self.cs = cs

        self.sigma = []
        t_l = [1, 2, 4, 8]
        count = 0
        for t in time:
            r1 = cs*t/AU
            t_lodato = t_l[count]
            count += 1
            self.sigma.append([])
            for i in range(len(self.r)):
                sigma_i = C*t_lodato**(-3/2)/(3*np.pi*self.r[i]) * np.exp(-(self.r[i]/r1)/t_lodato)
                self.sigma[-1].append(sigma_i)
    
    def get_sigma(self):
        return self.sigma

sample_pts = 1000
sample_range = 4        # in interval [1e(-sample_range), 1e(sample_range)]

# z=1/x
#z_g1 = np.logspace(-sample_range, 0, sample_pts)
z_g1 = [2/(i+1) for i in np.logspace(sample_range, 0, sample_pts)]  # also x_s1
z_s1 = np.logspace(0, sample_range, sample_pts)                     # also x_g1

########### Shu 1977 ##########
cs = 0.26*1e5    # sound speed in Visser (2009) [cm/s]
t_acc = 2.5*1e5 * 31556926 # t_acc in Visser (2009) [s]
time = [0.05*t_acc, 0.1*t_acc, 0.2*t_acc, 0.3*t_acc, 0.5*t_acc, 0.7*t_acc, t_acc] # time in Visser (2009) [s]

s = shu(z_g1=z_g1, z_s1=z_s1, cs=cs, time=time)
a_shu_s1 = s.get_a_s1()
v_shu_s1 = s.get_v_s1()
a_shu_g1 = s.get_a_g1()
v_shu_g1 = s.get_v_g1()
dadz = s.get_dadz()
dvdz = s.get_dvdz()

# dimensionless x, dimensionless a and v
# plt.plot([1/i for i in s.get_z_g1()], s.get_a_g1())
# plt.plot([1/i for i in s.get_z_s1()], s.get_a_s1())

# plt.plot([1/i for i in s.get_z_g1()], [-i for i in s.get_v_g1()])
# plt.plot([1/i for i in s.get_z_s1()], [-i for i in s.get_v_s1()])


color = ['blue', 'orange', 'green', 'red', 'purple', 'cyan', 'grey']
# r [cm], a [g/cm^-3] and v [cm/s]
for i in range(len(s.get_r_g1())):
    i = 0
    # density
    # plt.plot([j/AU for j in s.get_r_g1()[i]] + [j/AU for j in s.get_r_s1()[i]], s.get_rho_g1()[i]+ s.get_rho_s1()[i], color=color[i], label=f't={time[i]/1e12}e12')
    #plt.plot([j/AU for j in s.get_r_s1()[i]], s.get_rho_s1()[i], color=color[i])

    # velocity
    # plt.plot(s.get_r_g1()[i], s.get_u_g1()[i])
    # plt.plot(s.get_r_s1()[i], s.get_u_s1()[i], label=f't={time[i]/1e12}e12')

##### Shu (1977) 2D Graph #####
# Replace this with your actual 1-D data
fig, axes = plt.subplots(2, 4)
axes_flat = axes.flatten()

common_vmin = min(np.min(image) for image in s.get_rho_g1())
common_vmax = max(np.max(image) for image in s.get_rho_s1())

fig, axes = plt.subplots(2, 4)
axes_flat = axes.flatten()

for i, ax in enumerate(axes_flat):
    if i < len(time):
        radii = np.linspace(0, 2000, 2000)  # Distances from the origin
        molecule_counts = s.get_rho_s1()[i][::-1] + s.get_rho_g1()[i][::-1]

        # Create a grid
        grid_size = sample_pts  # Adjust as needed
        x = np.linspace(0, 2000, 2000)
        y = np.linspace(0, 2000, 2000)
        xx, yy = np.meshgrid(x, y)
        r = np.sqrt(xx**2 + yy**2)

        # Interpolate data onto the grid
        interpolated_data = np.interp(r.ravel(), radii, molecule_counts, left=0, right=0)
        interpolated_data = interpolated_data.reshape(xx.shape)

        # Logarithmic colormap
        norm = LogNorm(vmin=common_vmin, vmax=common_vmax)

        # Plot the radially symmetrical image
        im = ax.imshow(interpolated_data, cmap='viridis', origin='lower', extent=[0, 2000, 0, 2000], norm=norm)

        ax.set_title(f't = {time[i]/t_acc}')

        # Add contour line at the value 10^(-18)
        contour_value = 1e-18
        ax.contour(xx, yy, interpolated_data, levels=[contour_value], colors='red', linewidths=2)

    else:
        ax.axis('off')
        cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
        cbar.set_label('Intensity')



########### Ulrich 1976 ##########
# u_s1 = ulrich(x=z_g1, cs=cs, time=time, shu=True)
# u_g1 = ulrich(x=z_s1, cs=cs, time=time, shu=True)

# # density
# plt.plot(u_s1.get_r()[i], u_s1.get_a()[i][1], color='purple')
# plt.plot(u_g1.get_r()[i], u_g1.get_a()[i][1], color='purple')

# velocity
# plt.plot(u_s1.get_r()[i], u_s1.get_v_r()[i][1], color='purple')
# plt.plot(u_g1.get_r()[i], u_g1.get_v_r()[i][1], color='purple')

########### Terebey 1984 ##########
# t = terebey(z_g1=z_g1, z_s1=z_s1, a_shu_s1=a_shu_s1, a_shu_g1=a_shu_g1, v_shu_s1=v_shu_s1, store_dadz=dadz, store_dvdz=dvdz, cs=cs, time=time)
# tsc_quad_a_s1 = t.get_quad_a_s1()
# tsc_quad_a_g1 = t.get_quad_a_g1()
# tsc_mono_a_s1 = t.get_mono_a_s1()
# tsc_mono_a_g1 = t.get_mono_a_g1()

# tsc_quad_v_s1 = t.get_quad_v_s1()
# tsc_quad_v_g1 = t.get_quad_v_g1()
# tsc_quad_w_s1 = t.get_quad_w_s1()
# tsc_quad_w_g1 = t.get_quad_w_g1()

# tsc_mono_v_s1 = t.get_mono_v_s1()

# tsc_density_g1 = t.get_density_g1(theta=np.pi/4)
# tsc_density_s1 = t.get_density_s1(theta=np.pi/4)

# # graph with dimension
# count = 0
# for i in range(len(t.get_r_g1())):
#     plt.plot([j/AU for j in t.get_r_g1()[i]], tsc_density_g1[i], color=color[i], label=f't={time[i]/1e12}e12')
#     plt.plot([j/AU for j in t.get_r_s1()[i]], tsc_density_s1[i], color=color[i])

# graph in the paper
# graph_type = 0
# if graph_type == 1:
#     # Shu
#     plt.plot([1/i for i in z_g1], [i for i in a_shu_g1], color='orange', label='a shu', linewidth=2)
#     plt.plot([1/i for i in z_s1], [i for i in a_shu_s1], color='orange', linewidth=2)
#     # quadrupolar
#     plt.plot([1/i for i in z_g1], [i for i in tsc_quad_a_g1], label=r"$-\alpha_Q$", color='red', linewidth=2)
#     plt.plot([1/i for i in z_s1], [abs(i) for i in tsc_quad_a_s1], color='red', linewidth=2)                 
#     # monopolar
#     plt.plot([1/i for i in z_g1], [i for i in tsc_mono_a_g1], label=r"$a_m$", color='blue', linewidth=2)
#     plt.plot([1/i for i in z_s1], [i for i in tsc_mono_a_s1], color='blue', linewidth=2)

# if graph_type == 2:
#     plt.plot([1/i for i in z_g1], [i for i in tsc_quad_v_g1], label=r"$v_Q$", color= 'blue', linewidth=2)
#     plt.plot([1/i for i in z_s1], [i for i in tsc_quad_v_s1], color='blue', linewidth=2)       
    
#     plt.plot([1/i for i in z_g1], [i for i in tsc_quad_w_g1], label=r"$-w_Q$", color='green', linewidth=2)
#     plt.plot([1/i for i in z_s1], [-i for i in tsc_quad_w_s1], color='green', linewidth=2)

# if graph_type == 3:
#     plt.plot([1/i for i in z_s1], [-i for i in v_shu_s1], color='orange', label='-v shu', linewidth=2)
#     plt.plot([1/i for i in z_s1], [i for i in tsc_mono_v_s1], label = r"$v_m$", linewidth=2)

########### Lodato 2007 ##########
# l_s1 = lodato(x=z_g1, time=time)
# l_g1 = lodato(x=z_s1, time=time)

# l_s1_sigma = l_s1.get_sigma()
# l_g1_sigma = l_g1.get_sigma()
# for i in range(len(l_s1_sigma)):
#     plt.plot([x*cs*time[i]/AU for x in z_g1], l_s1_sigma[i], label=f't={time[i]}', color='red')
# for j in range(len(l_g1_sigma)):
#      plt.plot([x*cs*time[j]/AU for x in z_s1], l_g1_sigma[j], label=f't={time[j]}', color='red')


# plt.xscale('log')
# plt.yscale('log')
plt.show()