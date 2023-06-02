import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

#universal constant in cgs
AU = 1.496e+13
G = 6.67e-8
M = 1.989e+33 # 1 solar mass

class penston69():
    # notice that penston is dimensionless
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

class shu76():
    # a= speed of sound; t= time interval to calculate; p0= initial condition 
    def __init__(self, a=0.2*1e5, t=[1e12, 2e12, 4e12, 8e12], p0=[-5.44, 71.5], num_data_point=20):
        def ode(p, x):
            v = p[0]
            a = p[1]

            dvdx = (a*(x-v)-(2/x))*(x-v)/((x-v)**2-1)
            dadx = (a*(a-(2/x)*(x-v))*(x-v))/((x-v)**2-1)

            return [dvdx, dadx]
        
        x = [i/num_data_point for i in range(1, num_data_point)]
        self.t = t
        self.r = [[i*a*t[0] for i in x], [i*a*t[1] for i in x], [i*a*t[2] for i in x], [i*a*t[3] for i in x]]

        sol = odeint(ode, p0, x)

        v = sol[:, 0]
        self.u = [-(a*i) for i in v]
        alpha = sol[:, 1]
        self.rho = [[(i/(4*np.pi*6.67e-8*t[0]**2)) for i in alpha], [(i/(4*np.pi*6.67e-8*t[1]**2)) for i in alpha], [(i/(4*np.pi*6.67e-8*t[2]**2)) for i in alpha], [(i/(4*np.pi*6.67e-8*t[3]**2)) for i in alpha]]

    def zero_to_nan(self, values):
        """Replace every 0 with 'nan' and return a copy."""
        return [float('nan') if x==0 else x for x in values]
    
    def get_t(self):
        return self.t
    
    def get_r(self):
        return self.r
    
    def get_density(self):
        return self.rho
    
    def get_v_r(self):
        return self.u



t = [1e12, 2e12, 4e12, 8e12]
x = [i/20 for i in range(1, 20)]
angle = [np.pi/6, np.pi/4, np.pi/3, np.pi/2]

class ulrich76():
    #set shu=True when compare to shu
    def __init__(self, a=0.2*1e5, t=[1e12, 2e12, 4e12, 8e12], M=1.989e+33 ,angle=[np.pi/6, np.pi/4, np.pi/3, np.pi/2], num_data_point=20, shu=False):
        self.angle = angle
        self.M = M
        self.t = t
        x = [i/num_data_point for i in range(1, num_data_point)]
        self.r = [[i*a*t[0] for i in x], [i*a*t[1] for i in x], [i*a*t[2] for i in x], [i*a*t[3] for i in x]]
        #structure of v_r
        #[[[value at different r]=angle]=time]; [time][angle]=[value at different r]
        self.v_r = [[[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))]]
        self.v_theta = [[[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))]]
        self.v_phi = [[[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))]]
        self.rho = [[[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))], [[] for i in range(len(angle))]]

        self.density(shu)
        self.vel_r()
        self.vel_theta()
        self.vel_phi()

    def density(self, shu=False):
        for time in range(len(t)):
            for j in range(len(self.angle)):
                theta = self.angle[j]
                for i in self.r[time]:
                    theta_0 = np.arccos(np.roots([50*AU/i, 0, 1-(50*AU/i), -np.cos(theta)])[2])
            # dm for ulrich
                    if shu:
                        dm = 2e-6*M/3.154e7
                    if not shu:
                        dm = 1e-7*M/3.154e+7*np.exp(-np.cos(theta_0)/2)

            #print("dm=", dm, "theta=", theta, "theta_0=", theta_0)
                    d = dm/(8*np.pi*np.sqrt(G*M*i**3))
                    d *= (1+np.cos(theta)/np.cos(theta_0))**(-1/2)
                    d *= (np.cos(theta)/(2*np.cos(theta_0))+(50*AU/i)*np.cos(theta_0)**2)**(-1)

                    self.rho[time][j].append(d)
    
    def vel_r(self):
        for time in range(len(t)):
        # len(t) = len(r)
            for j in range(len(self.angle)):
                theta = self.angle[j]
                for i in self.r[time]:
                    theta_0 = np.arccos(np.roots([50*AU/i, 0, 1-(50*AU/i), -np.cos(theta)])[2])
                    self.v_r[time][j].append(np.sqrt(M*G/i)*np.sqrt(1+np.cos(theta)/np.cos(theta_0)))
    
    def vel_theta(self):
        for time in range(len(t)):
            for j in range(len(self.angle)):
                theta = self.angle[j]
                for i in self.r[time]:
                    theta_0 = np.arccos(np.roots([50*AU/i, 0, 1-(50*AU/i), -np.cos(theta)])[2])
                    self.v_theta[time][j].append(np.sqrt(M*G/i)*(np.cos(theta_0)-np.cos(theta))*np.sqrt((np.cos(theta_0)+np.cos(theta))/(np.cos(theta_0)*np.sin(theta))))

    def vel_phi(self):
        for time in range(len(t)):
            for j in range(len(self.angle)):
                theta = self.angle[j]
                for i in self.r[time]:
                    theta_0 = np.arccos(np.roots([50*AU/i, 0, 1-(50*AU/i), -np.cos(theta)])[2])
                    self.v_phi[time][j].append(np.sqrt(M*G/i)*(np.sin(theta_0)/np.sin(theta_0))*np.sqrt(1-np.cos(theta)/np.cos(theta_0)))

    
    def get_r(self):
        return self.r

    def get_t(self):
        return self.t
    
    def get_density(self):
        return self.rho

    def get_v_r(self):
        return self.v_r
    
    def get_v_theta(self):
        return self.v_theta

    def get_v_phi(self):
        return self.v_phi

class tereby84():
    def __init__(self, K=1.5*10**(-3), num_data_point=20):        
        # for x>1
        def quadrupolar_greater1_ode(param, z):
            a = param[0]
            v = param[1]
            w = param[2]
            q = param[3]
            p = param[4]

            dadz = 1/(z**2-1) * ((-12*z**2*w) + 2*z**2*(-a/z + 2*v -3*q*z**4 + 2*p/z))
            dvdz = 1/(z**2-1) * (1/z*(-a/z + 2*v -3*q*z**4 + 2*p/z) - 6*z*w)
            dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
            dqdz = -5*a/z**6
            dpdz = a/5/z

            if (dadz <0):
                dadz = 0

            return [dadz, dvdz, dwdz, dqdz, dpdz]
                
        self.z = [i/1000 for i in range(1, 1000, 990//num_data_point)]
        param0 = [0, 0, 0, K, 0]

        sol = odeint(quadrupolar_greater1_ode, param0, z)

        self.a = sol[:, 0]
        self.v = sol[:, 1]

    def get_z(self):
        return self.z

    def get_quadr_g1_density(self):
        return self.a
    
    def get_quadru_g1_vel(self):
        return self.v



shu = shu76()
shu_r = shu.get_r()
shu_t = shu.get_t()
shu_rho = shu.get_density()
shu_v_r = shu.get_v_r()

ulrich = ulrich76(shu=True)
ulrich_r = ulrich.get_r()
ulrich_t = ulrich.get_t()
ulrich_rho = ulrich.get_density()
ulrich_v_r = ulrich.get_v_r()
ulrich_v_theta = ulrich.get_v_theta()
ulrich_v_phi = ulrich.get_v_phi()

marker = ["o", "^", "s", "*"]
angle_name = [30, 45, 60, 90]
color = ['red', 'blue', 'green', 'purple']
for time in range(len(t)):
    plt.plot([i/AU for i in ulrich_r[time]], ulrich_rho[time][1], marker[time], label=f'{t[time]/1e12}e12 {angle_name[1]}', color=color[time])
    plt.plot([i/AU for i in shu_r[time]], shu_rho[time])

    #plt.plot([i/AU for i in shu_r[time]], shu_v_r, label=f'shu, {t[time]/1e12}e12')

#plt.plot([i/AU for i in ulrich_r[time]], ulrich_v_r[time][1], marker="o", label=r"$v_r$")

#plt.plot([i/AU for i in ulrich_r[time]], ulrich_v_theta[time][1], marker="^", label=r"$v_\theta$")

#plt.plot([i/AU for i in ulrich_r[time]], ulrich_v_phi[time][1], marker="s", label=r"$v_\phi$")

plt.xscale('log')
plt.yscale('log')
plt.ylabel('log|density| (cgs)', fontsize=20)
plt.xlabel('log(AU)', fontsize=20)
plt.legend(loc='lower left')
plt.show()
