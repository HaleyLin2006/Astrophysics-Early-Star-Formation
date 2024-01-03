import numpy as np
from scipy.integrate import quad
import os
import sys

sys.path.insert(1, '/Users/haleylin/Desktop/ASTR Sphinx/ASTR_code/')
from parameter import param
from data_handler import data

param = param()
data = data()

class star():
    def __init__(self):
        pass
    
    def get_r_star(self):
        r_star = []
        for t in [i/31556926 for i in param.time]:
                if t<20000:
                    r_star.append(5*param.AU)
                elif t>20000 and t<20100:
                    print("here")
                    r_star.append(5*param.AU*(1-np.sqrt((t-20000)/100)) + param.r_ps)
                else:
                    r_star.append(param.r_ps)
    
        # plt.plot([i/param.sun_R for i in r_star])
        return r_star

    def luminosity_accretion(self, mass):
        def m_d_ode(u):
            return np.sqrt(1-u) / u**(4/3)
        
        r_star = self.get_r_star()
        centrifugal_radius = [1/16 * param.cs * (0.975**3) * t**3 * param.omega**2 for t in param.time]
        
        l_acc = []
        for t in range(len(param.time)):
            dmdt = 0.975 * param.cs**3 / param.G
            l0 = param.G * mass[t] * dmdt / r_star[t]
            u_star = r_star[t] / centrifugal_radius[t]

            eta_star    = 0.5
            eta_d       = 0.75

            # print(quad(m_d_ode, u_star, 1))
            m_d = 1/3 * u_star**(1/3) * quad(m_d_ode, u_star, 1)[0]
            # print(quad(m_d_ode, u_star, 1))

            luminosity = 1/(6*u_star) * (6*u_star - 2 + (2-5*u_star)*np.sqrt(1-u_star))
            luminosity += eta_star/2 * (1 - (1-eta_d)*m_d) * (1 - (1-eta_d)*np.sqrt(1-u_star))
            luminosity *= l0

            l_acc.append(luminosity)
            
        # print([i/(3.846*1e33) for i in l_acc])
        
        return l_acc
    
    def luminosity(self, mass):
        return self.luminosity_accretion(mass)

    def disk_inner_radius(self, mass):
        t_evap = 2000
        l_acc = self.luminosity(mass)

        r = []
        for t in range(len(param.time)):
            r.append(np.sqrt(l_acc[t] / (4*np.pi * param.bolzmann * t_evap**4)))
        
        # print([i/param.AU for i in r])
        # plt.plot([i/param.AU for i in r])

        return r