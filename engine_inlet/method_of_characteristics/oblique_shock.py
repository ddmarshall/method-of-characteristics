import math 
import scipy.optimize 
"""
This module contains functions useful for calculating flow properties of shock waves
"""
class Oblique_Shock:

    def __init__(self, M1, gam, R=None, thet=None, beta=None):
        
        self.M1 = M1 
        self.gam = gam 
        if R is not None: self.R = R
        if thet is not None: self.thet = thet 
        if beta is not None: self.beta = beta

        if beta is None and thet is None: 
            return 

        if beta is not None and thet is not None: 
            self.solve_weak_oblique_shock(M1, thet, gam, R=R, beta=beta)
            return 

        if beta is not None and thet is None: 
            self.get_flow_deflection(beta, M1, gam)
        
        self.solve_weak_oblique_shock(M1, self.thet, gam, R=R)



    def get_shock_wave_angle(self, M, thet, gam): 
        """
        calculates the shock wave angle required for a given upstream mach number, flow deflection from upstream angle
        Inputs
            M: upstream mach number 
            thet: (rad) flow deflection (positive number) from upstream direction
            gam: ratio of specific heats
        Return: 
            beta: (rad) shock wave angle 
        """
        #check to see if given thet and M are allowable
        mu = math.asin(1/M) #mach angle 
        
        """
        res = scipy.optimize.minimize(lambda beta, M, gam: -1*self.get_flow_deflection(beta, M, gam), x0=0.5*(math.pi/2 + mu), args=(M,gam))
        thetMax = self.get_flow_deflection(float(res.x), M, gam)
        if thet > thetMax:
            raise ValueError("Upstream Mach Number and Deflection Angle Will Result in Detached Shock")
        """
        def thetBetaM(beta, thet, M, gam):
            return (2/math.tan(beta))*(M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta)) + 2) - math.tan(thet)

        beta_weak = scipy.optimize.root_scalar(thetBetaM, args=(thet,M,gam), method='bisect', bracket=[mu, 0.5*(mu + math.pi/2)])
        beta_strong = scipy.optimize.root_scalar(thetBetaM, args=(thet,M,gam), method='bisect', bracket=[0.5*(mu + math.pi/2), math.pi/2])

        return beta_weak.root, beta_strong.root



    def get_flow_deflection(self, beta, M, gam): 
        """
        gives the flow deflection angle (relative to initial flow direction) given a shock wave angle and upstream mach number
        """
        thet = math.atan(2*((M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta) + 2)))/math.tan(beta))
        self.thet = thet



    def solve_weak_oblique_shock(self, M1, thet, gam, R=None, beta=None): 

        if beta is None: 
            self.beta_w,self.beta_s = self.get_shock_wave_angle(M1, thet, gam) #weak shock solution

        else: 
            self.beta_w = beta

        Mn1 = M1*math.sin(self.beta_w)
        Mn2 = math.sqrt((1 + 0.5*(gam-1)*Mn1**2)/(gam*Mn1**2 - 0.5*(gam-1)))
        self.M2 =  Mn2/math.sin(self.beta_w - thet)

        self.p2_p1 = 1 + 2*gam*(Mn1**2 - 1)/(gam+1)
        self.rho2_rho1 = (gam+1)*Mn1**2/(2+(gam-1)*Mn1**2)
        self.T2_T1 = (1+2*gam*(Mn1**2 - 1)/(gam+1))*(2 + (gam-1)*Mn1**2)/((gam+1)*Mn1**2)
        
        if R is not None: 
            c_p = R/(1-1/gam)
            self.deltaS = c_p*math.log(self.T2_T1)

        p01_p1 = (1 + 0.5*(gam-1)*M1**2)**(gam/(gam-1))
        p02_p2 = (1 + 0.5*(gam-1)*self.M2**2)**(gam/(gam-1))
        self.p02_p01 = p02_p2*self.p2_p1/p01_p1