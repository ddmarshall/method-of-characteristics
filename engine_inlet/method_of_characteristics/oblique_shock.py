import math 
import scipy.optimize 
"""
This module contains functions useful for calculating flow properties of oblique shock waves
"""
class Oblique_Shock:

    def __init__(self, M1, gam, R=None, deflec=None, beta=None):
        """
        M1: upstream mach number
        gam: specific heat ratio (calorically perfect)
        R: ideal gas constant 
        deflec: flow deflection angle from upstream velocity 
        beta: shock wave angle relative to upstream velocity
        """
        self.M1 = M1 
        self.gam = gam 
        if R is not None: self.R = R
        if deflec is not None: self.deflec = deflec 
        if beta is not None: self.beta = beta

        if beta is None and deflec is None: 
            return 

        if beta is not None and deflec is not None: 
            self.solve_weak_oblique_shock(M1, deflec, gam, R=R, beta=abs(beta))
            return 

        if beta is not None and deflec is None: 
            self.get_flow_deflection(abs(beta), M1, gam)
            if beta < 0: self.deflec = self.deflec*-1
        
        self.solve_weak_oblique_shock(M1, self.deflec, gam, R=R)

    def get_shock_wave_angle(self, M, deflec, gam): 
        """
        calculates the shock wave angle required for a given upstream mach number, flow deflection from upstream angle
        Inputs
            M: upstream mach number 
            deflec: (rad) flow deflection from upstream direction
            gam: ratio of specific heats
        Return: 
            beta: (rad) shock wave angle 
        """
        #check to see if given thet and M are allowable
        mu = math.asin(1/M) #mach angle 
        
        func = lambda beta, M, gam: -1*self.get_flow_deflection(beta, M, gam)
        res = scipy.optimize.minimize(func, x0=0.5*(math.pi/2 + mu), args=(M,gam))
        betaThetMax = float(res.x) #shock wave angle for max deflection
        deflecMax = self.get_flow_deflection(float(res.x), M, gam)
        if deflec > deflecMax:
            raise ValueError("Upstream Mach Number and Deflection Angle Will Result in Detached Shock")
        
        def thetBetaM(beta, deflec, M, gam):
            return (2/math.tan(beta))*(M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta)) + 2) - math.tan(deflec)

        beta_weak = scipy.optimize.root_scalar(thetBetaM, args=(deflec,M,gam), method='bisect', bracket=[mu, betaThetMax])
        beta_strong = scipy.optimize.root_scalar(thetBetaM, args=(deflec,M,gam), method='bisect', bracket=[betaThetMax, math.pi/2])

        return beta_weak.root, beta_strong.root

    def get_flow_deflection(self, beta, M, gam): 
        """
        gives the flow deflection angle (relative to initial flow direction) given a shock wave angle and upstream mach number
        beta: shock wave angle
        M: upstream mach number 
        gam: specific heat ratio
        """
        deflec = math.atan(2/math.tan(beta))*(M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta)) + 2)
        self.deflec = abs(deflec)
        return abs(deflec) 

    def solve_weak_oblique_shock(self, M1, deflec, gam, R=None, beta=None): 
        """
        calculates property ratios across oblique shock wave
        M1: upstream mach number 
        deflec: flow deflection from upstream velocity 
        gam: specific heat ratio 
        R: ideal gas constant 
        beta: shock wave angle relative to upstream velocity
        """
        if beta is None: 
            self.beta_w, self.beta_s = self.get_shock_wave_angle(M1, abs(deflec), gam) #weak shock solution
            if deflec < 0: 
                self.beta_w = self.beta_w*-1
                self.beta_s = self.beta_s*-1
        else: 
            self.beta_w = beta

        Mn1 = M1*math.sin(abs(self.beta_w))
        Mn2 = math.sqrt((1 + 0.5*(gam-1)*Mn1**2)/(gam*Mn1**2 - 0.5*(gam-1)))
        self.M2 =  Mn2/math.sin(abs(self.beta_w - deflec))

        self.p2_p1 = 1 + 2*gam*(Mn1**2 - 1)/(gam+1)
        self.rho2_rho1 = (gam+1)*Mn1**2/(2+(gam-1)*Mn1**2)
        self.T2_T1 = (1+2*gam*(Mn1**2 - 1)/(gam+1))*(2 + (gam-1)*Mn1**2)/((gam+1)*Mn1**2)
        
        if R is not None: 
            c_p = R/(1-1/gam)
            self.deltaS = c_p*math.log(self.T2_T1)

        p01_p1 = (1 + 0.5*(gam-1)*M1**2)**(gam/(gam-1))
        p02_p2 = (1 + 0.5*(gam-1)*self.M2**2)**(gam/(gam-1))
        self.p02_p01 = p02_p2*self.p2_p1/p01_p1