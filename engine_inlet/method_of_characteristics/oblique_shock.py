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
            if beta*deflec < 0: raise ValueError("beta and deflec must have same signs") 
            self.solve_weak_oblique_shock()
            return 

        if beta is not None and deflec is None: 
            self.deflec = self.get_flow_deflection(self.beta, self.M1, self.gam)
            if beta < 0: self.deflec = self.deflec*-1
        
        self.solve_weak_oblique_shock()

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

        k = 1
        if deflec < 0: k = -1
        
        deflec = abs(deflec)
        
        #check if deflection is greater than maximum possible deflection 
        #func = lambda beta, M, gam: -1*self.get_flow_deflection(beta, M, gam)
        #res = scipy.optimize.minimize(func, x0=0.5*(math.pi/2 + mu), args=(M,gam))
        #betaThetMax = float(res.x) #shock wave angle for max deflection
        #deflecMax = self.get_flow_deflection(float(res.x), M, gam)
        term = (0.5*(gam+1) - math.cos(2*mu)) - math.sqrt(gam+1)*math.sqrt((0.5*(gam+1) - math.cos(2*mu))**2 + gam*(3-gam)/4)
        beta_max = math.acos(term/gam)/2
        deflecMax = self.get_flow_deflection(beta_max, M, gam)

        if deflec > abs(deflecMax):
            print(f"Warning: For Upstream Mach Number ({M}), Deflection Angle ({math.degrees(deflec)} deg) is greater than max deflection: ({math.degrees(deflecMax)} deg). Returning 90 deg wave angle.")
            return k*math.pi/2, None
        
        #calculate strong and weak shock solutions
        #def thetBetaM(beta, deflec, M, gam):
        #    return (2/math.tan(beta))*(M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta)) + 2) - math.tan(deflec)
        #beta_weak = scipy.optimize.root_scalar(thetBetaM, args=(abs(deflec),M,gam), method='bisect', bracket=[mu, betaThetMax])
        #beta_strong = scipy.optimize.root_scalar(thetBetaM, args=(abs(deflec),M,gam), method='bisect', bracket=[betaThetMax, math.pi/2])
        delta = 1
        lam = math.sqrt((M**2 - 1)**2 - 3*(1 + 0.5*(gam-1)*M**2)*(1 + 0.5*(gam+1)*M**2)*math.tan(deflec)**2)
        chi = ( (M**2 - 1)**3  - 9*(1 + 0.5*(gam-1)*M**2)*(1 + 0.5*(gam-1)*M**2 + 0.25*(gam+1)*M**4)*math.tan(deflec)**2 )/(lam**3)
        beta_weak = math.atan((M**2 - 1 + 2*lam*math.cos((4*math.pi*delta + math.acos(chi))/3))/(3*(1 + 0.5*(gam-1)*M**2)*math.tan(deflec)))
        delta = 0 
        beta_strong = math.atan((M**2 - 1 + 2*lam*math.cos((4*math.pi*delta + math.acos(chi))/3))/(3*(1 + 0.5*(gam-1)*M**2)*math.tan(deflec)))

        return k*beta_weak, k*beta_strong

    def get_flow_deflection(self, beta, M, gam): 
        """
        gives the flow deflection angle (relative to initial flow direction) given a shock wave angle and upstream mach number
        """
        absBeta = abs(beta) 
        deflec = math.atan(2*((M**2*math.sin(absBeta)**2 - 1)/(M**2*(gam + math.cos(2*absBeta)) + 2))/math.tan(absBeta))
        return deflec

    def solve_weak_oblique_shock(self): 
        """
        calculates property ratios across oblique shock wave
        """
        M1 = self.M1
        deflec = self.deflec
        gam = self.gam 

        if hasattr(self, 'beta') == False: 
            #calculate wave angles
            self.beta_w, self.beta_s = self.get_shock_wave_angle(M1, deflec, gam) #weak shock solution
            self.beta = self.beta_w 

        Mn1 = M1*math.sin(abs(self.beta))
        Mn2 = math.sqrt((1 + 0.5*(gam-1)*Mn1**2)/(gam*Mn1**2 - 0.5*(gam-1)))
        self.M2 =  Mn2/math.sin(abs(self.beta - deflec))

        self.p2_p1 = 1 + 2*gam*(Mn1**2 - 1)/(gam+1)
        self.rho2_rho1 = (gam+1)*Mn1**2/(2+(gam-1)*Mn1**2)
        self.T2_T1 = (1+2*gam*(Mn1**2 - 1)/(gam+1))*(2 + (gam-1)*Mn1**2)/((gam+1)*Mn1**2)
        
        if hasattr(self, 'R'): 
            c_p = self.R/(1-1/gam)
            self.deltaS = c_p*math.log(self.T2_T1)

        p01_p1 = (1 + 0.5*(gam-1)*M1**2)**(gam/(gam-1))
        p02_p2 = (1 + 0.5*(gam-1)*self.M2**2)**(gam/(gam-1))
        self.p02_p01 = p02_p2*self.p2_p1/p01_p1 

pass 