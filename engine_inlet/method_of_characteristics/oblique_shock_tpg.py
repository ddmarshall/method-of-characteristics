# THINGS TO ADD TO CALL FUNCTION
# gas states, Tv, delta


import math 
import scipy.optimize 
from scipy.optimize import fsolve
from . import oblique_shock as obs_cpg 
import numpy as np 
from . import gas_state_TPG as gas_state
GS = gas_state.ThermallyPerfectGas
"""
This module contains functions useful for calculating flow properties of oblique shock waves 
under the thermally perfect gas assumption 
"""
class Oblique_Shock:
   

    def __init__(self, M1, gasProp, TOL=1e-7, deflec=None, beta=None):
        """
        M1: upstream mach number
        gam: specific heat ratio (calorically perfect)
        R: ideal gas constant 
        deflec: flow deflection angle from upstream velocity 
        beta: shock wave angle relative to upstream velocity
        gas_state: the current state of the gas (Temp, Dens, Press)
        """
        #unpack gasprops information and pull state values 
        self.M1 = M1 
        self.Tv = gasProp.Tv 
        self.delta = gasProp.delta
        self.R = gasProp.R
        self.TOL = TOL
        self.p01 = gasProp.p0_p0f * gasProp.p0f
        self.T0 = gasProp.T0
        self.rho01 = self.p01 / (gasProp.R * self.T0)


        self.T1 = GS._T_given_T0_M(gasProp.R, self.Tv, self.T0, self.M1)
        self.rho1 = GS._isen_rho_given_T(self.Tv, self.T1, self.T0, self.rho01)
        self.p1 = self.rho1 * self.R * self.T1
        

        
        self.gasProp1 = gasProp #need for input into CPG obs 
        
        if gasProp.R is not None: self.R = gasProp.R
        if deflec is not None: self.deflec = deflec 
        if beta is not None: self.beta = beta

        if beta is None and deflec is None: 
            return 

        if beta is not None and deflec is not None:
            if beta*deflec < 0: raise ValueError("beta and deflec must have same signs") 
            self.solve_weak_oblique_shock()
            return 

#        if beta is not None and deflec is None: 
#            self.deflec = self.get_flow_deflection(self.beta, self.M1, self.gam)
#            if beta < 0: self.deflec = self.deflec*-1
        
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

        self.e1 = self.R * ( (3 + 2 * self.delta) / 2 * self.T1 + self.delta * (self.Tv / (np.exp(self.Tv/self.T1) - 1)) )
        
        self.Cv = self.R * ( (3 + 2 * self.delta) / 2 + self.delta * (self.Tv / (2 * self.T1 * math.sinh(self.Tv / (2 * self.T1)))**2 ))
        self.Cp = self.Cv + self.R
        self.gam = self.Cp / self.Cv
        self.a1 = math.sqrt(self.gam * self.R * self.T1)
        self.velocity1 = self.M1 * self.a1
        self.alpha_v1 = 1 / self.T1
        self.gasProp1.gam = self.gam

        if hasattr(self, 'beta') == True: 
            #calculate wave angles
            self.solve_normal_shock(self.beta)
            return

        
        if hasattr(self, 'deflec') == True:
            #calculate 4x4 Newton-Raphson method
            self.solve_weak_shock_deflec(self.deflec)
            return


        
        #if hasattr(self, 'R'): 
        #    c_p = self.R/(1-1/gam)
        #    self.deltaS = c_p*math.log(self.T2_T1)


       # self.p02_p01 = p02_p2*self.p2_p1/p01_p1 

    def solve_normal_shock(self, beta):

        #instantiate for Newton-Raphson Method
        count = 1
        #stop function if counter exceeds number
        break_point = 200
        TOL = self.TOL
        delta_x = np.array([1, 1, 1])


        #Solve CPG shock equations in order to get initial conditions 
        obs_init = obs_cpg.Oblique_Shock(self.M1, self.gasProp1, beta=beta)
        x = np.array([obs_init.p2_p1, obs_init.rho2_rho1, obs_init.T2_T1])


        while(max(abs(delta_x)) > TOL and count <= break_point):

            #Find all the second state values from gases 
            #TODO: CHANGE THESE TO A SEPARATE CLASS
            self.alpha_v2 = 1 / (x[2] * self.T1)
            self.k_T2 = 1 / (x[0] * self.p1)
            self.c_v2 = self.R * ( (3 + 2 * self.delta) / 2 + self.delta * ( 0.5 * self.Tv / (x[2] * self.T1 * np.sinh(0.5 * self.Tv / (x[2] * self.T1)) ) )**2)


            #create the jacobian Matrix 
            j11 = 0.5 * (self.p1 / (self.rho1 * self.e1)) * (1/x[1] - 1)
            j12 = (self.p1 / (self.rho1 * self.e1)) * (1 / x[1]**2) * ( (x[0] - 1) / 2 - (self.T1 / self.p1) * (x[2] * self.alpha_v2) / self.k_T2)
            j13 = (self.T1 / self.e1) * self.c_v2
            j21 = 1.0
            j22 = - (self.rho1 * self.velocity1 ** 2 * math.sin(beta) ** 2 / self.p1) * (1 / x[1] ** 2)
            j23 = 0.0
            j31 = 1.0
            j32 = - (1 / self.p1) * (1 / (x[1] * self.k_T2) )
            j33 = - (self.T1 / self.p1) * (self.alpha_v2 / self.k_T2)



            Jacobian = np.array([[j11, j12, j13], 
                                 [j21, j22, j23], 
                                 [j31, j32, j33]])

            #Create the f(x) vector 
            f1 = self.energy2_eq(x[1] * self.rho1, x[2] * self.T1) / self.e1 - 1 + (self.p1 / (2 * self.rho1 * self.e1) * (x[0] + 1) * (1/x[1] -1))
            f2 = x[0] - 1 + (self.rho1 * self.velocity1 **2 * math.sin(beta)**2 / self.p1) * (1 / x[1] - 1)
            f3 = x[0] - (x[1] * self.rho1 * self.R * x[2] * self.T1) / self.p1

            f = np.array([f1, f2, f3]) 

            #Solve the matrix equation A * delta_x  =  x for delta_x 
            delta_x = np.linalg.solve(Jacobian, f)

            #new x 
            x = x - delta_x 
            #increase counter
            count += 1

        #Pull out the values from x vector
        self.p2_p1 = float(x[0])
        self.rho2_rho1 = float(x[1])
        self.T2_T1 = float(x[2])
        #With Beta known, solve for theta
        #Define f4 and its respective derivative
        x = np.array([[self.p2_p1], [self.rho2_rho1], [self.T2_T1]])
        def f4(deflec):
            return  np.array((x[1] - 1) * math.sin(2 * self.beta - deflec) - (x[1] + 1) * math.sin(deflec))
           
        def df4(deflec):
            return -(x[1] - 1) * math.cos(2 * self.beta - deflec) - (x[1] + 1) * math.cos(deflec)
        

        initial_guess = float(obs_init.deflec)
        deflec = fsolve(f4, initial_guess, fprime=df4)
        self.deflec = deflec.item()



        h0 = GS._h(self.R, self.Tv, self.T1) + self.velocity1**2 / 2
        self.velocity2 = np.sqrt(2 * (h0 - GS._h(self.R, self.Tv, self.T2_T1 * self.T1))) 
        self.a2 = np.sqrt(GS._gamma(self.R, self.Tv, self.T2_T1 * self.T1) * self.R * self.T2_T1 * self.T1)
        self.M2 = self.velocity2 / self.a2
        self.Mn2 = self.M2 * math.sin(self.beta - self.deflec)
        #pull out the final values from x and add them to the object
        rho02 = GS._isen_rho_given_T(self.Tv, self.T0, self.T2_T1 * self.T1, self.rho2_rho1 * self.rho1)
        self.rho02_rho01 = rho02 / self.rho01
        self.p02_p01 = self.rho02_rho01
        self.T2_T0 = self.T2_T1 * self.T1 / self.T0

    
            

    def energy2_eq(self, rho2, T2): #put it into the gas properties class
        
        return self.R * ( (3 + 2 * self.delta) / 2 * (T2) + self.delta * self.Tv / (math.exp(self.Tv / T2) - 1))

    def solve_weak_shock_deflec(self, deflec):

        #instantiate for Newton-Raphson Method
        count = 1
        #stop function if counter exceeds number
        break_point = 200
        TOL = self.TOL
        delta_x = np.array([1, 1, 1, 1])


        #Solve CPG shock equations in order to get initial conditions 
        obs_init = obs_cpg.Oblique_Shock(self.M1, self.gasProp1, deflec=deflec)
        x = np.array([obs_init.p2_p1, obs_init.rho2_rho1, obs_init.T2_T1, obs_init.beta])


        while(max(abs(delta_x)) > TOL and count <= break_point):

            #Find all the second state values from gases 
            #TODO: CHANGE THESE TO A SEPARATE CLASS
            self.alpha_v2 = 1 / (x[2] * self.T1)
            self.k_T2 = 1 / (x[0] * self.p1)
            self.c_v2 = self.R * ( (3 + 2 * self.delta) / 2 + self.delta * ( 0.5 * self.Tv / (x[2] * self.T1 * np.sinh(0.5 * self.Tv / (x[2] * self.T1)) ) )**2)

            #create the jacobian Matrix 
            j11 = 0.5 * (self.p1 / (self.rho1 * self.e1)) * (1/x[1] - 1)
            j12 = (self.p1 / (self.rho1 * self.e1)) * (1 / x[1]**2) * ( (x[0] - 1) / 2 - (self.T1 / self.p1) * (x[2] * self.alpha_v2) / self.k_T2)
            j13 = (self.T1 / self.e1) * self.c_v2
            j14 = 0
            j21 = 1
            j22 = - (self.rho1 * self.velocity1 ** 2 * math.sin(x[3]) ** 2 / self.p1) * (1 / x[1] ** 2)
            j23 = 0
            j24 = (self.rho1 * self.velocity1 **2 * math.sin(2 * x[3]) / self.p1) * (1/x[1] - 1)
            j31 = 1 
            j32 = - (1 / (self.p1 * x[1] * self.k_T2) )
            j33 = - (self.T1 / self.p1) * (self.alpha_v2 / self.k_T2)
            j34 = 0
            j41 = 0
            j42 = math.sin(2 * x[3] - deflec) - math.sin(deflec)
            j43 = 0
            j44 = 2 * (x[1] - 1) * math.cos(2 * x[3] - deflec)


            Jacobian = np.array([[j11, j12, j13, j14], [j21, j22, j23, j24], [j31, j32, j33, j34], [j41, j42, j43, j44]])

            #Create the f(x) vector 
            f1 = self.energy2_eq(x[1] * self.rho1, x[2] * self.T1) / self.e1 - 1 + (self.p1 / (2 * self.rho1 * self.e1) * (x[0] + 1) * (1/x[1] -1))
            f2 = x[0] - 1 + (self.rho1 * (self.velocity1 **2 * math.sin(x[3])**2) / self.p1) * (1 / x[1] - 1)
            f3 = x[0] - (x[1] * self.rho1 * self.R * x[2] * self.T1) / self.p1
            f4 = (x[1] - 1) * math.sin(2 * x[3] - deflec) - (x[1] + 1) * math.sin(deflec)

            f = np.array([f1, f2, f3, f4]) 

            #Solve the matrix equation A * delta_x  =  x for delta_x 
            delta_x = np.linalg.solve(Jacobian, f)

            #new x 
            x = x - delta_x
            #increase counter
            count += 1
        #With Beta known, solve for theta
        #Define f4 and its respective derivative

        #pull out the final values from x and add them to the object
        self.p2_p1 = x[0].item()
        self.rho2_rho1 = x[1].item()
        self.T2_T1 = x[2].item()
        self.beta = x[3].item()

        h0 = GS._h(self.R, self.Tv, self.T1) + self.velocity1**2 / 2
        self.velocity2 = np.sqrt(2 * (h0 - GS._h(self.R, self.Tv, self.T2_T1 * self.T1))) 
        self.a2 = np.sqrt(GS._gamma(self.R, self.Tv, self.T2_T1 * self.T1) * self.R * self.T2_T1 * self.T1)
        self.M2 = self.velocity2 / self.a2
        self.Mn2 = self.M2 * math.sin(self.beta - self.deflec)

        rho02 = GS._isen_rho_given_T(self.Tv, self.T0, self.T2_T1 * self.T1, self.rho2_rho1 * self.rho1)
        self.rho02_rho01 = rho02 / self.rho01
        self.p02_p01 = self.rho02_rho01
        self.T2_T0 = self.T2_T1 * self.T1 / self.T0











pass 