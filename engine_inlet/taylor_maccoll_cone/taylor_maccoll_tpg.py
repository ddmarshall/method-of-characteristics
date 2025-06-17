import scipy.integrate
import scipy.interpolate
import scipy.optimize 
import math
import numpy as np
from engine_inlet.method_of_characteristics.gas_state_TPG import ThermallyPerfectGas as GS
import engine_inlet.method_of_characteristics.oblique_shock_tpg as obs_tpg

from scipy.optimize import least_squares



class TaylorMaccoll_Cone:
    """
    This class calculates the Taylor-Maccoll flowfield for an infinite, straight cone while 
    under Thermally Perfect Gas assumption 
    """
    def __init__(self, cone_ang, M_inf, gasProps):
        """
        cone_ang = cone half angle in radians 
        M_inf = free-stream Mach number 
        gasProps = object contiaining gas properties 
        """
        self.cone_ang = cone_ang #cone half angle
        self.M_inf = M_inf #freestream mach number 
        self.gam = gasProps.gam #specific heat ratio 
        self.R = gasProps.R #ideal gas constant
        self.T0 = gasProps.T0 #freestream stagnation temperature 
        self.Tv = gasProps.Tv #vibrational characteristic temperature
        self.gas_state = gasProps


        self.h0 = GS._h(self.R, self.Tv, self.T0)


        #solve tmc flow and obtain ratios
        self.solve_TMC_flow()

        self.p02_p01 = self.shock_values.p02_p01
     #   self.obtain_flow_properties()
     #   self.convert_velocity_to_cartesian()
       
    def solve_TMC_flow(self): 
            
        def TMC_flow(thet, V, R, Tv, h0, T0):
            #Specifies system of ODEs for numerical solver to integrate (verified, don't touch)
            #Stop Condition triggered when V_thet --> 0
            V_r, V_thet = V

            #Find h given the velocity
            #NOTE: max is there to catch cases when the solver goes past the root and 
            #needs points that would result in crossing h0 before it gets to the stop condition
            h = max(h0 - (V_r**2 + V_thet**2) / 2, 1e-4) 

            def f(T):
                return h - GS._h(R, Tv, T)
        
            initial_guess = float(T0 / 2)
            bounds = (0, T0)
            root = least_squares(f, initial_guess, bounds=bounds)
            T = root.x[0]

            #Overflow Error Avoidance 
            if T <=300 :
                gam = 1.4
            else:
                gam = GS._gamma(R,Tv,T)
            a_sq = gam * R * T

            eq1 = V_thet
            eq2 = (2 * a_sq * V_r + a_sq * V_thet / np.tan(thet) - V_thet**2 * V_r) / (V_thet**2 - a_sq)

            return [eq1, eq2]

        def TMC_cone_guess(shock_ang, cone_ang, M_inf, gasProps, R, ret):
            """
            TODO: Fix issue where IVP fails to capture cone surface with certain inputs (make process more robust)
            Solves TMC cone flow using prescribed shock angle. Returns cone angle error or solution
            shock_ang: prescribed shock wave angle (rads)
            cone_ang: true cone half-angle (rads)
            M_inf: freestream Mach number 
            gam: specific heat ratio 
            """

            #get conditions directly after shock
            #Needed values across the shock: Mn1 
            

            Mn1 = M_inf*math.sin(shock_ang) #normal freestream mach component
    
            if Mn1 <= 1: #if shock angle and freestream mach number result in a subsonic or sonic normal mach component abort and return nothing 
                return None

            #Using guessed shock angle, solve for the deflection and the other values across the shock 
            shock_values = obs_tpg.Oblique_Shock(M_inf, gasProps, beta=shock_ang)

            #pull out flow deflection 
            flow_deflec = shock_values.deflec

            V_thet_init = -shock_values.velocity2*math.sin(shock_ang-flow_deflec)
            V_r_init = shock_values.velocity2*math.cos(shock_ang-flow_deflec)

            #initial conditions necessary for the ODE solver
            y0 = [V_r_init, V_thet_init]

            final_angle = 1e-4 #angle for solver to integrate to (should never reach angle = 0)

            #Solver's real stop condition (when V_thet --> 0)
            def stop_condition(t, y, R, Tv, h0, T0): 
                return y[1] - 1e-4
            stop_condition.terminal = True

            sol = scipy.integrate.solve_ivp(TMC_flow, (shock_ang, final_angle), y0, args=[self.R, self.Tv, self.h0, self.T0], dense_output=True, events=stop_condition) #dense output turned on
            lastpoint = sol.y[1].max()

            
            if sol.y[1].max() < 0:
                #Solution scheme reached the end of the integration without seeing the theoretical cone angle
                 
                return None
                raise ValueError("IVP Solve Failed To Capture Theoretical Cone Surface")

            func = lambda thet: sol.sol(thet)[1] #returns V_theta for a given theta
            cone_ang_exp = scipy.optimize.bisect(func, sol.t.min(), sol.t.max()) #find cone angle based on shock angle 
            
            if ret=="error":
                return abs(cone_ang_exp - cone_ang)
            elif ret=="solution":
                return sol, shock_values
            else: 
                raise ValueError("Invalid ret specified")

        def TMC_cone_flow(cone_ang, M_inf, R):
            """
            Computes the flow solution for taylor maccoll flow around an infinite cone
            cone_ang = cone half angle (rads)
            M_inf = free stream Mach number 
            gam = specific heat ratio
            plotting = turn solver output plotting on (set to True) 
            """
            alpha_inf = math.asin(1/M_inf)
            #shock_ang_est = 1*(cone_ang + 0.5*alpha_inf) #initial guess shock angle 
            shock_ang_est = 1.1*alpha_inf

            #TODO: root_scalar? 
            fsolve_output = scipy.optimize.fsolve(TMC_cone_guess, x0=shock_ang_est, args=(cone_ang, M_inf, self.gas_state, R, "error"), full_output=True)
            shock_ang = float(fsolve_output[0])    
            
            #run function with correct shock angle:
            solution, shock_values = TMC_cone_guess(shock_ang, cone_ang, M_inf, self.gas_state,  R, "solution")
            return [solution, shock_ang, shock_values]

        self.numerical_solution, self.shock_ang, self.shock_values = TMC_cone_flow(self.cone_ang, self.M_inf, self.R)
    
    def obtain_flow_properties(self):
        """
        TODO add functions to find flow properties for variable angle theta 
        Obtains isentropic pressure, density, and temperature relations given a velocity solution
        Equations from Anderson Intro to Aerodynamics Ch 8
        """

        shock_values = obs_tpg.Oblique_Shock_tpg(self.M_inf, self.gas_state, R=self.R, beta=self.shock_ang)
        #Mach Number on Cone Surface and directly behind shock: 
        def Mach(Vrp, Vthetp): 
            h = self.h0 - (Vrp**2 + Vthetp**2) / 2
            def f(T):
                return h - GS._h(self.R, self.Tv, T)
        
            initial_guess = float(self.T0 / 2)
            bounds = (0, self.T0)
            root = least_squares(f, initial_guess, bounds=bounds)
            T = root.x

            a_sq = GS._gamma(self.R, self.Tv, T) * self.R * T
            return np.sqrt((Vrp**2 + Vthetp**2) / a_sq)
        
        def Mach_thet(thet):
            Vrp, Vthetp = self.numerical_solution.sol(thet)
            return Mach(Vrp, Vthetp, self.gam)

        self.f_mach = Mach_thet #function to continuously get mach number

        [V_R_c, V_thet_c] = self.numerical_solution.sol(self.cone_ang)
        M_c = Mach(V_R_c, V_thet_c, self.gam) 

        [V_R_2, V_thet_2] = self.numerical_solution.sol(self.shock_ang)
        M_2 = Mach(V_R_2, V_thet_2, self.gam)
        

        ##TODO: Ask marshall about this function 
        shock_turn_ang = shock_values.deflec #shock turn angle 

        #temperature 
        def T0_T_thet(thet): #function to continuously get stagnation to static temperature ratio 
            def f(T):
                return h - GS._h(self.R, self.Tv, T)
            Vrp, Vthetp = self.numerical_solution.sol(thet)
            h = self.h0 - np.sqrt(Vrp**2 + Vthetp**2) / 2
            initial_guess = float(self.T0 / 2)
            bounds = (0, self.T0)
            root = least_squares(f, initial_guess, bounds=bounds)
            T = root.x
            return self.T0 / T
        
        self.f_T0_T = T0_T_thet

        
        T02_T2 = self.T0 / (shock_values.T2_T1 * self.gas_state.T1)
        T2_T1 = shock_values.T2_T1
        T0c_Tc = T0_T_thet([0])
        Tc_T1 = (1/T0c_Tc)*T02_T2*T2_T1


        #Pressure: 
        p0_p = lambda M, gam: (1 + M**2*(gam-1)/2)**(gam/(gam-1)) #isentropic total:static pressure ratio
        def p0_p_thet(thet):
            M = Mach_thet(thet)
            return p0_p(M, self.gam)
        self.f_p0_p = p0_p_thet


        p01_p1       = p0_p(self.M_inf, self.gam)
        p02_p2       = p0_p(M_2, self.gam)
        #M1_n         = self.M_inf*math.sin(self.shock_ang)#get freestream normal mach component 
        p2_p1        = shock_values.p2_p1#static pressure ratio across shock
        p02_p01      = p02_p2*p2_p1/p01_p1 #stagnation pressure ratio across shock 
        p0c_pc       = p0_p(M_c, self.gam)
        pc_p1        = (1/p0c_pc)*p02_p2*p2_p1 #cone surface static pressure vs freestream static pressure 
        p0c_p01      = pc_p1*p0c_pc*(1/p01_p1) #cone surface total pressure vs freestream total pressure 


        #density
        rho0_rho = lambda M, gam: (1 + (gam-1)*M**2/2)**(1/(gam-1))
        def rho0_rho_thet(thet):
            M = Mach(thet)
            return rho0_rho(M, self.gam)
        self.f_rho0_rho = rho0_rho_thet

        rho2_rho1_normal = lambda M, gam: (gam+1)*M**2/(2 + (gam-1)*M**2)

        rho02_rho2 = rho0_rho(M_2, self.gam)
        rho2_rho1 = shock_values.rho2_rho1
        rho0c_rhoc = rho0_rho(M_c, self.gam)
        rhoc_rho1 = (1/rho0c_rhoc)*rho02_rho2*rho2_rho1

        #store values 
        self.M_c = M_c
        self.shock_turn_ang = shock_turn_ang
        self.p2_p1 = p2_p1
        self.p02_p01 = p02_p01
        self.pc_p1 = pc_p1
        self.p0c_p01 = p0c_p01
        self.rho2_rho1 = rho2_rho1
        self.rhoc_rho1 = rhoc_rho1
        self.T2_T1 = T2_T1
        self.Tc_T1 = Tc_T1


    def f_veloc_uv(self, thet): 
        [Vrp, Vthetp] = self.numerical_solution.sol(thet)
        V = np.sqrt(Vrp**2 + Vthetp**2)
        phi = math.atan(Vthetp/Vrp)
        u = V*math.cos(phi + thet)
        v = V*math.sin(phi + thet)
        return u,v

    