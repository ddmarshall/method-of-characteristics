import scipy.integrate
import scipy.interpolate
import scipy.optimize 
import math
import matplotlib.pyplot as plt
import numpy as np

class TaylorMaccoll_Cone:
    """
    This class calculates the Taylor-Maccoll flowfield for an infinite, straight cone 
    """
    def __init__(self, cone_ang, M_inf, gam, R, T0):
        """
        cone_ang = cone half angle in radians 
        M_inf = free-stream Mach number 
        gam = specific heat ratio (calorically perfect)
        """
        self.cone_ang = cone_ang #cone half angle
        self.M_inf = M_inf #freestream mach number 
        self.gam = gam #specific heat ratio 
        self.R = R #ideal gas constant
        self.T0 = T0 #freestream stagnation temperature 

        #solve tmc flow and obtain ratios
        self.solve_TMC_flow()
        self.obtain_flow_properties()
        self.convert_velocity_to_rectangular()
       
    def solve_TMC_flow(self): 

        def TMC_flow(thet, V, gam):
            #Specifies system of ODEs for numerical solver to integrate (verified, don't touch)
            V_r, V_thet = V
            eq1 = V_thet
            k = (gam-1)/2
            eq2 = -(k*(2*V_r**3 - 2*V_r + 2*V_r*V_thet**2 + (V_thet*V_r**2 - V_thet + V_thet**3)/math.tan(thet)) + V_r*V_thet**2)/(k*(V_r**2+V_thet**2-1) + V_thet**2)

            return [eq1, eq2]

        def TMC_cone_guess(shock_ang, cone_ang, M_inf, gam, ret):
            """
            TODO: Fix issue where IVP fails to capture cone surface with certain inputs (make process more robust)
            Solves TMC cone flow using prescribed shock angle. Returns cone angle error or solution
            shock_ang: prescribed shock wave angle (rads)
            cone_ang: true cone half-angle (rads)
            M_inf: freestream Mach number 
            gam: specific heat ratio 
            """
            #get conditions directly after shock
            Mn1 = M_inf*math.sin(shock_ang) #normal freestream mach component
            Mn2 = math.sqrt((1 + 0.5*(gam-1)*Mn1**2)/(gam*Mn1**2 - 0.5*(gam-1))) #normal shock relation 
            flow_deflec = math.atan((2*(1/math.tan(shock_ang))*((M_inf**2)*(math.sin(shock_ang)**2) - 1))/((M_inf**2)*(gam + math.cos(2*shock_ang)) + 2))

            M2 = Mn2/math.sin(shock_ang - flow_deflec)

            V_nondim = math.sqrt((((gam-1)/2)*(M2**2))/(1 + ((gam-1)/2)*(M2**2)))

            V_thet_init = -V_nondim*math.sin(shock_ang-flow_deflec)
            V_r_init = V_nondim*math.cos(shock_ang-flow_deflec)

            y0 = [V_r_init, V_thet_init]
            gam = 1.4 #specific heat ratio 
            final_angle = 2*cone_ang/3 #angle for solver to integrate to (must be beyond cone angle)

            sol = scipy.integrate.solve_ivp(TMC_flow, (shock_ang, final_angle), y0, args=[gam], dense_output=True) #dense output turned on

            if sol.y[1].min() > 0 or sol.y[1].max() < 0:
                #Check to see if solution will lead to a theoretical cone angle 
                #return None
                raise ValueError("IVP Solve Failed To Capture Theoretical Cone Surface")

            func = lambda thet: sol.sol(thet)[1] #returns V_theta for a given theta
            cone_ang_exp = scipy.optimize.root_scalar(func, method='bisect', bracket=(sol.t.min(), sol.t.max())) #find cone angle based on shock angle 

            if ret=="error":
                return abs(cone_ang_exp.root - cone_ang)
            elif ret=="solution":
                return sol
            else: 
                raise ValueError("Invalid ret specified")

        def TMC_cone_flow(cone_ang, M_inf, gam):
            """
            Computes the flow solution for taylor maccoll flow around an infinite cone
            cone_ang = cone half angle (rads)
            M_inf = free stream Mach number 
            gam = specific heat ratio
            plotting = turn solver output plotting on (set to True) 
            """
            alpha_inf = math.asin(1/M_inf)
            shock_ang_est = cone_ang + 0.5*alpha_inf #initial guess shock angle 
            fsolve_output = scipy.optimize.fsolve(TMC_cone_guess, x0=shock_ang_est, args=(cone_ang, M_inf, gam, "error"), full_output=True)
            shock_ang = float(fsolve_output[0])          
            #TODO switch to root_scalar 

            print(f"found shock angle: {round(math.degrees(shock_ang),3)} (deg)")
            
            #run function with correct shock angle:
            solution = TMC_cone_guess(shock_ang, cone_ang, M_inf, gam, "solution")
            return [solution, shock_ang]

        self.numerical_solution, self.shock_ang = TMC_cone_flow(self.cone_ang, self.M_inf, self.gam)
    
    def obtain_flow_properties(self):
        """
        TODO add functions to find flow properties for variable angle theta 
        Obtains isentropic pressure, density, and temperature relations given a velocity solution
        Equations from Anderson Intro to Aerodynamics Ch 8
        """
        #Mach Number on Cone Surface and directly behind shock: 
        Mach = lambda V_R, V_thet, gam: math.sqrt((2/(gam-1))/(-1 + 1/(math.sqrt(V_R**2 + V_thet**2)**2)))#anderson eq 13.81 rearranged
        
        def Mach_thet(thet):
            Vrp, Vthetp = self.numerical_solution.sol(thet)
            return Mach(Vrp, Vthetp, self.gam)

        self.f_mach = Mach_thet #function to continuously get mach number

        [V_R_c, V_thet_c] = self.numerical_solution.sol(self.cone_ang)
        M_c = Mach(V_R_c, V_thet_c, self.gam) 

        [V_R_2, V_thet_2] = self.numerical_solution.sol(self.shock_ang)
        M_2 = Mach(V_R_2, V_thet_2, self.gam)
        
        shock_turn_ang = math.atan(2/math.tan(self.shock_ang)*(self.M_inf**2*math.sin(self.shock_ang)**2 - 1)/(self.M_inf**2*(self.gam + math.cos(2*self.shock_ang)) + 2))#shock turn angle 

        #Pressure: 
        p0_p = lambda M, gam: (1 + M**2*(gam-1)/2)**(gam/(gam-1)) #isentropic total:static pressure ratio
        def p0_p_thet(thet):
            M = Mach_thet(thet)
            return p0_p(M, self.gam)
        self.f_p0_p = p0_p_thet

        p2_p1_normal = lambda M, gam: 1 + (2*gam/(gam+1))*(M**2 - 1) #normal shock static pressure ratio 
        p01_p1 = p0_p(self.M_inf, self.gam)
        p02_p2 = p0_p(M_2, self.gam)
        M1_n = self.M_inf*math.sin(self.shock_ang)#get freestream normal mach component 
        p2_p1 = p2_p1_normal(M1_n, self.gam) #static pressure ratio across shock
        p02_p01 = p02_p2*p2_p1/p01_p1 #stagnation pressure ratio across shock 
        p0c_pc = p0_p(M_c, self.gam)
        pc_p1 = (1/p0c_pc)*p02_p2*p2_p1 #cone surface static pressure vs freestream static pressure 
        p0c_p01 = pc_p1*p0c_pc*(1/p01_p1) #cone surface total pressure vs freestream total pressure 

        #temperature 
        T0_T = lambda M, gam: 1 + (gam-1)*M**2/2
        def T0_T_thet(thet): #function to continuously get stagnation to static temperature ratio 
            M = Mach_thet(thet)
            return T0_T(M, self.gam)
        self.f_T0_T = T0_T_thet

        T2_T1_normal = lambda M,gam: (1+2*gam*(M**2-1)/(gam+1))*(2 + (gam-1)*M**2)/((gam+1)*M**2)
        
        T02_T2 = T0_T(M_2, self.gam)
        T2_T1 = T2_T1_normal(M1_n, self.gam)
        T0c_Tc = T0_T(M_c, self.gam)
        Tc_T1 = (1/T0c_Tc)*T02_T2*T2_T1

        #density
        rho0_rho = lambda M, gam: (1 + (gam-1)*M**2/2)**(1/(gam-1))
        def rho0_rho_thet(thet):
            M = Mach(thet)
            return rho0_rho(M, self.gam)
        self.f_rho0_rho = rho0_rho_thet

        rho2_rho1_normal = lambda M, gam: (gam+1)*M**2/(2 + (gam-1)*M**2)

        rho02_rho2 = rho0_rho(M_2, self.gam)
        rho2_rho1 = rho2_rho1_normal(M1_n, self.gam)
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

    def convert_velocity_to_rectangular(self):
        """
        TODO clean up this function (stagnation temperature doesn't change from the shock)
        returns a function which finds the rectangular velocity components u and v at any point in the flow field
        """ 
        gam = self.gam
        R = self.R
        shock_ang = self.shock_ang
        #determine temperature just upstream of cone shock assuming isentropic acceleration
        T1 = self.T0/(1 + 0.5*(gam-1)*self.M_inf**2) #static temperature of freestream
        T2 = T1*self.T2_T1 #temperature immediately behind shock
        M2 = self.f_mach(shock_ang) #mach number immediately behind shock 
        T0 = T2*(1 + 0.5*(gam-1)*M2**2) #get stagnation temperature

        #define function to convert velocities
        def get_veloc_uv(thet): 
            [Vrp, Vthetp] = self.numerical_solution.sol(thet)
            M = self.f_mach(thet)
            T = T0/(1 + 0.5*(gam-1)*M**2)
            a = math.sqrt(gam*R*T)
            V = M*a
            Vmax = math.sqrt(2*(a**2/(gam-1) + V**2/2))
            Vr, Vthet = Vrp*Vmax, Vthetp*Vmax
            u = Vr*math.cos(thet) - Vthet*math.sin(thet)
            v = Vr*math.sin(thet) + Vthet*math.cos(thet) 
            return u,v

        self.f_veloc_uv = get_veloc_uv
        
if __name__ == "__main__":
    """
    TODO: Get rid of this section when module is functional
    """
    cone_ang = math.radians(30)
    M_inf = 3
    gam, R, T0 = 1.4, 287.05, 288.15

    def plot_TMC_flow(cone_obj):
        """
        Plots the continuous solution between the shock and cone surface
        cone_obj: taylor maccoll solved flow object
        """
        dense_output = cone_obj.numerical_solution.sol
        #create a plot of V_theta, V_R vs theta
        n = 100 #number of theta slices 
    
        theta_arr = np.linspace(cone_obj.cone_ang, cone_obj.shock_ang, n)

        V_R_arr = np.array([dense_output(thet)[0] for thet in theta_arr])
        V_thet_arr = np.array([dense_output(thet)[1] for thet in theta_arr])

        plt.figure(), plt.grid()
        plt.plot([math.degrees(theta) for theta in theta_arr], V_R_arr, label="V_R")
        plt.plot([math.degrees(theta) for theta in theta_arr], V_thet_arr, label="V_thet")
        plt.show()

    cone = TaylorMaccoll_Cone(cone_ang, M_inf, gam, R, T0) 
    