import scipy.integrate
import scipy.interpolate
import scipy.optimize 
import math
import matplotlib.pyplot as plt
import numpy as np

class TaylorMaccoll_Cone:
    """
    docstring
    """
    def __init__(self, cone_ang, M_inf, gam):
        self.cone_ang = cone_ang 
        self.M_inf = M_inf
        self.gam = gam 
       
    def solve_TMC_flow(self, t_step): 

        self.t_step = t_step

        def TMC_flow(thet, V, gam):
            #Specifies system of ODEs for numerical solver to integrate
            V_r, V_thet = V
            eq1 = V_thet
            k = (gam-1)/2
            eq2 = -(k*(2*V_r**3 - 2*V_r + 2*V_r*V_thet**2 + (V_thet*V_r**2 - V_thet + V_thet**3)/math.tan(thet)) + V_r*V_thet**2)/(k*(V_r**2+V_thet**2-1) + V_thet**2)

            return [eq1, eq2]

        def TMC_cone_guess(shock_ang, cone_ang, M_inf, gam, ret, t_eval=None, plotSol=None,):
            """
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

            sol = scipy.integrate.solve_ivp(TMC_flow, (shock_ang, final_angle), y0, args=[gam], t_eval=t_eval, dense_output=True)

            if sol.y[1].min() > 0 or sol.y[1].max() < 0:
                #Check to see if solution will lead to a theoretical cone angle 
                #return None
                raise ValueError("IVP Solve Failed To Capture Theoretical Cone Surface")

            if plotSol: 
                #if plotting is enabled, plot solver solution 
                plt.figure()
                plt.plot([math.degrees(x) for x in sol.t], sol.y[0], label="V'_r"), plt.plot([math.degrees(x) for x in sol.t], sol.y[1], label="V'_\u03B8")
                plt.legend(), plt.grid()
                plt.xlabel("\u03B8 (deg)"), plt.ylabel("V'")
                plt.show() 

            #interpolate to find expected cone-half angle (theta where y_2=0)
            interpfunc = scipy.interpolate.interp1d(sol.t, sol.y[1])
            cone_ang_exp = scipy.optimize.root_scalar(interpfunc, method='bisect', bracket=(sol.t.min(), sol.t.max()))
            #TODO use dense output here, no interpolation 


            if ret=="error":
                return abs(cone_ang_exp.root - cone_ang)
            elif ret=="solution":
                return sol
            else: 
                raise ValueError("Invalid ret specified")

        def TMC_cone_flow(cone_ang, M_inf, gam, plotting=None):
            """
            Computes the flow solution for taylor maccoll flow around an infinite cone
            cone_ang = cone half angle (rads)
            M_inf = free stream Mach number 
            gam = specific heat ratio
            plotting = turn solver output plotting on (set to True) 
            """
            alpha_inf = math.asin(1/M_inf)
            shock_ang_est = cone_ang + 0.5*alpha_inf #initial guess shock angle 
            shock_ang = float(scipy.optimize.fsolve(TMC_cone_guess, x0=shock_ang_est, args=(cone_ang, M_inf, gam, "error")))
            print(f"found shock angle: {round(math.degrees(shock_ang),2)} (deg)")
            
            #run function with correct shock angle:
            theta_points = np.flip(np.arange(2*cone_ang/3, shock_ang, self.t_step))
            solution = TMC_cone_guess(shock_ang, cone_ang, M_inf, gam, "solution", t_eval=theta_points, plotSol=plotting)
            return [solution, shock_ang]

        [solution, shock_ang] = TMC_cone_flow(self.cone_ang, self.M_inf, self.gam)
        self.shock_ang = float(shock_ang)

        #Write portion of solution between cone surface and shock 
        self.numerical_Solution = solution 
        self.theta = [theta for theta in solution.t if theta >= self.cone_ang]
        self.Vprim_R = [Vp for Vp in solution.y[0] if list(solution.y[0]).index(Vp) < len(self.theta)]
        self.Vprim_thet = [Vp for Vp in solution.y[1] if list(solution.y[1]).index(Vp) < len(self.theta)]
        #self.theta = solution.t
        #self.Vprim_R = solution.y[0]
        #self.Vprim_thet = solution.y[1]

    def obtain_flow_properties(self):
        """
        Obtains isentropic pressure, density, and temperature relations given a velocity solution
        Equations from Anderson Intro to Aerodynamics Ch 8
        """
        if hasattr(self, 'numerical_Solution') == False:
            #Make sure object has necessary attributes 
            raise ValueError("Object Missing Flow Solution") 

        #Obtain V_prime magnitude
        self.Vprim = [math.sqrt(x1**2 + x2**2) for x1, x2 in zip(self.Vprim_R, self.Vprim_thet)]
        #Obtain Mach Magnitude 
        self.M = [math.sqrt((2/(self.gam-1))/(-1 + 1/(V**2))) for V in self.Vprim]

        #Obtain pressure ratio 
        self.p_p0 = [1/(1 + 0.5*(self.gam-1)*m**2)**(self.gam/(self.gam-1)) for m in self.M]
        
        #Obtain density ratio 
        self.rho_rho0 = [1/(1 + 0.5*(self.gam-1)*m**2)**(1/(self.gam-1)) for m in self.M]
        
        #Obtain temperature ratio 
        self.T_T0 = [1/(1 + 0.5*(self.gam-1)*m**2) for m in self.M]

    def obtain_ZH_flow_properties(self):

        if hasattr(self, 'numerical_Solution') == False:
            #Make sure object has necessary attributes 
            raise ValueError("Object Missing Flow Solution")

        #Obtain Mach Magnitude 
        self.M = [math.sqrt((2/(self.gam-1))/(-1 + 1/(V**2))) for V in self.Vprim]

        #?convert to Z&H mach 
        Mstr = [math.sqrt((gam+1)*M**2/(2+(gam-1)*M**2)) for M in self.M]
        
        #?Obtaining Pressure Ratios (P SURE THIS IS WRONG)
        p0_p = [(1 + 0.5*(self.gam-1)*m**2)**(self.gam/(self.gam-1)) for m in Mstr]

        #?Obtain density ratio (P SURE THIS IS WRONG)
        rho0_rho = [(1 + 0.5*(self.gam-1)*m**2)**(1/(self.gam-1)) for m in Mstr]
        pass
        
if __name__ == "__main__":
    #The following code will only execute of this file is run directly:
    #Testing out functions in class: 
    cone_ang = math.radians(30)
    M_inf = 3
    gam = 1.4
    cone_Flow = TaylorMaccoll_Cone(cone_ang, M_inf, gam) #constructing object

    theta_step = math.radians(0.1)
    cone_Flow.solve_TMC_flow(theta_step) #solving IVP 

    cone_Flow.obtain_flow_properties() #obtaining isentropic properties in flowfield
    cone_Flow.obtain_ZH_flow_properties()

    #Creating Mach vs Theta Plot
    plt.figure(), plt.plot([math.degrees(x) for x in cone_Flow.theta], cone_Flow.M)
    plt.vlines(math.degrees(cone_Flow.cone_ang), min(cone_Flow.M), max(cone_Flow.M), 
        linestyles='dashed', label='Cone Angle')
    plt.vlines(math.degrees(cone_Flow.shock_ang), min(cone_Flow.M), max(cone_Flow.M),
        colors='r', linestyles='dashed', label='Shock Angle')
    plt.grid(), plt.xlabel('\u03B8 (deg)'), plt.ylabel('Mach'), plt.legend()

    #Creating Plot for Isentropic Property Ratios
    plt.figure(), plt.plot([math.degrees(x) for x in cone_Flow.theta], cone_Flow.T_T0, label="T/T0")
    plt.plot([math.degrees(x) for x in cone_Flow.theta], cone_Flow.p_p0, label="p/p0")
    plt.plot([math.degrees(x) for x in cone_Flow.theta], cone_Flow.rho_rho0, label="\u03C1/\u03C10")
    plt.grid(), plt.xlabel('\u03B8 (deg)'), plt.legend(), plt.show() 

    pass 