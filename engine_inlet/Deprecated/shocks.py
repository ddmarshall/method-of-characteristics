import math 
import scipy.optimize as sci_opt

"""
Various Shock Wave Subroutines for Inlet MOC code 
"""

#From Anderson, J.D.A. Modern Compressible Flow page 135-6
def oblique_shock_MCF(M1, thet, gam):
    #M      : upstream mach number 
    #thet   : flow deflection angle (rad) 
    #gam    : specific heat ratio 

    #returns: [downstream mach number, shock angle (rad)]

    #first calculate the shock angle from the theta-beta-M equation 
    func = lambda beta, M, thet, gam : 2*(1/math.tan(beta))*((M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta)) + 2)) - math.tan(thet)

    beta = sci_opt.fsolve(func, thet, args=(M1, thet, gam))

    #get the normal mach number component: 
    M_n1 = M1*math.sin(beta)

    #calculate normal mach number after shock using normal shock relation: 
    M_n2 = math.sqrt((M_n1**2 + (2/(gam-1)))/((2*gam/(gam-1))*M_n1**2 - 1))

    #get downstream mach number 
    M2 = M_n2/math.sin(beta-thet)

    #return results 
    return [M2, beta]

#[M2, beta] = oblique_shock_MCF(1.7, math.radians(15), 1.4)

def conical_oblique_shock(M1, c1, thet, gam):
    """
    Numerically solves for the initial data line following the oblique shock on a cone 
    Method from Anderson, B.H NASA technical note 1969
    """

    #estimate shock angle (beta_s) based on freestream conditions 

    b = -(M1**2 + 2)/M1**2 - gam*math.sin(thet)**2
    c = (2*M1**2 + 1)/M1**4 + (((gam+1)/2)**2 + ((gam-1)/M1**2))*math.sin(thet)**2
    d = -math.cos(thet)**2/M1**4

    func = lambda beta_s, b, c, d : math.sin(beta_s)**6 + b*math.sin(beta_s)**4 + c*math.sin(beta_s)**2 + d
    beta_s = sci_opt.fsolve(func, thet, args=(b,c,d))

    #Get flow properties immediately behind shock using oblique shock relations 

    sig2 = (gam-1)/(gam+1)
    q0 = M1*c1

    u = (1 - sig2)*q0*math.cos(beta_s)**2 + 1/q0
    v = (q0 - u)/math.tan(beta_s)
    dvdu = v/(u - q0)

    #Numerically integrate governing equation

    

pass 