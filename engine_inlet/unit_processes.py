import math
import numpy as np

"""
Testing out Marshall's system of equations approach to solving for the flow field properties of an interior point
Assumptions: supersonic, irrotational, isentropic flow 
"""

class operator_functions: 
    """
    returns repeatedly-used functions which are necessary for MOC unit processes 
    """
    def __init__(self):
        self.a0 = lambda gam, R, T0: math.sqrt(gam*R*T0)
        self.a = lambda a0, gam, u, v : math.sqrt(a0**2 - 0.5*(gam-1)*(u**2 + v**2))
        self.S = lambda delta, a, v, y : delta*(a**2*v/y)
        self.Q = lambda u, a : u**2 - a**2
        self.R = lambda u, v, Q, lam: 2*u*v - Q*lam
        self.lam = lambda lam1, lam2 : 0.5*(lam1 + lam2)
        self.lam_min = lambda u, v, a : (u*v - a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
        self.lam_plus = lambda u, v, a : (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)

def interior_point(pt1, pt2, gasProps, delta, vel_TOL, funcs):
    """
    MOC interior point solution using irrotational, isentropic equations
    """
    #Unpack input objects
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y #point 1 data
    u2, v2, x2, y2 = pt2.u, pt2.v, pt2.x, pt2.y #point 2 data
    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties

    #initial estimates for first iteration: 
    u13, v13 = u1, v1 
    y13 = 1 if y1 == 0 else y1  
    u23, v23 = u2, v2
    y23 = 1 if y2 == 0 else y2

    #Calculate coefficients
    a0 = funcs.a0(gam, R_gas, T0) #stagnation speed of sound
    a1 = funcs.a(a0, gam, u1, v1)
    lam13 = funcs.lam_plus(u1, v1, a1)
    a13 = funcs.a(a0, gam, u13, v13)
    S13, Q13 = funcs.S(delta, a13, v13, y13), funcs.Q(u13, a13)
    R13 = funcs.R(u13, v13, Q13, lam13)

    a2 = funcs.a(a0, gam, u2, v2)
    lam23 = funcs.lam_min(u2, v2, a2)
    a23 = funcs.a(a0, gam, u23, v23)
    S23, Q23 = funcs.S(delta, a23, v23, y23), funcs.Q(u23, a23)
    R23 = funcs.R(u23, v23, Q23, lam23)

    #Construct system and invert to find solution vector 
    coeffMat = np.array([[lam13, -1, 0, 0],[lam23, -1, 0, 0],[-S13, 0, Q13, R13],[-S23, 0, Q23, R23]])
    RHSvec = np.array([lam13*x1-y1, lam23*x2-y2, -S13*x1+Q13*u1+R13*v1, -S23*x2+Q23*u2+R23*v2])
    [x3,y3,u3,v3] = np.linalg.solve(coeffMat, RHSvec)
    
    #Iterate above process until values converge
    def solve_interior_point(y3, u3, v3): 

        #Calculate Coefficients 
        u13, v13, y13 = 0.5*(u1 + u3), 0.5*(v1 + v3), 0.5*(y1 + y3)

        a3 = funcs.a(a0, gam, u3, v3)
        lam13 = funcs.lam(funcs.lam_plus(u1, v1, a1), funcs.lam_plus(u3, v3, a3))
        a13 = funcs.a(a0, gam, u13, v13)
        S13, Q13 = funcs.S(delta, a13, v13, y13), funcs.Q(u13, a13)
        R13 = funcs.R(u13, v13, Q13, lam13)

        u23, v23, y23 = 0.5*(u2 + u3), 0.5*(v2 + v3), 0.5*(y2 + y3)

        lam23 = funcs.lam(funcs.lam_min(u2, v2, a2), funcs.lam_min(u3, v3, a3))
        a23 = funcs.a(a0, gam, u23, v23)
        S23, Q23 = funcs.S(delta, a23, v23, y23), funcs.Q(u23, a23)
        R23 = funcs.R(u23, v23, Q23, lam23)

        #Construct System 
        coeffMat = np.array([[lam13, -1, 0, 0],[lam23, -1, 0, 0],[-S13, 0, Q13, R13],[-S23, 0, Q23, R23]])
        RHSvec = np.array([lam13*x1-y1, lam23*x2-y2, -S13*x1+Q13*u1+R13*v1, -S23*x2+Q23*u2+R23*v2])
        
        #Solve System and return solution vector
        return np.linalg.solve(coeffMat, RHSvec)

    delta_vel = vel_TOL
    while delta_vel >= vel_TOL: 
        #print("\niterating...\n")
        y3_old, u3_old, v3_old = y3, u3, v3
        [x3, y3, u3, v3] = solve_interior_point(y3_old, u3_old, v3_old)
        delta_vel = max([abs(u3 - u3_old), abs(v3 - v3_old)])

    return [x3, y3, u3, v3]

"""
def direct_wall(pt1, y, dydx, gasProps, delta, funcs):
    #MOC direct wall solution using irrotational, isentropic equations 
    #currently only works for wall above 
    
    #unpacking input data 
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y
    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties 

    #initial value computations 
    u13, v13 = u1, v1 
    y13 = 1 if y1 == 0 else y1

    a0 = funcs.a0(gam, R_gas, T0) #stagnation speed of sound

    a1 = funcs.a(a0, gam, u1, v1)
    a13 = funcs.a(a0, gam, u13, v13)
    lam13 = funcs.lam_plus(u1, v1, a1)
    S13 = funcs.S(delta, a13, v13, y13)
    Q13 = funcs.Q(u13, a13)
    R13 = funcs.R(u13, v13, Q13, lam13)

    #computing the x&y location of 3' for first iteration (2 methods suggested)
        #method 1: go straight up from x1 to get x_3' 
        x3p = x1
        y3p = #?Hmmm

        #method 2: use lam13 to project to wall to get (x,y)_3'

    #computing partial derivatives for taylor series expansion using 2nd order centered differences
    pfpx_3p =  (f(x+delx, y)-f(x-delx,y))/(2*delx) #partial f_w / partial x at 3' 
    pfpy_3p = (f(x, y+dely)-f(x, y-dely))/(2*dely) #partial f_w / partial y at 3'

    #construct system and solve
    coeffMat = np.array([[lam13, -1, 0, 0],[pfpx_3p, pfpy_3p, 0, 0],[-S13, 0, Q13, R13],[0, 0, pfpx_3p, pfpy_3p]])
    RHSvec = np.array([lam13*x1-y1, pfpx_3p*x3p+pfpy_3p*y3p-f_w3p, -S13*x1+Q13*u1+R13*v1, 0])
    [x3,y3,u3,v3] = np.linalg.solve(coeffMat, RHSvec)
"""
"""
def symmetry_boundary(pt2, gasProps, delta, funcs):
    #MOC symmetry boundary solution for irrotational, isentropic equations
    #currently only works for above symmetry plane

    #unpack inputs: 
    u2, v2, x2, y2 = pt2.u, pt2.v, pt2.x, pt2.y #point 2 data
    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties
    
    #make initial estimates for first iteration
    u23, v23 = u2, v2
    y23 = 1 if y2 == 0 else y2

    a0 = funcs.a0(gam, R_gas, T0) 
    a2 = funcs.a(a0, gam, u2, v2)
    a23 = funcs.a(a0, gam, u23, v23)
    lam23 = funcs.lam_min(u2, v2, a2)
    S23 = funcs.S(delta, a23, v23, y23)
    Q23 = funcs.Q(u23, a23)
    R23 = funcs.R(u23, v23, Q23, lam23)

    #xonstruct system and solve
    coeffMat = np.array([[lam23, -1, 0, 0],[0, 1, 0, 0],[-S23, 0, Q23, R23],[0, 0, 0, 1]])
    RHSvec = np.array([lam23*x2-y2, 0, -S23*x2+Q23*u2+R23*v2, 0])
    [x3,y3,u3,v3] = np.linalg.solve(coeffMat, RHSvec)

    #iterate until convergence is attained 
    #TODO
"""

if __name__  == "__main__":
    #Initial Values: 
    class point_data: 
        def __init__(self, u, v, x, y):
            self.u, self.v, self.x, self.y = u, v, x, y
    class gasProps:
        def __init__(self, gam, R, T0): 
            self.gam, self.R, self.T0 = gam, R, T0
    #Point 2 (m/s, m)
    pt2 = point_data(2473.4, 812.8, 0.13146, 0.040118)
   
    #point 1 (m/s, m)
    pt1 = point_data(2502.8, 737.6, 0.135683, 0.037123)
    
    #initial settings
    delta = 1 #axisymmetric flow
    gas = gasProps(1.2, 320, 3000)

    moc_funcs = operator_functions() #loading in heavily used functions in moc calculations

    [x3, y3, u3, v3] = interior_point(pt1, pt2, gas, delta, 0.001, moc_funcs)

pass