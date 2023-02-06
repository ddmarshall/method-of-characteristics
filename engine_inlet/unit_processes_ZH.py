import numpy as np 
import scipy.optimize as scp_opt 
import math

"""
contains various unit process subroutines for supersonic flow - method of characteristics computations
Methods are from Gas Dynamics Vol 2 by Zucrow & Hoffman Chapter 16 - the method of characteristics applied to steady two-dimensional irrotational supersonic flow 
"""
#Interior Operator-------------------------------------------------------------------------------------------------------- 
def interior(pt1, pt2, gasProps, v_TOL): 
    
    #Interior point using euler predictor-corrector method on steady, axisymmetric, irrotational supersonic flow
    #Using a^2 = gamma*R*T - (gamma-1)*V^2/2 for speed of sound 

    #Unpacking Inputs
    x1 = pt1.x
    x2 = pt2.x
    y1 = pt1.y
    y2 = pt2.y
    u1 = pt1.u
    u2 = pt2.u
    v1 = pt1.v
    v2 = pt2.v

    gam = gasProps.gamma 
    R = gasProps.R
    T = gasProps.T_0 

    #speed of sound equation (perfect gas)
    def a(gamma, R, T, V):
        return math.sqrt(gamma*R*T - V**2*(gamma-1)/2)

    #Calculate Coefficients for the Predictor 
    u_p = u2
    v_p = v2
    y_p = y2

    u_m = u1
    v_m = v1
    y_m = y1

    #+ characteristic
    V_p = math.sqrt(u_p**2 + v_p**2)
    thet_p = math.atan(v_p/u_p)
    a_p = a(gam, R, T, V_p)
    alf_p = math.asin(a_p/V_p)

    #- characteristic
    V_m = math.sqrt(u_m**2 + v_m**2)
    thet_m = math.atan(v_m/u_m)
    a_m = a(gam, R, T, V_m)
    alf_m = math.asin(a_m/V_m)

    #+characteristic
    lam_p = math.tan(thet_p + alf_p)
    Q_p = u_p**2 - a_p**2 
    R_p = 2*u_p*v_p - Q_p*lam_p
    S_p = a_p**2*v_p/y_p

    #-characteristic
    lam_m = math.tan(thet_m - alf_m)
    Q_m = u_m**2 - a_m**2 
    R_m = 2*u_m*v_m - Q_m*lam_m
    S_m = a_m**2*v_m/y_m

    #determination of x4, y4, u4, and v4 for predictor 

    [y4, x4] = np.linalg.solve(np.array([[1,-lam_p],[1,-lam_m]]), np.array([y2 - lam_p*x2, y1 - lam_m*x1]))

    T_p = S_p*(x4 - x2) + Q_p*u2 + R_p*v2
    T_m = S_m*(x4 - x1) + Q_m*u1 + R_m*v1
    
    #getting velocity at pt4
    [u4, v4] = np.linalg.solve(np.array([[Q_p, R_p],[Q_m, R_m]]), np.array([T_p, T_m]))

    #Iteration of corrector: 
    def corrector(y1, y2, y4, u1, u2, u4, v1, v2, v4):

        u_p = 0.5*(u2+u4)
        v_p = 0.5*(v2+v4)
        y_p = 0.5*(y2+y4)
        u_m = 0.5*(u1+u4)
        v_m = 0.5*(v1+v4)
        y_m = 0.5*(y1+y4)

        #+ characteristic
        V_p = math.sqrt(u_p**2 + v_p**2)
        thet_p = math.atan(v_p/u_p)
        a_p = a(gam, R, T, V_p)
        alf_p = math.asin(a_p/V_p)

        #- characteristic
        V_m = math.sqrt(u_m**2 + v_m**2)
        thet_m = math.atan(v_m/u_m)
        a_m = a(gam, R, T, V_m)
        alf_m = math.asin(a_m/V_m)

        #+characteristic
        lam_p = math.tan(thet_p + alf_p)
        Q_p = u_p**2 - a_p**2 
        R_p = 2*u_p*v_p - Q_p*lam_p
        S_p = a_p**2*v_p/y_p

        #-characteristic
        lam_m = math.tan(thet_m - alf_m)
        Q_m = u_m**2 - a_m**2 
        R_m = 2*u_m*v_m - Q_m*lam_m
        S_m = a_m**2*v_m/y_m

        #determination of x4, y4, u4, and v4 for predictor 

        [y4, x4] = np.linalg.solve(np.array([[1,-lam_p],[1,-lam_m]]), np.array([y2 - lam_p*x2, y1 - lam_m*x1]))

        T_p = S_p*(x4 - x2) + Q_p*u2 + R_p*v2
        T_m = S_m*(x4 - x1) + Q_m*u1 + R_m*v1
        
        #getting velocity at pt4
        [u4, v4] = np.linalg.solve(np.array([[Q_p, R_p],[Q_m, R_m]]), np.array([T_p, T_m]))

        return [u4, v4, x4, y4]
    
    TOL = v_TOL
    while TOL >= v_TOL:
        print("iterating...")
        u4_old = u4
        v4_old = v4
        [u4, v4, x4, y4] = corrector(y1, y2, y4, u1, u2, u4, v1, v2, v4)
        TOL = max([abs(u4-u4_old), abs(v4-v4_old)])
    
    return [u4, v4, x4, y4]

#Axis of Symmetry Operator-----------------------------------------------------------------------------------------------
def axis(pt1, gasProps, v_TOL):
    
    #Axis of Symmetry Point Calculation using Euler Predictor Corrector Method
    #pt1: object containing x,y,u,v components
    #gasProps: object containing relevant thermodynamic gas properties
    
    #Unpacking Inputs
    x1 = pt1.x
    y1 = pt1.y
    u1 = pt1.u
    v1 = pt1.v

    gam = gasProps.gamma 
    R = gasProps.R
    T = gasProps.T_0 

    #speed of sound equation (perfect gas)
    def a(gamma, R, T, V):
        return math.sqrt(gamma*R*T - V**2*(gamma-1)/2)

    #calculate flow properties at initial-value point
    V1 = math.sqrt(u1**2 + v1**2)
    thet1 = math.atan(v1/u1)
    a1 = a(gam, R, T, V1)
    alf1 = math.asin(a1/V1)

    #calculation of the coefficients for the predictor
    lam_m = math.tan(thet1 - alf1)
    Q_m = u1**2 - a1**2
    R_m = 2*u1*v1 - (u1**2 - a1**2)*lam_m
    delta = 1 # set to 0 for non axisymmetric 
    S_m = delta*a1**2*v1/y1

    #determination of x4 and u4 for predictor algorithm
    y4 = 0 #y-location of axis of symmetry 
    x4 = (y1 - lam_m*x1)/(-lam_m) - y4

    #determination of T constant for predictor 
    T_m = S_m*(x4-x1) + Q_m*u1 + R_m*v1

    #determination of u4 noting that v4 = 0
    v4 = 0 
    u4 = T_m/Q_m

    def corrector(y1, y4, u1, u4, v1, v4): 

        u_m = 0.5*(u1+u4)
        v_m = 0.5*(v1+v4)
        y_m = 0.5*(y1+y4)

        #- characteristic
        V_m = math.sqrt(u_m**2 + v_m**2)
        thet_m = math.atan(v_m/u_m)
        a_m = a(gam, R, T, V_m)
        alf_m = math.asin(a_m/V_m)

        #-characteristic
        lam_m = math.tan(thet_m - alf_m)
        Q_m = u_m**2 - a_m**2 
        R_m = 2*u_m*v_m - Q_m*lam_m
        S_m = a_m**2*v_m/y_m

        y4 = 0 #y-location of axis of symmetry 
        x4 = (y1 - lam_m*x1)/(-lam_m) - y4

        T_m = S_m*(x4-x1) + Q_m*u1 + R_m*v1
        u4 = T_m/Q_m

        return [x4, u4]

    TOL = v_TOL
    while TOL >= v_TOL:
        print("iterating...")
        u4_old = u4
        [x4, u4] = corrector(y1, y4, u1, u4, v1, v4)
        TOL = abs(u4-u4_old)

    return [x4, u4]

#Direct Wall Point Operator---------------------------------------------------------------------------------------------- 
def dir_wall_above(pt2, wall, wall_xBounds, wall_dydx, gasProps, v_TOL):
    
    #direct wall point calculation, assuming wall above initial point, using Euler predictor-corrector algorithm
    #pt_i: object containing intial point position and velocity components 
    #wall: function defining wall surface (give x, return y in [m])
    #wall_xBounds: iterable containing minimum and maximum x-values of surface (x_min, x_max)
    #wall_dydx: function defining wall slope (give x, return slope)
    #gasProps: object containing relevant thermodynamic properties of the gas 
    #v_TOL: velocity tolerance (m/s) for convergence evaluation 
    
    """
    #getting wall position and first derivative functions for wall boundary:
    def wall_y_interp():
        interpFunc = scp_int.interp1d(wall[0], wall[1], kind="cubic") 
        return interpFunc

    def wall_dydx(wall_x, interpFunc):
        #pretty jerry-rigged function (return to this when improving)
        x_step = (max(wall[0]) - min(wall[0]))/10000 #one thousandth of the x interval 
        
        y1 = interpFunc(wall_x - x_step)
        y2 = interpFunc(wall_x + x_step)
        return (y2 - y1)/(2*x_step)
    """

    #speed of sound equation
    def a(gamma, R, T, V):
        return math.sqrt(gamma*R*T - V**2*(gamma-1)/2)

    #unpacking input parameters
    x2 = pt2.x
    y2 = pt2.y
    u2 = pt2.u
    v2 = pt2.v
    gam = gasProps.gamma 
    R = gasProps.R
    T = gasProps.T_0 

    #calculate flow properties at the initial point (assuming boundary above)
    V2 = math.sqrt(u2**2 + v2**2)
    thet2 = math.atan(v2/u2)
    a2 = a(gam, R, T, V2)
    alf2 = math.asin(a2/V2)
    
    #determining coefficients for the predictor 
    lam_p = math.tan(thet2 + alf2)
    Q_p = u2**2 - a2**2
    R_p = 2*u2*v2 - (u2**2 - a2**2)*lam_p
    S_p = a2**2*v2/y2
    
    #Determination of x4, y4 for predictor
    def func(x4, wall, lam_p, y2, x2):
        RHS = y2 - lam_p*x2
        LHS = wall(x4) - lam_p*x4

        return abs(RHS - LHS)

    res = scp_opt.minimize(func, [0], args=(wall, lam_p, y2, x2), bounds=[(wall_xBounds[0], wall_xBounds[1])])
    
    x4 = res.x
    y4 = wall(x4)

    #getting u4 and v4 velocity for predictor
    T_p = S_p*(x4-x2) + Q_p*u2 + R_p*v2
    m4 = wall_dydx(x4) #getting wall slope at point 4
    u4 = T_p/(Q_p + R_p*m4)
    v4 = m4*u4        
    
    def corrector(x2, y2, u2, v2, x4, y4, u4, v4, gam, R, T):
        #getting input values by averaging between pts 2 and 4
        x2 = 0.5*(x2 + x4)
        y2 = 0.5*(y2 + y4)
        u2 = 0.5*(u2 + u4)
        v2 = 0.5*(v2 + v4)

        #calculate coefficients for corrector 
        V2 = math.sqrt(u2**2 + v2**2)
        thet2 = math.atan(v2/u2)
        a2 = a(gam, R, T, V2)
        alf2 = math.asin(a2/V2)
        
        #determining coefficients for the corrector 
        lam_p = math.tan(thet2 + alf2)
        Q_p = u2**2 - a2**2
        R_p = 2*u2*v2 - (u2**2 - a2**2)*lam_p
        S_p = a2**2*v2/y2
        
        #Determination of x4, y4
        res = scp_opt.minimize(func, [0], args=(wall, lam_p, y2, x2), bounds=[(wall_xBounds[0], wall_xBounds[1])])
        
        x4 = res.x
        y4 = wall(x4)

        #getting u4 and v4 velocity for corrector
        T_p = S_p*(x4-x2) + Q_p*u2 + R_p*v2
        m4 = wall_dydx(x4) #getting wall slope at point 4
        u4 = T_p/(Q_p + R_p*m4)
        v4 = m4*u4 

        return [x4, y4, u4, v4]

    TOL = v_TOL
    while TOL >= v_TOL: 
        print("iterating...")
        u4_old = u4
        v4_old = v4
        x4_old = x4
        y4_old = y4
        [x4, y4, u4, v4] = corrector(x2, y2, u2, v2, x4, y4, u4, v4, gam, R, T)

        TOL = max([abs(u4-u4_old), abs(v4-v4_old)])

    return [x4, y4, u4, v4]

#Inverse Wall Point Operator---------------------------------------------------------------------------------------------
def inv_wall_above(pt1, pt3, pt4, gasProps, v_TOL):

    #unpacking input parameters
    x1 = pt1.x
    y1 = pt1.y
    u1 = pt1.u
    v1 = pt1.v 

    x3 = pt3.x
    y3 = pt3.y
    u3 = pt3.u
    v3 = pt3.v 

    x4 = pt4.x
    y4 = pt4.y

    gam = gasProps.gamma 
    R = gasProps.R
    T = gasProps.T_0 

    #speed of sound equation 
    def a(gamma, R, T, V):
        return math.sqrt(gamma*R*T - V**2*(gamma-1)/2)
    
    #location of point 2 for the predictor 
    u2 = u3 #first pass assumption
    v2 = v3 #first pass assumption 

    thet2 = math.atan(v2/u2)
    V2 = math.sqrt(u2**2 + v2**2)
    a2 = a(gam, R, T, V2)
    alf2 = math.asin(a2/V2)
    lam_p = math.tan(thet2 + alf2)

    #solving for x2 and y2 
    m_13 = (y3-y1)/(x3-x1) #slope of line connecting points 1 and 3
    [x2, y2] = np.linalg.solve(np.array([[lam_p, -1],[m_13, -1]]), np.array([lam_p*x4 - y4, m_13*x1 - y1]))

    #interpolating along line 13 to get u2 and v2
    u2 = u1 + (x2-x1)/(x3-x1)*(u3-u1)
    v2 = v1 + (x2-x1)/(x3-x1)*(v3-v1)
    
    #calculation of the coefficients for the predictor with new values @ point 2 
    thet2 = math.atan(v2/u2)
    V2 = math.sqrt(u2**2 + v2**2)
    a2 = a(gam, R, T, V2)
    alf2 = math.asin(a2/V2)

    lam_p = math.tan(thet2 + alf2)
    Q_p = u2**2 - a2**2 
    R_p = 2*u2*v2 - Q_p*lam_p
    S_p = a2**2*v2/y2

    #determination of u4 and v4 for the predictor
    T_p = S_p*(x4-x2) + Q_p*u2 + R_p*v2
    thet4 = math.radians(25) #HARD CODED VALUE (remove when done testing)

    [u4, v4] = np.linalg.solve(np.array([[Q_p, R_p],[math.tan(thet4), -1]]), np.array([T_p, 0]))

    #application of the corrector
    def corrector(x1, y1, u2, v2, x3, x4, y4, u4, v4, gam, R ,T):
        
        #updating input velocity components between points 2 and 4
        u2 = 0.5*(u2 + u4)
        v2 = 0.5*(v2 + v4)

        #calculating flow properties for corrector
        thet2 = math.atan(v2/u2)
        V2 = math.sqrt(u2**2 + v2**2)
        a2 = a(gam, R, T, V2)
        alf2 = math.asin(a2/V2)
        lam_p = math.tan(thet2 + alf2)
        #solving for x2 and y2 
        [x2, y2] = np.linalg.solve(np.array([[lam_p, -1],[m_13, -1]]), np.array([lam_p*x4 - y4, m_13*x1 - y1]))

        #interpolating along line 13 to get u2 and v2
        u2 = u1 + (x2-x1)/(x3-x1)*(u3-u1)
        v2 = v1 + (x2-x1)/(x3-x1)*(v3-v1)
        
        #calculation of the coefficients for the corrector
        thet2 = math.atan(v2/u2)
        V2 = math.sqrt(u2**2 + v2**2)
        a2 = a(gam, R, T, V2)
        alf2 = math.asin(a2/V2)

        lam_p = math.tan(thet2 + alf2)
        Q_p = u2**2 - a2**2 
        R_p = 2*u2*v2 - Q_p*lam_p
        S_p = a2**2*v2/y2

        #determination of u4 and v4 for the corrector
        T_p = S_p*(x4-x2) + Q_p*u2 + R_p*v2
        thet4 = math.radians(25) #HARD CODED VALUE (remove when done testing)
        [u4, v4] = np.linalg.solve(np.array([[Q_p, R_p],[math.tan(thet4), -1]]), np.array([T_p, 0]))

        return [u4, v4]

    TOL = v_TOL
    while TOL >= v_TOL:
        print("Iterating...")
        u4_old = u4
        v4_old = v4
        [u4, v4] = corrector(x1, y1, u2, v2, x3, x4, y4, u4, v4, gam, R ,T)

        TOL = max([abs(u4-u4_old), abs(v4-v4_old)])

    return [u4, v4] 

if __name__ == "__main__":
    #execute the following code when this file is run directly 
    class pointData: 
        def __init__(self, u, v, x, y):
            self.u, self.v, self.x, self.y = u, v, x, y

    class gasProps: 
        def __init__(self, R_gas, T0, gam):
            self.R, self.T_0, self.gamma = R_gas, T0, gam

    pt2 = pointData(1967.1, 1141.3, 0.06048, 0.059625)
    gasProps = gasProps(320, 3000, 1.2)
    wall = lambda x : 0.0221852 + 0.71568*x - 1.0787*x**2 #wall y function 
    wall_dydx = lambda x : 0.71568 - 2*1.0787*x #wall slope function 
    wall_xBounds = (0, 0.25) #wall minimum and maximum x (m)
    [x4, y4, u4, v4] = dir_wall_above(pt2, wall, wall_xBounds, wall_dydx, gasProps, 0.001) 
    pass 