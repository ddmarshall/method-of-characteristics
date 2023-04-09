import numpy as np 
from scipy.optimize import fsolve 
#import scipy.optimize as scp_opt
import math

"""
contains various unit process subroutines for supersonic flow - method of characteristics computations
Methods are from Gas Dynamics Vol 2 by Zucrow & Hoffman Chapter 16 - the method of characteristics applied to steady two-dimensional irrotational supersonic flow 
"""
class operator_funcs: 
    """
    generates repeatedly-used functions object which are necessary for MOC unit processes
    """
    def __init__(self):
        self.a = lambda a0, gam, u, v : math.sqrt(a0**2 - 0.5*(gam-1)*(u**2 + v**2))
        self.S = lambda delta, a, v, y : delta*(a**2*v/y)
        self.Q = lambda u, a : u**2 - a**2
        self.R = lambda u, v, Q, lam: 2*u*v - Q*lam
        self.lam = lambda lam1, lam2 : 0.5*(lam1 + lam2)
        self.lam_min = lambda u, v, a : (u*v - a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
        self.lam_plus = lambda u, v, a : (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)


def get_percent_changes(pt_old, pt_new):
    """
    Evaluates percent change in velocity and position over one iteration
    """
    vel_old = math.sqrt(pt_old[2]**2 + pt_old[3]**2)
    pos_old = math.sqrt(pt_old[0]**2 + pt_old[1]**2)
    delVel = abs(math.sqrt(pt_old[2]**2 + pt_old[3]**2) - math.sqrt(pt_new[2]**2 + pt_new[3]**2))#change in velocity
    delPos = abs(math.sqrt(pt_old[0]**2 + pt_old[1]**2) - math.sqrt(pt_new[0]**2 + pt_new[1]**2))#change in position

    pcVel = delVel/vel_old
    pcPos = delPos/pos_old

    return pcVel, pcPos


def interior_point(pt2, pt1, gasProps, delta, pcTOL, funcs):
    """
    #Interior point using euler predictor-corrector method on steady, axisymmetric, irrotational supersonic flow
    #Using a^2 = gamma*R*T - (gamma-1)*V^2/2 for speed of sound 
    pt1: above point (source of negative characteristic)
    pt2: below point (source of positive characteristic)
    """
    #Unpacking Inputs
    x1,x2 = pt1.x, pt2.x
    y1,y2 = pt1.y, pt2.y
    u1,u2 = pt1.u, pt2.u
    v1,v2 = pt1.v, pt2.v

    gam = gasProps.gam 
    R = gasProps.R
    T0 = gasProps.T0
    a0 = math.sqrt(gam*R*T0)

    #Calculate Coefficients for the Predictor 
    u_p, v_p, y_p = u2, v2, y2
    u_m, v_m, y_m = u1, v1, y1

    a_p = funcs.a(a0, gam, u_p, v_p)    #+ characteristic
    a_m = funcs.a(a0, gam, u_m, v_m)    #- characteristic

    #+characteristic
    lam_p = funcs.lam_plus(u_p, v_p, a_p)
    Q_p = funcs.Q(u_p, a_p)
    R_p = funcs.R(u_p, v_p, Q_p, lam_p)
    S_p = funcs.S(delta, a_p, v_p, y_p)

    #-characteristic
    lam_m = funcs.lam_min(u_m, v_m, a_m)
    Q_m = funcs.Q(u_m, a_m)
    R_m = funcs.R(u_m, v_m, Q_m, lam_m)
    S_m = funcs.S(delta, a_m, v_m, y_m)

    #determination of x4, y4, u4, and v4 for predictor 
    [y4, x4] = np.linalg.solve(np.array([[1,-lam_p],[1,-lam_m]]), np.array([y2 - lam_p*x2, y1 - lam_m*x1]))

    T_p = S_p*(x4 - x2) + Q_p*u2 + R_p*v2
    T_m = S_m*(x4 - x1) + Q_m*u1 + R_m*v1
    
    #getting velocity at pt4
    [u4, v4] = np.linalg.solve(np.array([[Q_p, R_p],[Q_m, R_m]]), np.array([T_p, T_m]))

    #Iteration of corrector: 
    def corrector(y1, y2, y4, u1, u2, u4, v1, v2, v4):

        u_p, v_p, y_p = 0.5*(u2+u4), 0.5*(v2+v4), 0.5*(y2+y4)
        u_m, v_m, y_m = 0.5*(u1+u4), 0.5*(v1+v4), 0.5*(y1+y4)

        a_p = funcs.a(a0, gam, u_p, v_p)#+ characteristic
        a_m = funcs.a(a0, gam, u_m, v_m)#- characteristic

        #+characteristic
        lam_p = funcs.lam_plus(u_p, v_p, a_p)
        Q_p = funcs.Q(u_p, a_p)
        R_p = funcs.R(u_p, v_p, Q_p, lam_p)
        S_p = funcs.S(delta, a_p, v_p, y_p)

        #-characteristic
        lam_m = funcs.lam_min(u_m, v_m, a_m)
        Q_m = funcs.Q(u_m, a_m)
        R_m = funcs.R(u_m, v_m, Q_m, lam_m)
        S_m = funcs.S(delta, a_m, v_m, y_m)

        #determination of x4, y4, u4, and v4 for predictor 
        [y4, x4] = np.linalg.solve(np.array([[1,-lam_p],[1,-lam_m]]), np.array([y2 - lam_p*x2, y1 - lam_m*x1]))
        T_p = S_p*(x4 - x2) + Q_p*u2 + R_p*v2
        T_m = S_m*(x4 - x1) + Q_m*u1 + R_m*v1
        
        #getting velocity at pt4
        [u4, v4] = np.linalg.solve(np.array([[Q_p, R_p],[Q_m, R_m]]), np.array([T_p, T_m]))

        return [u4, v4, x4, y4]
    
    pc_it = pcTOL
    while pc_it >= pcTOL:
        x4_old, y4_old = x4, y4
        u4_old, v4_old = u4, v4
        [u4, v4, x4, y4] = corrector(y1, y2, y4, u1, u2, u4, v1, v2, v4)
        pcVel, pcPos = get_percent_changes([x4_old, y4_old, u4_old, v4_old],[x4, y4, u4, v4])
        pc_it = max([pcVel, pcPos])

    return [x4, y4, u4, v4]


def direct_wall(pt2, y_x, dydx, gasProps, delta, pcTOL, funcs, charDir):
    
    #unpacking input parameters
    x2 = pt2.x
    y2 = pt2.y
    u2 = pt2.u
    v2 = pt2.v
    gam = gasProps.gam
    R = gasProps.R
    T0 = gasProps.T0 
    a0 = math.sqrt(gam*R*T0)

    #calculate flow properties at the initial point (assuming boundary above)
    a2 = funcs.a(a0, gam, u2, v2)
    
    #determining coefficients for the predictor
    if charDir == "pos":  
        lam = funcs.lam_plus(u2, v2, a2)
    elif charDir == "neg":
        lam = funcs.lam_min(u2, v2, a2)
     
    
    Q = funcs.Q(u2, a2)
    R = funcs.R(u2, v2, Q, lam)
    S = funcs.S(delta, a2, v2, y2)

    #Determination of x4, y4 for predictor
    func = lambda x4: y_x(x4) - lam*x4 - y2 + lam*x2

    x4 = fsolve(func, x0=x2)[0]
    y4 = y_x(x4)

    #getting u4 and v4 velocity for predictor
    T = S*(x4-x2) + Q*u2 + R*v2
    u4 = T/(Q + R*dydx(x4))
    v4 = dydx(x4)*u4        
    
    def corrector(x2, y2, u2, v2, x4, y4, u4, v4):
        #getting input values by averaging between pts 2 and 4
        x24, y24 = 0.5*(x2 + x4), 0.5*(y2 + y4)
        u24, v24 = 0.5*(u2 + u4), 0.5*(v2 + v4)

        #calculate coefficients for corrector 
        a24 = funcs.a(a0, gam, u24, v24)
        
        #determining coefficients for the corrector
        if charDir == "pos": 
            lam = funcs.lam_plus(u24, v24, a24)
        elif charDir == "neg":
            lam = funcs.lam_min(u24, v24, a24)

        Q = funcs.Q(u24, a24)
        R = funcs.R(u24, v24, Q, lam)
        S = funcs.S(delta, a24, v24, y24)
        
        #Determination of x4, y4
        x4 = fsolve(func, x0=x2)[0]
        y4 = y_x(x4)

        #getting u4 and v4 velocity for corrector
        T_p = S*(x4-x2) + Q*u2 + R*v2
        u4 = T_p/(Q + R*dydx(x4))
        v4 = dydx(x4)*u4 

        return [x4, y4, u4, v4]

    pc_it = pcTOL
    while pc_it >= pcTOL: 
        u4_old, v4_old = u4, v4
        x4_old, y4_old = x4, y4
        [x4, y4, u4, v4] = corrector(x2, y2, u2, v2, x4, y4, u4, v4)
        pcVel, pcPos = get_percent_changes([x4_old, y4_old, u4_old, v4_old],[x4, y4, u4, v4])
        pc_it = max([pcVel, pcPos])

    return [x4, y4, u4, v4]

"""
DEPRECATED FUNCTIONS (Uncomment to use)
def interior_point_old(pt1, pt2, gasProps, v_TOL): 
    
    #Interior point using euler predictor-corrector method on steady, axisymmetric, irrotational supersonic flow
    #Using a^2 = gamma*R*T - (gamma-1)*V^2/2 for speed of sound 
    #pt1: above point (source of negative characteristic)
    #pt2: below point (source of positive characteristic)
    
    #Unpacking Inputs
    x1,x2 = pt1.x, pt2.x
    y1,y2 = pt1.y, pt2.y
    u1,u2 = pt1.u, pt2.u
    v1,v2 = pt1.v, pt2.v

    gam = gasProps.gam 
    R = gasProps.R
    T = gasProps.T0 

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


def dir_wall_above_old(pt2, wall, wall_xBounds, wall_dydx, gasProps, v_TOL):
    
    #direct wall point calculation, assuming wall above initial point, using Euler predictor-corrector algorithm
    #pt_i: object containing intial point position and velocity components 
    #wall: function defining wall surface (give x, return y in [m])
    #wall_xBounds: iterable containing minimum and maximum x-values of surface (x_min, x_max)
    #wall_dydx: function defining wall slope (give x, return slope)
    #gasProps: object containing relevant thermodynamic properties of the gas 
    #v_TOL: velocity tolerance (m/s) for convergence evaluation 
    

    #speed of sound equation
    def a(gamma, R, T, V):
        return math.sqrt(gamma*R*T - V**2*(gamma-1)/2)

    #unpacking input parameters
    x2 = pt2.x
    y2 = pt2.y
    u2 = pt2.u
    v2 = pt2.v
    gam = gasProps.gam
    R = gasProps.R
    T = gasProps.T0 

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
        u4_old = u4
        v4_old = v4
        x4_old = x4
        y4_old = y4
        [x4, y4, u4, v4] = corrector(x2, y2, u2, v2, x4, y4, u4, v4, gam, R, T)

        TOL = max([abs(u4-u4_old), abs(v4-v4_old)])

    return [x4, y4, u4, v4]


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
"""

if __name__ == "__main__":
    #execute the following code when this file is run directly 
    class pointData: 
        def __init__(self, u, v, x, y):
            self.u, self.v, self.x, self.y = u, v, x, y

    class gasProps: 
        def __init__(self, R_gas, T0, gam):
            self.R, self.T0, self.gam = R_gas, T0, gam

    gasProps = gasProps(320, 3000, 1.2)
    funcs = operator_funcs() 
    delta = 1 #axisymmetric 

    #interior point solution: 
    pt1 = pointData(2473.4, 812.8, 0.13146, 0.040118)
    pt2 = pointData(2502.8, 737.6, 0.135683, 0.037123)
    [x4, y4, u4, v4] = interior_point(pt2, pt1, gasProps, delta, 0.00001, funcs)
    print(f"interior solution: u4,v4,x4,y4 = {u4, v4, x4, y4}")
    

    #direct wall above solution: 
    pt2 = pointData(1967.1, 1141.3, 0.06048, 0.059625)
    wall = lambda x : 0.0221852 + 0.71568*x - 1.0787*x**2 #wall y function 
    wall_dydx = lambda x : 0.71568 - 2*1.0787*x #wall slope function 
    wall_xBounds = (0, 0.25) #wall minimum and maximum x (m)
    #[x4, y4, u4, v4] = dir_wall_above_old(pt2, wall, wall_xBounds, wall_dydx, gasProps, 0.0001) 
    #print(f"old direct wall solution: u4,v4,x4,y4 = {u4, v4, x4, y4}")
    [x4, y4, u4, v4] = direct_wall(pt2, wall, wall_dydx, gasProps, delta, 0.000001, funcs, "pos")
    print(f"above direct wall solution: u4,v4,x4,y4 = {u4, v4, x4, y4}") 

    #direct wall below solution: 
    pt2.y, pt2.v = pt2.y*-1, pt2.v*-1
    wall_mod = lambda x: wall(x)*-1
    wall_dydx_mod = lambda x: wall_dydx(x)*-1
    [x4, y4, u4, v4] = direct_wall(pt2, wall_mod, wall_dydx_mod, gasProps, delta, 0.000001, funcs, "neg")
    print(f"below direct wall solution: u4,v4,x4,y4 = {u4, v4, x4, y4}") 