import numpy as np 
import math
import matplotlib.pyplot as plt
"""
Running some tests to compare Marshall's iterative system solve and Zucrow & Hoffman's Euler predictor-corrector 
for an interior point solution. Currently looks like Marshall's method converges on different values when inputs 
are flipped over the horizontal axis (differing in magnitude, not just in sign as expected). Not sure what the issue
is. 
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

def interior_point_mod(pt1, pt2, gasProps, delta, pcTOL, funcs):
    """
    MOC interior point solution using irrotational, isentropic equations
    """
    #Unpack input objects
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y #point 1 data
    u2, v2, x2, y2 = pt2.u, pt2.v, pt2.x, pt2.y #point 2 data
    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties

    #getting speed of sound 
    a0 = gasProps.a0 #stagnation speed of sound

    #initial estimates for first iteration: 
    u13, v13 = u1, v1 
    y13 = 1 if y1 == 0 else y1  
    u23, v23 = u2, v2
    y23 = 1 if y2 == 0 else y2

    uList, vList = [],[]
    xList, yList = [],[]

    def solve_interior_point(u13, v13, y13, u23, v23, y23, first_iter=None): 
        #Calculate coefficients
        a1 = funcs.a(a0, gam, u1, v1)

        if first_iter: 
            lam13 = funcs.lam_plus(u1, v1, a1)
        else: 
            a3 = funcs.a(a0, gam, u3, v3)
            lam13 = 0.5*(funcs.lam_plus(u1, v1, a1) + funcs.lam_plus(u3, v3, a3))

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
        return np.linalg.solve(coeffMat, RHSvec) #[x3, y3, u3, v3]
    
    #first iteration
    [x3, y3, u3, v3] = solve_interior_point(u13, v13, y13, u23, v23, y23, first_iter=True)
    uList.append(u3), vList.append(v3)
    xList.append(x3), yList.append(y3)

    #Iterate above process until values converge
    pc_it = pcTOL
    while pc_it >= pcTOL:
        x3_old, y3_old, u3_old, v3_old = x3, y3, u3, v3
        u13, v13, y13 = 0.5*(u1 + u3), 0.5*(v1 + v3), 0.5*(y1 + y3)
        u23, v23, y23 = 0.5*(u2 + u3), 0.5*(v2 + v3), 0.5*(y2 + y3)
        [x3, y3, u3, v3] = solve_interior_point(u13, v13, y13, u23, v23, y23)
        uList.append(u3), vList.append(v3)
        xList.append(x3), yList.append(y3)

        pcVel, pcPos = get_percent_changes([x3_old, y3_old, u3_old, v3_old],[x3, y3, u3, v3]) #get percent change across iteration
        pc_it = max([pcVel, pcPos])

    return [x3, y3, u3, v3, xList, yList, uList, vList]

def interior_point_ZH(pt1, pt2, gasProps, v_TOL): 
    
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
    xList, yList = [x4],[y4]

    T_p = S_p*(x4 - x2) + Q_p*u2 + R_p*v2
    T_m = S_m*(x4 - x1) + Q_m*u1 + R_m*v1
    
    #getting velocity at pt4
    [u4, v4] = np.linalg.solve(np.array([[Q_p, R_p],[Q_m, R_m]]), np.array([T_p, T_m]))
    uList, vList = [u4], [v4]

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
        xList.append(x4), yList.append(y4), uList.append(u4), vList.append(v4)
        TOL = max([abs(u4-u4_old), abs(v4-v4_old)])
    
    return [x4, y4, u4, v4, xList, yList, uList, vList]

def interior_convergence_study(pt1s, pt2s, gasProps, delta, pcTOL, funcs): 

    if len(pt1s) != len(pt2s):
        raise ValueError("length of pt1s and pt2s must be the same")
    
    fig = plt.figure()
    xPlot = fig.add_subplot(221)
    yPlot = fig.add_subplot(222)
    uPlot = fig.add_subplot(223)
    vPlot = fig.add_subplot(224)
    for i, pt1 in enumerate(pt1s): 
        pt2 = pt2s[i]
        [_,_,_,_,xList,yList,uList,vList] = interior_point_mod(pt1, pt2, gasProps, delta, pcTOL, funcs)
        xPlot.plot(range(len(xList)), [abs(x) for x in xList], '-o', label=str(i))
        yPlot.plot(range(len(yList)), [abs(x) for x in yList], '-o', label=str(i))
        uPlot.plot(range(len(uList)), [abs(x) for x in uList], '-o', label=str(i))
        vPlot.plot(range(len(vList)), [abs(x) for x in vList], '-o', label=str(i))

    xPlot.set_ylabel("x-pos (m)"), yPlot.set_ylabel("absolute y-pos (m)")
    uPlot.set_ylabel("x-vel (m/s)"), vPlot.set_ylabel("absolute y-vel (m/s)")
    xPlot.grid(), yPlot.grid(), uPlot.grid(), vPlot.grid() 
    xPlot.legend()
    
    plt.show()

def interior_convergence_study_ZH(pt1s, pt2s, gasProps, v_TOL):
    
    if len(pt1s) != len(pt2s):
        raise ValueError("length of pt1s and pt2s must be the same")
    
    fig = plt.figure()
    xPlot = fig.add_subplot(221)
    yPlot = fig.add_subplot(222)
    uPlot = fig.add_subplot(223)
    vPlot = fig.add_subplot(224)
    for i, pt1 in enumerate(pt1s): 
        pt2 = pt2s[i]
        [_,_,_,_,xList,yList,uList,vList] = interior_point_ZH(pt2, pt1, gasProps, v_TOL)
        xPlot.plot(range(len(xList)), [abs(x) for x in xList], '-o', label=str(i))
        yPlot.plot(range(len(yList)), [abs(x) for x in yList], '-o', label=str(i))
        uPlot.plot(range(len(uList)), [abs(x) for x in uList], '-o', label=str(i))
        vPlot.plot(range(len(vList)), [abs(x) for x in vList], '-o', label=str(i))

    xPlot.set_ylabel("x-pos (m)"), yPlot.set_ylabel("absolute y-pos (m)")
    uPlot.set_ylabel("x-vel (m/s)"), vPlot.set_ylabel("absolute y-vel (m/s)")
    xPlot.grid(), yPlot.grid(), uPlot.grid(), vPlot.grid() 
    xPlot.legend()
    
    plt.show()


if __name__ == "__main__":
    import copy
    class gasProps: 
        def __init__(self, gam, R, T0):
            self.gam, self.R, self.T0 = gam, R, T0
            self.a0 = math.sqrt(gam*R*T0)
    class pointData:
        def __init__(self, u, v, x, y, thet=None):
            self.u, self.v, self.x, self.y = u, v, x, y
            if thet is not None: self.thet = thet

    gas = gasProps(1.2, 320, 3000)  
    funcs = operator_funcs()
    delta = 1

    #Testing Interior Operator 
    pt2_1 = pointData(2473.4, 812.8, 0.13146, 0.040118) #Point 2 (m/s, m)
    pt1_1 = pointData(2502.8, 737.6, 0.135683, 0.037123) #point 1 (m/s, m

    pt2_2 = copy.copy(pt1_1)
    pt1_2 = copy.copy(pt2_1)

    pt2_2.v, pt2_2.y = pt2_2.v*-1, pt2_2.y*-1
    pt1_2.v, pt1_2.y = pt1_2.v*-1, pt1_2.y*-1

    pt1s = [pt1_1, pt1_2]
    pt2s = [pt2_1, pt2_2]

    #interior_convergence_study(pt1s, pt2s, gas, delta, 0.0001, funcs)
    interior_convergence_study_ZH(pt1s, pt2s, gas, 0.001)