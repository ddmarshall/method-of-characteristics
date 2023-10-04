import math
import numpy as np
import scipy.optimize

"""
Method of characteristics operator functions for irrotational, isentropic axisymmetric/2D flow. 
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


def interior_point(pt1, pt2, gasProps, delta, pcTOL, funcs):
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

    def solve_interior_point(u13, v13, y13, u23, v23, y23, first_iter=None): 
        #Calculate coefficients
        a1 = funcs.a(a0, gam, u1, v1)
        a2 = funcs.a(a0, gam, u2, v2)
        if first_iter: 
            lam13 = funcs.lam_plus(u1, v1, a1)
            lam23 = funcs.lam_min(u2, v2, a2)
        else: 
            a3 = funcs.a(a0, gam, u3, v3)
            lam13 = 0.5*(funcs.lam_plus(u1, v1, a1) + funcs.lam_plus(u3, v3, a3))
            lam23 = 0.5*(funcs.lam_min(u2, v2, a2) + funcs.lam_min(u3, v3, a3))

        a13 = funcs.a(a0, gam, u13, v13)
        S13, Q13 = funcs.S(delta, a13, v13, y13), funcs.Q(u13, a13)
        R13 = funcs.R(u13, v13, Q13, lam13)

        a23 = funcs.a(a0, gam, u23, v23)
        S23, Q23 = funcs.S(delta, a23, v23, y23), funcs.Q(u23, a23)
        R23 = funcs.R(u23, v23, Q23, lam23)

        #Construct system and invert to find solution vector 
        coeffMat = np.array([[lam13, -1, 0, 0],[lam23, -1, 0, 0],[-S13, 0, Q13, R13],[-S23, 0, Q23, R23]])
        RHSvec = np.array([lam13*x1-y1, lam23*x2-y2, -S13*x1+Q13*u1+R13*v1, -S23*x2+Q23*u2+R23*v2])
        return np.linalg.solve(coeffMat, RHSvec) #[x3, y3, u3, v3]
    
    #first iteration
    [x3, y3, u3, v3] = solve_interior_point(u13, v13, y13, u23, v23, y23, first_iter=True)

    #Iterate above process until values converge
    pc_it = pcTOL
    while abs(pc_it) >= pcTOL:
        x3_old, y3_old, u3_old, v3_old = x3, y3, u3, v3
        u13, v13, y13 = 0.5*(u1 + u3), 0.5*(v1 + v3), 0.5*(y1 + y3)
        u23, v23, y23 = 0.5*(u2 + u3), 0.5*(v2 + v3), 0.5*(y2 + y3)
        [x3, y3, u3, v3] = solve_interior_point(u13, v13, y13, u23, v23, y23)

        pcVel, pcPos = get_percent_changes([x3_old, y3_old, u3_old, v3_old],[x3, y3, u3, v3]) #get percent change across iteration
        pc_it = max([pcVel, pcPos])

    return [x3, y3, u3, v3]


def direct_wall(pt1, y_x, dydx, gasProps, delta, pcTOL, funcs, charDir):
    """
    #MOC direct wall solution using irrotational, isentropic equations 
    pt1: interior point near wall 
    y_x: scalar position function defining wall
    dydx: scalar wall slope function 
    gasProps: gas properties object
    delta: 0 or 1 - 2D or axisymmetric 
    pcTOL: percent change tolerance between successive iterations 
    funcs: moc operator functions object 
    charDir: ("pos"/"neg") denotes what type of characteristics exists between pt1 and nearby wall
    """
    
    #unpacking input data 
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y
    gam, R_gas, T0 = gasProps.gam, gasProps.R, gasProps.T0 #gas properties 

    #getting speed of sound
    a0 = gasProps.a0 #stagnation speed of sound

    #initial value computations 
    u13, v13 = u1, v1 
    y13 = 1 if y1 == 0 else y1

    def solve_direct_wall(u13, v13, y13, first_iter=None):
        a1 = funcs.a(a0, gam, u1, v1)
        a13 = funcs.a(a0, gam, u13, v13)

        if first_iter: 
            if charDir=="pos": lam13 = funcs.lam_plus(u1, v1, a1)
            elif charDir=="neg": lam13 = funcs.lam_min(u1, v1, a1)
        else: 
            a3 = funcs.a(a0, gam, u3, v3)
            if charDir=="pos": lam13 = 0.5*(funcs.lam_plus(u1, v1, a1) + funcs.lam_plus(u3, v3, a3))
            elif charDir=="neg": lam13 = 0.5*(funcs.lam_min(u1, v1, a1) + funcs.lam_min(u3, v3, a3))

        S13 = funcs.S(delta, a13, v13, y13)
        Q13 = funcs.Q(u13, a13)
        R13 = funcs.R(u13, v13, Q13, lam13)

        #computing the x&y location of 3' for first iteration (2 methods suggested)
        try:
            #method 2: use lam13 to project to wall to get (x,y)_3'
            x3p = float(scipy.optimize.fsolve(lambda x: lam13*(x-x1) + y1 - y_x(x), x1))
            #sol = scipy.optimize.root_scalar(lambda x: lam13*(x-x1) + y1 - y_x(x), method='bisect', bracket=(0, 0.25)) #!Not solving with number as ZH book 
            #x3p = sol.root
            y3p = y_x(x3p)  
        except: 
            #method 1: go straight up from x1 to get x_3' if no downstream solution exists given lam13 and point 1
            x3p = x1
            y3p = y_x(x3p)

        #computing partial derivatives for taylor series expansion using 2nd order centered differences
        pfpx_3p = -dydx(x3p)
        pfpy_3p = 1 #!Hard Coded
        f_w3p = 0 #!Hard Coded 

        #construct system and solve
        coeffMat = np.array([[lam13, -1, 0, 0],[pfpx_3p, pfpy_3p, 0, 0],[-S13, 0, Q13, R13],[0, 0, pfpx_3p, pfpy_3p]])
        RHSvec = np.array([lam13*x1-y1, pfpx_3p*x3p+pfpy_3p*y3p-f_w3p, -S13*x1+Q13*u1+R13*v1, 0])
        return np.linalg.solve(coeffMat, RHSvec) #[x3, y3, u3, v3]

    #first iteration
    [x3, y3, u3, v3] = solve_direct_wall(u13, v13, y13, first_iter=True)

    #Iterate above process until values converge
    pc_it = pcTOL
    while abs(pc_it) >= pcTOL:
        x3_old, y3_old, u3_old, v3_old = x3, y3, u3, v3
        u13, v13, y13 = 0.5*(u1 + u3), 0.5*(v1 + v3), 0.5*(y1 + y3)
        [x3, y3, u3, v3] = solve_direct_wall(u13, v13, y13)

        pcVel, pcPos = get_percent_changes([x3_old, y3_old, u3_old, v3_old],[x3, y3, u3, v3]) #get percent change across iteration
        pc_it = max([pcVel, pcPos])

    return [x3, y3, u3, v3] 


def inverse_wall(pt1, pt2, pt3, gasProps, delta, pcTOL, funcs, charDir):
    """
    MOC inverse wall solution using irrotational isentropic equations
    1-2 is existing characteristic (1 downstream of 2)
    3 is wall point downstream of 1-2
    charDir: direction of characteristic from a to 3
    """
    #unpacking input data
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y
    u2, v2, x2, y2 = pt2.u, pt2.v, pt2.x, pt2.y
    x3, y3, thet3 = pt3.x, pt3.y, pt3.thet

    #gas properties 
    gam = gasProps.gam
    a0 = gasProps.a0 #stagnation speed of sound

    #compute values for C- characteristic (2 -> 1)
    u12, v12, y12 = 0.5*(u1+u2), 0.5*(v1+v2), 0.5*(y1+y2)
    a1 = funcs.a(a0, gam, u1, v1)
    a2 = funcs.a(a0, gam, u2, v2)
    a12 = funcs.a(a0, gam, u12, v12)
    if charDir == "pos":
        lam12 = 0.5*(funcs.lam_min(u1, v1, a1) + funcs.lam_min(u2, v2, a2))
    elif charDir == "neg":
        lam12 = 0.5*(funcs.lam_plus(u1, v1, a1) + funcs.lam_plus(u2, v2, a2))
    S12 = funcs.S(delta, a12, v12, y12)
    Q12 = funcs.Q(u12, a12)
    R12 = funcs.R(u12, v12, Q12, lam12)

    def solve_inverse_wall(ya, ua, va, u3, v3): 
        
        ua3, va3, ya3 = 0.5*(ua+u3), 0.5*(va+v3), 0.5*(ya+y3)
        aa3 = funcs.a(a0, gam, ua3, va3)

        aa = funcs.a(a0, gam, ua, va)
        a3 = funcs.a(a0, gam, u3, v3)
        if charDir == "pos":
            lama3 = 0.5*(funcs.lam_plus(ua, va, aa) + funcs.lam_plus(u3, v3, a3))
        elif charDir == "neg":
            lama3 = 0.5*(funcs.lam_min(ua, va, aa) + funcs.lam_min(u3, v3, a3))

        Sa3 = funcs.S(delta, aa3, va3, ya3)
        Qa3 = funcs.Q(ua3, aa3)
        Ra3 = funcs.R(ua3, va3, Qa3, lama3)

        #calculate partial derivatives at point 3 
        pfpx_3 = -math.tan(thet3)
        pfpy_3 = 1

        coeffMat = np.array([[lam12,-1,0,0,0,0], [lama3,-1,0,0,0,0], [-S12,0,Q12,R12,0,0], [-Sa3,0,Qa3,Ra3,-Qa3,-Ra3], [v2-v1,0,0,-(x2-x1),0,0], [0,0,0,0,pfpx_3, pfpy_3]])
        RHSvec = np.array([lam12*x1-y1, lama3*x3-y3, -S12*x1+Q12*u1+R12*v1, -Sa3*x3, (v2-v1)*x1-(x2-x1)*v1,0])
        return np.linalg.solve(coeffMat, RHSvec) #[xa, ya ua, va, u3, v3]

    #first iteration values
    xa, ya, ua, va = 0.5*(x1+x2), 0.5*(y1+y2), 0.5*(u1+u2), 0.5*(v1+v2)
    u3, v3 = u2, v2

    #iterate until solution converges:
    pc_it = pcTOL
    while abs(pc_it) >= pcTOL:

        x3_old, y3_old, u3_old, v3_old = x3, y3, u3, v3
        [xa, ya, ua, va, u3, v3] = solve_inverse_wall(ya, ua, va, u3, v3)

        pcVel, pcPos = get_percent_changes([x3_old, y3_old, u3_old, v3_old],[x3, y3, u3, v3]) #get percent change across iteration
        pc_it = max([pcVel, pcPos])

    return [xa, ya, ua, va, u3, v3]
