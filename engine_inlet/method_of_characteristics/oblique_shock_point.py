import math
import scipy.optimize
import unit_processes as up
import numpy as np

"""
Currently this tests out the wall shock point calculation as described in B.H. Anderson Paper
"""
#Operator Functions
a = lambda a0, gam, u, v : math.sqrt(a0**2 - 0.5*(gam-1)*(u**2 + v**2))
S = lambda delta, a, v, y : delta*(a**2*v/y)
Q = lambda u, a : u**2 - a**2
R = lambda u, v, Q, lam: 2*u*v - Q*lam
lam = lambda lam1, lam2 : 0.5*(lam1 + lam2)
lam_min = lambda u, v, a : (u*v - a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
lam_plus = lambda u, v, a : (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)



class point:
    def __init__(self, u=None, v=None, x=None, y=None, T=None):
            if u is not None: self.u = u
            if v is not None: self.v = v
            if x is not None: self.x = x
            if y is not None: self.y = y
            if T is not None: self.T = T


def get_oblique_shock_angle(M, thet, gam): 
    """
    calculates the shock wave angle required for a given upstream mach number, flow deflection from upstream angle
    Inputs
        M: upstream mach number 
        thet: (rad) flow deflection (positive number) from upstream direction
        gam: ratio of specific heats
    Return: 
        beta: (rad) shock wave angle 
    """
    def thetBetaM(beta, thet, M, gam):
        return (2/math.tan(beta))*(M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta)) + 2) - math.tan(thet)

    alpha = math.asin(1/M)
    beta_weak = scipy.optimize.root_scalar(thetBetaM, args=(thet,M,gam), method='bisect', bracket=[alpha, 0.5*(alpha + math.pi/2)])
    beta_strong = scipy.optimize.root_scalar(thetBetaM, args=(thet,M,gam), method='bisect', bracket=[0.5*(alpha + math.pi/2), math.pi/2])

    return beta_weak.root, beta_strong.root



def get_downstream_properties(u1, v1, T1, thet, beta, gam):

    M1 = math.sqrt((u1**2 + v1**2)/(gam*R*T1))
    Mn1 = M1*math.sin(beta)
    Mn2 = math.sqrt((1 + 0.5*(gam-1)*Mn1**2)/(gam*Mn1**2 - 0.5*(gam-1)))
    M2 = Mn2/math.sin(beta - thet)

    u2 = u1*math.cos(thet) - v1*math.sin(thet)
    v2 = u1*math.sin(thet) + v1*math.cos(thet)
    T2 = T1*(1 + 2*gam/(gam+1)*(M1**2 - 1))*(2 + (gam-1)*M1**2)/(M1**2*(gam+1))

    k = math.sqrt((M2**2*(gam*R*T2))/((u1*math.cos(thet) - v1*math.sin(thet))**2 + (u1*math.sin(thet) + v1*math.cos(thet))**2))
    u2 = k*u2
    v2 = k*v2 

    return u2, v2, T2



def get_flow_deflection(M, beta, gam): 
    """
    gives the flow deflection angle given a shock wave angle and upstream mach number
    """
    return math.atan(2*((M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta) + 2)))/math.tan(beta))



def shock_point(pt1, pt2, beta1, pt4_up, a0, delta):
    """
    compute downstream shock point from previous upstream shock point
    pt1: upstream shock point (x, y needed)
    pt2: downstream point along neg char from pt1 (x,y,u,v needed)
    beta1: shock angle at point 1
    pt4_up: upstream properties at new shock point
    a0: stagnation speed of sound
    """
    #unpacking 
    x1,y1 = pt1.x, pt1.y
    x2,y2,u2,v2 = pt2.x, pt2.y, pt2.u, pt2.v
    u4_up, v4_up, T4_up = pt4_up.u, pt4_up.v, pt4_up.T

    def solve_shock_point(beta4, returnSpec):
        """
        solve downstream-of-shock properties at shock point 4 given a assumed wave angle at point 4
        """
        #find location of shock point 
        K = math.tan(0.5*(beta1 + beta4))
        a2 = math.sqrt(a0**2 - 0.5*(gam-1)*(u2**2 + v2**2))
        lam_p = (u2*v2 + a2*math.sqrt(u2**2 + v2**2 - a2**2))/(u2**2 - a2**2) #first pass set lam_p = lam_p_2 instead of average between 2 and 4
        x4 = (K*x1 - y1 - lam_p*x2 + y2)/(K - lam_p)
        y4 = K*(x4 - x1) + y1

        #TODO interpolate along existing characteristic to get u, v, T at x4 and y4

        #get downstream properties 
        M4_up = math.sqrt((u4_up**2 + v4_up**2)/(gam*R*T4_up))
        thet4 = get_flow_deflection(M4_up, beta4,gam)
        u4, v4, T4 = get_downstream_properties(u4_up, v4_up, T4_up, thet4, beta4, gam)

        u24 = 0.5*(u2 + u4)
        v24 = 0.5*(v2 + v4)
        y24 = 0.5*(y2 + y4)


        a24 = math.sqrt(a0**2 - 0.5*(gam-1)*(u24**2 + v24**2))
        lam24 = 0.5*(lam2 + lam4)
        Q24 = u24**2 - a24**2
        R24 = 2*u24*v24 - Q24*lam24
        S24 = delta*(a24**2*v24/y24)
        
        if returnSpec == "zero":
            return Q24*(u4-u2) + R24*(v4-v2) - S24*(x4-x2) #compatibility relation along char 24
        if returnSpec is "sol": 
            return x4, y4, u4, v4, T4

    res = scipy.optimize.root_scalar(solve_shock_point, args=("zero"), method='bisect', bracket=[0,1]) #get beta using root finder
    beta_4 = res.root

    x4, y4, u4, v4, T4 = solve_shock_point(beta_4, "sol") #solve for shock point properties using solved-for wave angle
    return x4, y4, u4, v4, T4



def get_wall_shock_neg(pt_w_u, y_x, dy_dx, pt1, pt3, pcTOL, delta, gasProps):
    """
    gets the intial negative shock segment from wall flow deflection to first pos char 
    pt_w_u: upstream-of-shock properties at wall point 
    y_x: wall position function 
    dy_dx: wall slope function 
    pt1: upstream point of pos char 
    pt2: downstream point of pos char 
    """

    #unpacking
    gam, R, T0, a0 = gasProps.gam, gasProps.R, gasProps.T0, gasProps.a0
    M_w_i, T_w_i, u_w_i, v_w_i, x_w_i, y_w_i = pt_w_u.mach, pt_w_u.T, pt_w_u.u, pt_w_u.v, pt_w_u.x, pt_w_u.y
    u1,v1,x1,y1,T1 = pt1.u, pt1.v, pt1.x, pt1.y, pt1.T
    u3,v3,x3,y3,T3 = pt3.u, pt3.v, pt3.x, pt3.y, pt3.T
    
    #Get initial shock at wall point 
    thet = math.atan(v_w_i/u_w_i) #initial flow angle 
    wallDef = math.atan(dy_dx(x_w_i)) #wall angle 
    def_ = abs(thet - wallDef) #change in flow direction due to wall 
    beta_w, _ = get_oblique_shock_angle(M_w_i, def_, gam) #weak shock wave angle
    u_w, v_w, T_w = get_downstream_properties(u_w_i, v_w_i, T_w_i, def_, beta_w, gam) #downstream flow properties at wall shock point
    x_w, y_w = x_w_i, y_w_i

    a1 = a(a0, gam, u1, v1)
    a3 = a(a0, gam, u3, v3)
    lam1 = lam_plus(u1, v1, a1)
    lam3 = lam_plus(u3, v3, a3)
    lam13 = lam(lam1, lam3)

    def solve_shock(thet4, ret="thet"):
        
        lam_s = -math.tan(beta_w) #shock slope

        #find intersection point of shock and segment 1-3 (assuming char is a line with slope lam_13)
        x4 = (lam13*x3 - lam_s*x_w - y3 + y_w)/(lam13 - lam_s)
        y4 = lam_s*(x4 - x1) + y_w

        #linear interpolate to get velocity and temperature 
        linInt = lambda x, p1, p3: (p3-p1)/(x3-x1)*(x-x1) + p1
        u4_i = linInt(x4, u1, u3)
        v4_i = linInt(x4, v1, v3)
        T4_i = linInt(x4, T1, T3)
        M4_i = math.sqrt((u4_i**2 + v4_i**2)/(gam*R*T4_i))

        #apply oblique shock relations to get downstream conditions
        beta4,_ = get_oblique_shock_angle(M4_i, thet4, gam)
        u4, v4, T4 = get_downstream_properties(u4_i, v4_i, T4_i, thet4, beta4, gam)

        #get downstream wall point 
        pt4 = point(u4, v4, x4, y4)
        funcs = up.operator_funcs()
        [x3p, y3p, u3p, v3p] = up.direct_wall_abv(pt4, y_x, dy_dx, gasProps, delta, pcTOL, funcs, charDir="pos")

        #get upstream reference point
        [x_r, y_r, u_r, v_r] = up.direct_wall_abv(pt4, y_x, dy_dx, gasProps, delta, pcTOL, funcs, charDir="neg")

        #compute new shock point 
        pt_r = point(u_r, v_r, x_r, y_r)
        pt3p = point(u3p, v3p, x3p, y3p)

        [x4, y4, u4, v4] = up.interior_point(pt3p, pt_r, gasProps, 1, 0.0001, funcs)

        #get new flow deflection angle at 4
        thet4_old = thet4
        V4_i = np.array([u4_i, v4_i])
        V4 = np.array([u4, v4])

        thet4 = math.acos(np.dot(V4, V4_i)/(math.sqrt(u4_i**2 + v4_i**2)*math.sqrt(u4**2 + v4**2)))


        if ret=="thet": return thet4
        elif ret=="sol":

            pt4_dwn = point(x=x4, y=y4, u=u4, v=v4, T=T4)
            pt4_ups = point(u=u4_i, v=v4_i, T=T4_i)
            return [pt4_dwn, pt4_ups, thet4]

    thet4 = thet#intial guess
    err = 0.001
    while abs(err) >= 0.001: 
        thet4_old = thet4
        thet4 = solve_shock(thet4_old)     
        err = thet4-thet4_old
        print(f"error: {err} rad")
    
    return solve_shock(thet4, ret="sol")

if __name__ == "__main__":

    class Point: 
        def __init__(self, x=None, y=None, u=None, v=None, T=None):
            if x is not None: self.x = x
            if y is not None: self.y = y 
            if u is not None: self.u = u
            if v is not None: self.v = v
            if T is not None: self.T = T

    class GasProps: 
        def __init__(self, gam, R, T0):
            self.gam, self.R, self.T0 = gam, R, T0
            self.a0 = math.sqrt(gam*R*T0)

    #some test numbers: 
    u1, v1 = 500, 50
    T1 = 250 

    thet = math.radians(10)

    gam,R,T0 = 1.4, 287.05, 288.15
    gas = GasProps(gam, R, T0)

    pt_w_u = Point(u=546.4438657,v=41.15194476,x=2,y=0.996173133,T=214.7937135) #shock wall point 
    pt_w_u.mach = math.sqrt((pt_w_u.u**2 + pt_w_u.v**2)/(gas.gam*gas.R*pt_w_u.T))

    pt1 = Point(u=538.772962,v=57.72306326,x=2.05,y=0.869061914,T=215.6046032)
    pt3 = Point(u=539.5732147,v=55.48155192,x=2.156000458,y=0.935621478,T=215.5299505)

    def y_cowl(x):
        if x >= 2 and x <= 4.1:
            A,B,C,D = 0.014656593603382383, -0.155835602414445, 0.48384724402657875, 0.534568305777872
            return A*x**3 + B*x**2 + C*x + D
        else:
            return None 
    def dydx_cowl(x):
        if x >= 2 and x <= 4.1:
            A,B,C = 0.014656593603382383, -0.155835602414445, 0.48384724402657875
            return 3*A*x**2 + 2*B*x + C
        else:
            return None

    #M1 = math.sqrt(u1**2 + v1**2)/math.sqrt(gam*R*T1)
    #beta_w, beta_s = get_oblique_shock_angle(M1, thet, gam) #get shock angle 
    #u2, v2, T2 = get_downstream_properties(u1, v1, T1, thet, beta_w, gam) #compute downstream properties
    pt4_dwn, pt4_ups, thet4 = get_wall_shock_neg(pt_w_u, y_cowl, dydx_cowl, pt1, pt3, 0.0001, 1, gas)

    pass