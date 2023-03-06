import math
import scipy.optimize
import unit_processes as up
import numpy as np
import oblique_shock as obs

"""
Currently this tests out the wall shock point calculation as described in B.H. Anderson Paper
"""

class point:
    def __init__(self, u=None, v=None, x=None, y=None, T=None):
            if u is not None: self.u = u
            if v is not None: self.v = v
            if x is not None: self.x = x
            if y is not None: self.y = y
            if T is not None: self.T = T



def wall_shock_point(pt_w_u, y_x, dydx, pt1, pt3, pcTOL, delta, gasProps, shockDir):
    """
    gets the intial negative shock segment from wall flow deflection to first pos char 
    pt_w_u: upstream-of-shock properties at wall point 
    y_x: wall position function 
    dy_dx: wall slope function 
    pt1: upstream point of char to be intersected
    pt3: downstream point of char to be intersected 
    pcTOL: moc operator percent-change tolerance 
    delta: (0/1) 2D or axisymmetric 
    gasProps: gas properties object
    shockDir: ("neg"/"pos") direction of shock (use neg for wall above, pos for wall below)  
    """

    #unpacking
    gam, R, T0, a0 = gasProps.gam, gasProps.R, gasProps.T0, gasProps.a0
    M_w_i, T_w_i, u_w_i, v_w_i, x_w_i, y_w_i = pt_w_u.mach, pt_w_u.T, pt_w_u.u, pt_w_u.v, pt_w_u.x, pt_w_u.y
    u1,v1,x1,y1,T1 = pt1.u, pt1.v, pt1.x, pt1.y, pt1.T
    u3,v3,x3,y3,T3 = pt3.u, pt3.v, pt3.x, pt3.y, pt3.T

    #load in operator functions 
    f = up.operator_funcs()
    
    #Get initial shock at wall point 
    thet_i = math.atan(v_w_i/u_w_i) #initial flow angle 
    wallDef = math.atan(dydx(x_w_i)) #wall angle 
    delta = abs(wallDef - thet_i) #change in flow direction due to wall 
    
    shockObj = obs.Oblique_Shock(M_w_i, gam, thet=delta)
    beta_w = shockObj.beta_w


    T0w_Tw = 1 + 0.5*(gam - 1)*shockObj.M2**2 
    T_w = T0/T0w_Tw
    V_w = shockObj.M2*math.sqrt(gam*R*T_w)
    u_w = V_w*math.cos(wallDef)
    v_w = V_w*math.sin(wallDef)
    x_w, y_w = x_w_i, y_w_i
    ptw_dwn = point(u_w, v_w, x_w, y_w, T_w) #wall point downstream of shock 

    a1 = f.a(a0, gam, u1, v1)
    a3 = f.a(a0, gam, u3, v3)
    if shockDir=="neg":
        lam1, lam3 = f.lam_plus(u1, v1, a1), f.lam_plus(u3, v3, a3)
    elif shockDir=="pos":
        lam1, lam3 = f.lam_min(u1, v1, a1), f.lam_min(u3, v3, a3)
    
    lam13 = f.lam(lam1, lam3)

    def solve_shock(thet4, beta4, ret="thet"):
        
        lam_s = math.tan(0.5*(beta_w + beta4)) #shock slope
        if shockDir=="neg":
            lam_s = -1*lam_s

        #find intersection point of shock and segment 1-3 (assuming char is a line with slope lam_13)
        x4 = (lam13*x3 - lam_s*x_w - y3 + y_w)/(lam13 - lam_s)
        y4 = lam_s*(x4 - x_w) + y_w

        #linear interpolate to get velocity and temperature 
        linInt = lambda x, p1, p3: (p3-p1)/(x3-x1)*(x-x1) + p1
        u4_i = linInt(x4, u1, u3)
        v4_i = linInt(x4, v1, v3)
        T4_i = linInt(x4, T1, T3)
        M4_i = math.sqrt((u4_i**2 + v4_i**2)/(gam*R*T4_i))

        #apply oblique shock relations to get downstream conditions
        shockObj = obs.Oblique_Shock(M4_i, gam, thet=thet4)
        beta4 = shockObj.beta_w
        
        M4 = shockObj.M2
        T04_T4 = 1 + 0.5*(gam - 1)*M4**2 
        T4 = T0/T04_T4
        V4 = M4*math.sqrt(gam*R*T4)
        thet4_i = math.atan(v4_i/u4_i) #initial flow angle 

        if shockDir == "neg":thet4 = -1*thet4 

        fa4 = thet4_i + thet4 #downstream flow angle after shock 
        u4 = V4*math.cos(fa4)  
        v4 = V4*math.sin(fa4)

        pt4 = point(u4, v4, x4, y4)

        if shockDir == "neg":
            [x3p, y3p, u3p, v3p] = up.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="pos")#get downstream wall point 
            [x_r, y_r, u_r, v_r] = up.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="neg")#get upstream reference point
            pt_r = point(u_r, v_r, x_r, y_r)
            pt3p = point(u3p, v3p, x3p, y3p)
            pt1 = pt3p
            pt2 = pt_r 

        elif shockDir == "pos":
            [x3p, y3p, u3p, v3p] = up.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="neg")#get downstream wall point 
            [x_r, y_r, u_r, v_r] = up.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="pos")#get upstream reference point
            pt_r = point(u_r, v_r, x_r, y_r)
            pt3p = point(u3p, v3p, x3p, y3p)
            pt1 = pt_r
            pt2 = pt3p
            
        #compute new shock point 
        [x4, y4, u4, v4] = up.interior_point(pt1, pt2, gasProps, 1, 0.00001, f)

        #get new flow deflection angle at 4
        V4_i = np.array([u4_i, v4_i])
        V4 = np.array([u4, v4])

        thet4 = abs(math.acos(np.dot(V4, V4_i)/(math.sqrt(u4_i**2 + v4_i**2)*math.sqrt(u4**2 + v4**2))))

        if ret=="thet": return thet4, beta4
        elif ret=="sol":

            pt4_dwn = point(x=x4, y=y4, u=u4, v=v4, T=T4)
            pt4_ups = point(u=u4_i, v=v4_i, T=T4_i)
            return [pt4_dwn, pt4_ups, thet4, beta4, ptw_dwn]


    shock = obs.Oblique_Shock(M_w_i, gam, thet=delta)

    def err_func(thet_guess): 
        thet,_  = solve_shock(thet_guess, shock.beta_w)
        return thet - thet_guess

    res = scipy.optimize.fsolve(err_func, x0=delta)
    thet4 = res[0]
    beta4 = shock.beta_w
    
    return solve_shock(thet4, beta4, ret="sol") 


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

    pt4_dwn, pt4_ups, thet4, beta4, ptw_dwn = wall_shock_point(pt_w_u, y_cowl, dydx_cowl, pt1, pt3, 0.0001, 1, gas, "neg")

    #Now flip over x-axis to test positive shock direction
    v1 = -1*v1
    pt_w_u.v = pt_w_u.v*-1 
    pt_w_u.y = pt_w_u.y*-1
    pt1.y, pt1.v = pt1.y*-1, pt1.v*-1
    pt3.y, pt3.v = pt3.y*-1, pt3.v*-1

    y_cowl_inv = lambda x: y_cowl(x)*-1
    dydx_cowl_inv = lambda x: dydx_cowl(x)*-1
    pt4_dwn_inv, pt4_ups_inv, thet4_inv, beta4_inv, ptw_dwn_inv = wall_shock_point(pt_w_u, y_cowl_inv, dydx_cowl_inv, pt1, pt3, 0.0001, 1, gas, "pos") 