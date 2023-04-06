import math
import scipy.optimize
import unit_processes as up
import unit_processes_ZH as upzh
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


def wall_shock_point(pt_w_u, y_x, dydx, pt1, pcTOL, delta, gasProps, shockDir):
    """
    gets the intial negative shock segment from wall flow deflection to first pos char 
    pt_w_u: wall point where shock will originate
    y_x: wall position function 
    dy_dx: wall slope function 
    pt1: interior point to connect to
    pcTOL: moc operator percent-change tolerance 
    delta: (0/1) 2D or axisymmetric 
    gasProps: gas properties object
    shockDir: ("neg"/"pos") direction of shock (use neg for wall above, pos for wall below)  
    """
    #unpacking
    gam, R, T0, a0 = gasProps.gam, gasProps.R, gasProps.T0, gasProps.a0
    M_w_i, T_w_i, u_w_i, v_w_i, x_w_i, y_w_i = pt_w_u.mach, pt_w_u.T, pt_w_u.u, pt_w_u.v, pt_w_u.x, pt_w_u.y
    u1,v1,x1,y1,T1 = pt1.u, pt1.v, pt1.x, pt1.y, pt1.T

    #load in operator functions 
    f = up.operator_funcs()

    #generate pt 3
    if shockDir == "neg":
        #[x3, y3, u3, v3] = up.interior_point(pt1, pt_w_u, gasProps, delta, pcTOL, f)
        [x3, y3, u3, v3] = upzh.interior_point(pt1, pt_w_u, gasProps, delta, pcTOL, f)

        #print(f"interior point solution: x,y,u,v = {x3,y3,u3,v3}")
    elif shockDir == "pos":
        #[x3, y3, u3, v3] = up.interior_point(pt_w_u, pt1, gasProps, delta, pcTOL, f)
        [x3, y3, u3, v3] = upzh.interior_point(pt_w_u, pt1, gasProps, delta, pcTOL, f)

        #print(f"interior point solution: x,y,u,v = {x3,y3,u3,v3}")
    a3 = f.a(a0, gam, u3, v3)
    M3 = math.sqrt(u3**2 + v3**2)/a3
    T3 = T0/(1 + 0.5*(gam-1)*M3**2)

    #Get initial shock at wall point 
    thet_i = math.atan(v_w_i/u_w_i) #initial flow angle 
    wallDef = math.atan(dydx(x_w_i)) #wall angle 
    def_ = abs(wallDef - thet_i) #absolute change in flow direction due to wall 
    #print(f"wall flow deflection = {round(math.degrees(def_),3)}")
    
    shockObj = obs.Oblique_Shock(M_w_i, gam, thet=def_)
    beta_wall = shockObj.beta_w #shock wave angle at the wall
    #print(f"wall shock angle = {round(math.degrees(beta_wall),4)}\n\n")

    T0w_Tw = 1 + 0.5*(gam - 1)*shockObj.M2**2 
    T_w = T0/T0w_Tw
    V_w = shockObj.M2*math.sqrt(gam*R*T_w)
    u_w = V_w*math.cos(wallDef)
    v_w = V_w*math.sin(wallDef)
    x_w, y_w = x_w_i, y_w_i
    ptw_dwn = point(u=u_w, v=v_w, x=x_w, y=y_w, T=T_w) #wall point downstream of shock 

    a1 = f.a(a0, gam, u1, v1)
    a3 = f.a(a0, gam, u3, v3)
    if shockDir=="neg":
        lam1, lam3 = f.lam_plus(u1, v1, a1), f.lam_plus(u3, v3, a3)
    elif shockDir=="pos":
        lam1, lam3 = f.lam_min(u1, v1, a1), f.lam_min(u3, v3, a3)
    
    lam13 = f.lam(lam1, lam3) #average slope between point 1 and 3 

    def solve_shock(beta4, ret="posErr"):
        
        if beta4 < 0: 
            raise ValueError("beta (magnitude) must be greater than 0")

        lam_s = math.tan(0.5*(beta_wall + beta4)) #shock slope
        if shockDir == "neg":
            lam_s = -1*lam_s

        #find intersection point of shock and char 1-3 (assuming char is a line with slope lam_13)
        x4 = (lam_s*x_w - y_w - lam13*x1 + y1)/(lam_s - lam13)
        y4 = lam_s*(x4 - x_w) + y_w

        #linear interpolate to get velocity and temperature of upstream side of shock point
        linInt = lambda x, p1, p3: ((p3-p1)/(x3-x1))*(x-x1) + p1
        u4_i = linInt(x4, u1, u3)
        v4_i = linInt(x4, v1, v3)
        T4_i = linInt(x4, T1, T3)
        M4_i = math.sqrt((u4_i**2 + v4_i**2)/(gam*R*T4_i))
        #print(f"initial shock point location: x={round(x4,4)}, y={round(y4,4)}")
        #print(f"initial upstream velocity: u={round(u4_i,4)}, v={round(v4_i,4)}, M = {round(M4_i,3)}")

        #apply oblique shock relations to get downstream conditions at shock point
        shockObj = obs.Oblique_Shock(M4_i, gam, beta=beta4)
        thet4 = shockObj.thet #flow deflection magnitude
        #print(f"initial flow deflection = {round(math.degrees(thet4),3)}")
        
        M4 = shockObj.M2
        T04_T4 = 1 + 0.5*(gam - 1)*M4**2 
        T4 = T0/T04_T4
        V4 = M4*math.sqrt(gam*R*T4)
        thet4_i = math.atan(v4_i/u4_i) #initial flow angle 

        if shockDir == "neg":
            fa4 = thet4_i - thet4 #downstream flow angle after shock 

        elif shockDir == "pos":
            fa4 = thet4_i + thet4

        u4 = V4*math.cos(fa4)  
        v4 = V4*math.sin(fa4)
        #print(f"initial downstream velocity: u={round(u4,4)}, v={round(v4,4)}, M={round(M4,3)}")
    
        pt4 = point(u4, v4, x4, y4)

        if shockDir == "neg":
            #[x3p, y3p, u3p, v3p] = up.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="pos")#get downstream wall point 
            #[x_r, y_r, u_r, v_r] = up.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="neg")#get upstream reference point
            [x3p, y3p, u3p, v3p] = upzh.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="pos")#get downstream wall point 
            [x_r, y_r, u_r, v_r] = upzh.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="neg")#get upstream reference point
            pt_r = point(u_r, v_r, x_r, y_r)
            pt3p = point(u3p, v3p, x3p, y3p)
            pt1 = pt3p
            pt2 = pt_r 

        elif shockDir == "pos":
            #[x3p, y3p, u3p, v3p] = up.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="neg")#get downstream wall point 
            #[x_r, y_r, u_r, v_r] = up.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="pos")#get upstream reference point
            [x3p, y3p, u3p, v3p] = upzh.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="neg")#get downstream wall point 
            [x_r, y_r, u_r, v_r] = upzh.direct_wall(pt4, y_x, dydx, gasProps, delta, pcTOL, f, charDir="pos")#get upstream reference point
            pt_r = point(u_r, v_r, x_r, y_r)
            pt3p = point(u3p, v3p, x3p, y3p)
            pt1 = pt_r
            pt2 = pt3p
            
        #compute new shock point 
        x4_i, y4_i = x4, y4
        #[x4, y4, u4, v4] = up.interior_point(pt1, pt2, gasProps, delta, pcTOL, f)
        [x4, y4, u4, v4] = upzh.interior_point(pt1, pt2, gasProps, delta, pcTOL, f)

        #print(f"updated shock point location: x={round(x4,4)}, y={round(y4,4)}" )
        #print(f"updated downstream velocity: u={round(u4, 4)}, v={round(v4,4)}")

        if ret=="posErr": return math.sqrt((x4-x4_i)**2 + (y4-y4_i)**2)
        elif ret=="sol":

            #get new flow deflection angle at 4
            V4_i = np.array([u4_i, v4_i])
            V4 = np.array([u4, v4])
            thet4_upd = abs(math.acos(np.dot(V4, V4_i)/(math.sqrt(u4_i**2 + v4_i**2)*math.sqrt(u4**2 + v4**2))))
            #print(f"updated flow deflection = {round(math.degrees(thet4_upd),3)}\n")

            pt4_dwn = point(x=x4, y=y4, u=u4, v=v4, T=T4)
            pt4_ups = point(x=x4, y=y4, u=u4_i, v=v4_i, T=T4_i)
            return [pt4_dwn, pt4_ups, thet4_upd, beta4, ptw_dwn]

    def errFunc(beta):
        #print(f"\n**TESTING beta4 = {round(math.degrees(beta),4)}**")
        try: 
            posErr = solve_shock(beta)
            return posErr #deflection angle error
        except: 
            return None

    sol = scipy.optimize.minimize_scalar(errFunc, bounds=[0.8*beta_wall, 1.2*beta_wall])
    beta4 = float(sol.x)
    return solve_shock(beta4, ret="sol")


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
    gam,R,T0 = 1.4, 287.05, 288.15
    gas = GasProps(gam, R, T0)

    pt_w_u = Point(u=546.4438657,v=41.15194476,x=2,y=0.996173133,T=214.7937135) #shock wall point 
    pt_w_u.mach = math.sqrt((pt_w_u.u**2 + pt_w_u.v**2)/(gas.gam*gas.R*pt_w_u.T))
    pt1 = Point(u=538.772962,v=57.72306326,x=2.05,y=0.869061914,T=215.6046032)

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

    pt4_dwn, pt4_ups, thet4, beta4, ptw_dwn = wall_shock_point(pt_w_u, y_cowl, dydx_cowl, pt1, 0.000001, 1, gas, "neg")
    print(f"\nx={pt4_dwn.x}, \ty={pt4_dwn.y}, \tu={pt4_dwn.u}, \tv={pt4_dwn.v}")

    #Now flip over x-axis to test positive shock direction
    pt_w_u.v = pt_w_u.v*-1 
    pt_w_u.y = pt_w_u.y*-1
    pt1.y, pt1.v = pt1.y*-1, pt1.v*-1

    y_cowl_inv = lambda x: y_cowl(x)*-1
    dydx_cowl_inv = lambda x: dydx_cowl(x)*-1
    pt4_dwn_inv, pt4_ups_inv, thet4_inv, beta4_inv, ptw_dwn_inv = wall_shock_point(pt_w_u, y_cowl_inv, dydx_cowl_inv, pt1, 0.000001, 1, gas, "pos") 

    print(f"x={pt4_dwn_inv.x}, \ty={pt4_dwn_inv.y}, \tu={pt4_dwn_inv.u}, \tv={pt4_dwn_inv.v}\n")