import math
import scipy.optimize
#import unit_processes as up #irrotational moc
import method_of_characteristics.unit_processes as up #irrotational moc
import numpy as np
#import oblique_shock as obs
import method_of_characteristics.oblique_shock as obs

"""
Currently, this tests out the wall shock point calculation as described in B.H. Anderson Paper. Makes use of IRROTATIONAL unit processes
"""

class point:
    def __init__(self, u=None, v=None, x=None, y=None, T=None):
            if u is not None: self.u = u
            if v is not None: self.v = v
            if x is not None: self.x = x
            if y is not None: self.y = y
            if T is not None: self.T = T


def wall_shock_point(pt_w_ups, y_x, dydx, pt1, pcTOL, delta, gasProps, shockDir):
    """
    Generates a shock wave point in the interior of the flow from shock point on a solid boundary
    Inputs: 
        pt_w_ups: (point object) upstream shock point on solid boundary
        y_x: (function of x) continuous wall shape 
        dydx: (function of x) continous first derivative of y_x
        pt1: (point object) interior point upstream of shock wave and near boundary 
        pcTOL: percent change tolerance for iteration cutoff 
        delta: (0 or 1) 2D or Axisymmetric characteristic solutions 
        gasProps: (object) fluid properties 
        shockDir: ("neg" or "pos") denoting the direction of the shock wave
    Returns: 
        [pt4_dwn, pt4_ups, def_4, beta4, ptw_dwn, pt3p]
        pt4_dwn: (point object) downstream side of new shock point 
        pt4_ups: (point object) upstream side of new shock point 
        def_4: (rad) flow deflection angle at pt4
        beta4: (rad) shock wave angle at pt4
        ptw_dwn: (point object) downstream side of wall shock point 
        pt3p: (point object) interior point downstream of shock
    """
    #unpacking
    gam, R, T0, a0              = gasProps.gam, gasProps.R, gasProps.T0, gasProps.a0
    u_w_ups, v_w_ups, x_w_ups, y_w_ups  = pt_w_ups.u, pt_w_ups.v, pt_w_ups.x, pt_w_ups.y
    u1,v1,x1,y1                 = pt1.u, pt1.v, pt1.x, pt1.y

    #load in operator functions 
    f = up.operator_funcs()

    #generate pt 3
    if shockDir == "neg":
        [x3, _, u3, v3] = up.interior_point(pt1, pt_w_ups, gasProps, delta, pcTOL, f)
    elif shockDir == "pos":
        [x3, _, u3, v3] = up.interior_point(pt_w_ups, pt1, gasProps, delta, pcTOL, f)
    a3 = f.a(a0, gam, u3, v3)

    #Get initial shock at wall point 
    thet_w_i = math.atan(v_w_ups/u_w_ups) #initial flow angle 
    wallFlowAng = math.atan(dydx(x_w_ups)) #wall angle 
    def_w = wallFlowAng-thet_w_i # change in flow direction due to wall 
    print(f"\n\tflow deflection at the wall: {math.degrees(def_w)} deg")
    
    a_w_ups = f.a(a0, gam, u_w_ups, v_w_ups)
    M_w_ups = math.sqrt(u_w_ups**2 + v_w_ups**2)/a_w_ups
    shockObj = obs.Oblique_Shock(M_w_ups, gam, deflec=def_w)
    beta_wall = shockObj.beta + thet_w_i #shock wave angle at the wall
    print(f"\tshock wave angle at the wall: {math.degrees(beta_wall)} deg\n")

    T0w_Tw = 1 + 0.5*(gam - 1)*shockObj.M2**2 
    T_w = T0/T0w_Tw
    V_w = shockObj.M2*math.sqrt(gam*R*T_w)
    u_w = V_w*math.cos(wallFlowAng)
    v_w = V_w*math.sin(wallFlowAng)
    x_w, y_w = x_w_ups, y_w_ups
    ptw_dwn = point(u=u_w, v=v_w, x=x_w, y=y_w, T=T_w) #wall point downstream of shock 

    a1 = f.a(a0, gam, u1, v1)
    a3 = f.a(a0, gam, u3, v3)
    if shockDir=="neg":
        lam1, lam3 = f.lam_plus(u1, v1, a1), f.lam_plus(u3, v3, a3)
    elif shockDir=="pos":
        lam1, lam3 = f.lam_min(u1, v1, a1), f.lam_min(u3, v3, a3)
    
    lam13 = f.lam(lam1, lam3) #average slope between point 1 and 3 
   
    def solve_shock(def_4, ret="def"):

        #locate point 4
        lam_w = math.tan(beta_wall)
        a = np.array([[-lam_w, 1],[-lam13, 1]])
        b = np.array([y_w - lam_w*x_w, y1 - lam13*x1])
        x4, y4 = np.linalg.solve(a,b)
        #print(f"initial shock point x,y = {round(x4,4), round(y4,4)}")

        #linear interpolate to get velocity components
        u4_ups = linear_interpolate(x4, u1, u3, x1, x3)
        v4_ups = linear_interpolate(x4, v1, v3, x1, x3)
        a4_ups = f.a(a0, gam, u4_ups, v4_ups)
        M4_ups = math.sqrt((u4_ups**2 + v4_ups**2))/a4_ups        

        #oblique shock relations to get downstream conditions at point 4
        shockObj = obs.Oblique_Shock(M4_ups, gam, deflec=def_4)
        M4_dwn = shockObj.M2
        T04_T4 = 1 + 0.5*(gam - 1)*M4_dwn**2 
        T4 = T0/T04_T4
        V4 = M4_dwn*math.sqrt(gam*R*T4)
        thet4_ups = math.atan(v4_ups/u4_ups) #initial flow angle 
        thet4_dwn = thet4_ups + def_4#downstream flow angle after shock
        beta4 = shockObj.beta + thet4_ups
        u4 = V4*math.cos(thet4_dwn)  
        v4 = V4*math.sin(thet4_dwn)
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

        #check if reference point is upstream of shock wave 
        #!Delete this section if warning never comes up
        a_ref = f.a(a0, gam, pt_r.u, pt_r.v)
        M_ref = math.sqrt(pt_r.u**2 + pt_r.v**2)/a_ref
        mu_ref = math.asin(1/M_ref) + math.atan(pt_r.v/pt_r.u)
        if abs(mu_ref) <= abs(beta_wall - def_w): 
            print("WARNING: reference point upstream of shock wave")

        #compute new shock point 
        [x4, y4, u4, v4] = up.interior_point(pt1, pt2, gasProps, delta, pcTOL, f)
        #print(f"updated shock point x,y = {round(x4,4), round(y4,4)}")
        thet4_upd = math.atan(v4/u4)
        def_4_upd = thet4_upd - thet4_ups #get new estimate for the flow deflection
        
        if ret == "def":
            return def_4_upd 
        if ret == "sol":
            pt4_dwn = point(x=x4, y=y4, u=u4, v=v4)
            pt4_ups = point(x=x4, y=y4, u=u4_ups, v=v4_ups)
            return [pt4_dwn, pt4_ups, def_4_upd, beta4, ptw_dwn, pt3p]

    def_4 = def_w #initial guess deflection 
    defPercChange = pcTOL

    while abs(defPercChange) >= pcTOL:
        def_4_old = def_4
        def_4 = solve_shock(def_4)
        defPercChange = (def_4 - def_4_old)/def_4_old
     
    return solve_shock(def_4, ret="sol")


def interior_shock_point(pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a, pcTOL, delta, gasProps, shockDir):
    """
    Generates a shock wave point in the interior of the flow from an upstream inteior shock point
    Inputs: 
        pt_s_ups: (point object) upstream side of upstream interior shock point
        pt_s_dwn: (point object) downstream side of upstream interior shock point
        beta_s: (rad) shock wave angle at upstream interior shock point
        def_s: (rad) flow deflection through upstream interior shock point 
        pt1: (point object) interior point upstream of shock wave
        pt_a: (point object) interior point downstream of shock wave  
        pcTOL: percent change tolerance for iteration cutoff 
        delta: (0 or 1) 2D or Axisymmetric characteristic solutions 
        gasProps: (object) fluid properties 
        shockDir: ("neg" or "pos") denoting the direction of the shock wave
    Returns: 
        [pt4_dwn, pt4_ups, def4, beta4, pt3p]        
        pt4_dwn: (point object) downstream side of new shock point 
        pt4_ups: (point object) upstream side of new shock point 
        def_4: (rad) flow deflection angle at pt4
        beta4: (rad) shock wave angle at point 4
        pt3p: (point object) interior point downstream of shock
    """
    #unpacking
    gam, R, T0, a0      = gasProps.gam, gasProps.R, gasProps.T0, gasProps.a0
    x_s, y_s            = pt_s_ups.x, pt_s_ups.y
    u_s_dwn, v_s_dwn    = pt_s_dwn.u, pt_s_dwn.v
    u1,v1,x1,y1         = pt1.u, pt1.v, pt1.x, pt1.y
    x_a, u_a, v_a       = pt_a.x, pt_a.u, pt_a.v

    #load in operator functions 
    f = up.operator_funcs()
    a1 = f.a(a0, gam, u1, v1)
    
    #interior point solution to get point 3
    if shockDir == "neg":
        [x3, _, u3, v3] = up.interior_point(pt1, pt_s_ups, gasProps, delta, pcTOL, f)
        a3 = f.a(a0, gam, u3, v3)
        lam3, lam1 = f.lam_plus(u3, v3, a3), f.lam_plus(u1, v1, a1)
    elif shockDir == "pos":
        [x3, _, u3, v3] = up.interior_point(pt_s_ups, pt1, gasProps, delta, pcTOL, f)
        a3 = f.a(a0, gam, u3, v3)
        lam3, lam1 = f.lam_min(u3, v3, a3), f.lam_min(u1, v1, a1)
    
    lam13 = f.lam(lam1, lam3)

    def solve_shock(def_4, ret="def"):

        #locate point 4 using characteristic between 1 and 3. Gives upstream properties of point 4
        lam4s = math.tan(beta_s)

        a = np.array([[1, -lam13],[1, -lam4s]])
        b = np.array([y1 - lam13*x1, y_s - x_s*lam4s])
        y4,x4 = np.linalg.solve(a,b)
        u4_ups, v4_ups = linear_interpolate(x4, u1, u3, x1, x3), linear_interpolate(x4, v1, v3, x1, x3)
        thet4_ups = math.atan(v4_ups/u4_ups)
        a4_ups = f.a(a0, gam, u4_ups, v4_ups)
        pt4_ups = point(u=u4_ups, v=v4_ups, x=x4, y=y4)

        #get downstream properties at shock point 4
        M4_ups = math.sqrt(u4_ups**2 + v4_ups**2)/a4_ups
        shockObj = obs.Oblique_Shock(M4_ups, gam, deflec=def_4)
        M4_dwn = shockObj.M2
        beta4 = shockObj.beta + thet4_ups
        T04_T4 = 1 + 0.5*(gam - 1)*M4_dwn**2 #isentropic stagnation temperature ratio 
        T4_dwn = T0/T04_T4
        a4_dwn = math.sqrt(gam*R*T4_dwn) 
        V4_dwn = M4_dwn*a4_dwn 
        thet4_dws = thet4_ups + def_4

        u4_dwn, v4_dwn = V4_dwn*math.cos(thet4_dws), V4_dwn*math.sin(thet4_dws)
        pt4_dwn = point(x=x4, y=y4, u=u4_dwn, v=v4_dwn, T=T4_dwn)

        #interior point solution to get 3'
        if shockDir=="neg": 
            [x3p, y3p, u3p, v3p] = up.interior_point(pt4_dwn, pt_a, gasProps, delta, pcTOL, f)
        elif shockDir=="pos":
            [x3p, y3p, u3p, v3p] = up.interior_point(pt_a, pt4_dwn, gasProps, delta, pcTOL, f)

        pt3p = point(x=x3p, y=y3p, u=u3p, v=v3p)

        #get reference point
        a_a = f.a(a0, gam, u_a, v_a)
        a_s_dwn = f.a(a0, gam, u_s_dwn, v_s_dwn)
        if shockDir=="neg":
            lam4 = f.lam_min(u4_dwn, v4_dwn, a4_dwn)
            lam_as = f.lam(f.lam_plus(u_a, v_a, a_a), f.lam_plus(u_s_dwn, v_s_dwn, a_s_dwn))
        elif shockDir=="pos":
            lam4 = f.lam_plus(u4_dwn, v4_dwn, a4_dwn)
            lam_as = f.lam(f.lam_min(u_a, v_a, a_a), f.lam_min(u_s_dwn, v_s_dwn, a_s_dwn))

        a = np.array([[1, -lam4],[1, -lam_as]])
        b = np.array([y4 - lam4*x4, y_s - lam_as*x_s])
        y_ref, x_ref = np.linalg.solve(a,b)
        u_ref, v_ref = linear_interpolate(x_ref, u_s_dwn, u_a, x_s, x_a), linear_interpolate(x_ref, v_s_dwn, v_a, x_s, x_a)
        pt_ref = point(u=u_ref, v=v_ref, x=x_ref, y=y_ref)

        #interior point solution to get updated shock point 
        if shockDir == "neg":
            [x4_dwn, y4_dwn, u4_dwn, v4_dwn] = up.interior_point(pt3p, pt_ref, gasProps, delta, pcTOL, f)
        elif shockDir == "pos":
            [x4_dwn, y4_dwn, u4_dwn, v4_dwn] = up.interior_point(pt_ref, pt3p, gasProps, delta, pcTOL, f)

        thet4_upd = math.atan(v4_dwn/u4_dwn)
        def4_upd = thet4_upd - thet4_ups #get new estimate for the flow deflection
        
        if ret == "def":
            return def4_upd 
        if ret == "sol":
            pt4_dwn = point(x=x4_dwn, y=y4_dwn, u=u4_dwn, v=v4_dwn)
            pt4_ups = point(x=x4_dwn, y=y4_dwn, u=u4_ups, v=v4_ups)
            return [pt4_dwn, pt4_ups, def4_upd, beta4, pt3p]
    
    
    #iterative solution approach (seems to diverge occasionally)
    def_4 = def_s #initial guess deflection 
    defPercChange = pcTOL
    print("\nInterior Shock Point Solution:")
    while abs(defPercChange) >= pcTOL:
        def_4_old = def_4
        def_4 = solve_shock(def_4)
        defPercChange = (def_4 - def_4_old)/def_4_old

    print(f"\tconverged deflection: {math.degrees(def_4)} deg")
    return solve_shock(def_4, ret="sol")


def to_wall_shock_point(pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a, y_x, dydx, pcTOL, delta, gasProps, shockDir):
    """
    Generates a shock wave point on a solid boundary from an upstream interior point
    Inputs: 
        pt_s_ups: (point object) upstream side of upstream interior shock point
        pt_s_dwn: (point object) downstream side of upstream interior shock point
        beta_s: (rad) shock wave angle at upstream interior shock point
        def_s: (rad) flow deflection through upstream interior shock point 
        pt1: (point object) interior point upstream of shock wave
        pt_a: (point object) interior point downstream of shock wave
        y_x: (function of x) continous shape of solid boundary 
        dydx: (function of x) first derivative of y_x  
        pcTOL: percent change tolerance for iteration cutoff 
        delta: (0 or 1) 2D or Axisymmetric characteristic solutions 
        gasProps: (object) fluid properties 
        shockDir: ("neg" or "pos") denoting the direction of the shock wave
    Returns: 
        [pt4_dwn, pt4_ups, def_4, beta4, delta_thet_w, pt3p] 
        pt4_dwn: (point object) downstream side of new shock point 
        pt4_ups: (point object) upstream side of new shock point 
        def_4: (rad) flow deflection angle at pt4
        beta4: (rad) shock wave angle at pt4
        delta_thet_w: (rad) difference between downstream flow angle and wall flow angle (non zero value indicates a shock reflection)
        pt3p: (point object) interior point downstream of shock
    """
    #unpacking
    gam, R, T0, a0      = gasProps.gam, gasProps.R, gasProps.T0, gasProps.a0
    x_s, y_s            = pt_s_ups.x, pt_s_ups.y
    u_s_dwn, v_s_dwn    = pt_s_dwn.u, pt_s_dwn.v
    u1,v1,x1            = pt1.u, pt1.v, pt1.x
    x_a, u_a, v_a       = pt_a.x, pt_a.u, pt_a.v

    #load in operator functions 
    f = up.operator_funcs()

    if shockDir == "neg":
        [x3, y3, u3, v3] = up.direct_wall(pt_s_ups, y_x, dydx, gasProps, delta, pcTOL, f, "neg")
        
    elif shockDir == "pos":
        [x3, y3, u3, v3] = up.direct_wall(pt_s_ups, y_x, dydx, gasProps, delta, pcTOL, f, "pos")
        
    def solve_shock(def_4, ret="def"):

        #locate point 4 using characteristic between 1 and 3. Gives upstream properties of point 4
        lam_s = math.tan(beta_s)
        func = lambda x: y_x(x) - y_s - lam_s*(x - x_s)
        x4 = scipy.optimize.fsolve(func, x0=x_s)[0]
        y4 = y_x(x4)

        u4_ups, v4_ups = linear_interpolate(x4, u1, u3, x1, x3), linear_interpolate(x4, v1, v3, x1, x3)
        thet4_ups = dydx(x4)
        a4_ups = f.a(a0, gam, u4_ups, v4_ups)
        pt4_ups = point(u=u4_ups, v=v4_ups, x=x4, y=y4)

        #get downstream properties at shock point 4
        M4_ups = math.sqrt(u4_ups**2 + v4_ups**2)/a4_ups
        shockObj = obs.Oblique_Shock(M4_ups, gam, deflec=def_4)
        M4_dwn = shockObj.M2
        beta4 = shockObj.beta + thet4_ups
        T04_T4 = 1 + 0.5*(gam - 1)*M4_dwn**2 #isentropic stagnation temperature ratio 
        T4_dwn = T0/T04_T4
        a4_dwn = math.sqrt(gam*R*T4_dwn) 
        V4_dwn = M4_dwn*a4_dwn 
        thet4_dws = thet4_ups + def_4

        u4_dwn, v4_dwn = V4_dwn*math.cos(thet4_dws), V4_dwn*math.sin(thet4_dws)
        pt4_dwn = point(x=x4, y=y4, u=u4_dwn, v=v4_dwn, T=T4_dwn)

        #interior point solution to get 3'
        if shockDir=="neg": 
            [x3p, y3p, u3p, v3p] = up.interior_point(pt4_dwn, pt_a, gasProps, delta, pcTOL, f)
        elif shockDir=="pos":
            [x3p, y3p, u3p, v3p] = up.interior_point(pt_a, pt4_dwn, gasProps, delta, pcTOL, f)

        pt3p = point(x=x3p, y=y3p, u=u3p, v=v3p)

        #get reference point
        a_a = f.a(a0, gam, u_a, v_a)
        a_s_dwn = f.a(a0, gam, u_s_dwn, v_s_dwn)
        if shockDir=="neg":
            lam4 = f.lam_min(u4_dwn, v4_dwn, a4_dwn)
            lam_as = f.lam(f.lam_plus(u_a, v_a, a_a), f.lam_plus(u_s_dwn, v_s_dwn, a_s_dwn))
        elif shockDir=="pos":
            lam4 = f.lam_plus(u4_dwn, v4_dwn, a4_dwn)
            lam_as = f.lam(f.lam_min(u_a, v_a, a_a), f.lam_min(u_s_dwn, v_s_dwn, a_s_dwn))

        a = np.array([[1, -lam4],[1, -lam_as]])
        b = np.array([y4 - lam4*x4, y_s - lam_as*x_s])
        y_ref, x_ref = np.linalg.solve(a,b)
        u_ref, v_ref = linear_interpolate(x_ref, u_s_dwn, u_a, x_s, x_a), linear_interpolate(x_ref, v_s_dwn, v_a, x_s, x_a)
        pt_ref = point(u=u_ref, v=v_ref, x=x_ref, y=y_ref)

        #interior point solution to get updated shock point 
        if shockDir == "neg":
            [x4_dwn, y4_dwn, u4_dwn, v4_dwn] = up.interior_point(pt3p, pt_ref, gasProps, delta, pcTOL, f)
        elif shockDir == "pos":
            [x4_dwn, y4_dwn, u4_dwn, v4_dwn] = up.interior_point(pt_ref, pt3p, gasProps, delta, pcTOL, f)
        
        thet4_upd = math.atan(v4_dwn/u4_dwn)
        def4_upd = thet4_upd - thet4_ups #get new estimate for the flow deflection
        
        if ret == "def":
            return def4_upd 
        if ret == "sol":
            pt4_dwn = point(x=x4_dwn, y=y4_dwn, u=u4_dwn, v=v4_dwn)
            pt4_ups = point(x=x4_dwn, y=y4_dwn, u=u4_ups, v=v4_ups)

            #evalate if a reflection is expected to occur
            delta_thet_w = thet4_upd - math.atan(dydx(x4_dwn))
            reflec = False
            if shockDir == "pos" and delta_thet_w > 0: reflec = True
            elif shockDir == "neg" and delta_thet_w < 0: reflec = True 
            elif delta_thet_w == 0:
                reflec = None

            return [pt4_dwn, pt4_ups, def4_upd, beta4, reflec, pt3p]

    def_4 = def_s #initial guess deflection 
    defPercChange = pcTOL

    while abs(defPercChange) >= pcTOL:
        def_4_old = def_4
        def_4 = solve_shock(def_4)
        defPercChange = (def_4 - def_4_old)/def_4_old
        
    return solve_shock(def_4, ret="sol")        


def linear_interpolate(x, z1, z3, x1, x3):
    """
    simple 1-d linear interpolation between x1 and x3 with interpolated value between z1 and z3
    """
    return ((z3-z1)/(x3-x1))*(x-x1) + z1 #linear interpolation function


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
    
    def y_x_up(x):
        if x<1:
            return 1
        else: 
            return 1-math.tan(math.radians(5))*x
    def dydx_up(x):
        if x<1: 
            return 0
        else:
            return -math.tan(math.radians(5))
    def y_x_low(x):
        return 0.5
    def dydx_low(x):
        return 0

    pt_w_ups = Point(u=500, v=0, x=1, y=1)
    pt1 = Point(u=500, v=0, x=1, y=0.75)
    pcTOL = 0.000001 #percent change tolerance for iterative processes 
    delta = 0 #2D geometry 
    [pt4_dwn, pt4_ups, def_4, beta4, ptw_dwn, pt3p] = wall_shock_point(pt_w_ups, y_x_up, dydx_up, pt1, pcTOL, delta, gas, "neg")
    
    V4_ups = math.sqrt(pt4_ups.u**2 + pt4_ups.v**2)
    a4_ups = math.sqrt(gas.a0**2 - 0.5*(gas.gam-1)*V4_ups**2)
    M4_ups = V4_ups/a4_ups
    print("from wall solution:")
    print(f"upstream mach: {M4_ups}")
    print(f"flow deflection: {math.degrees(def_4)}")
    print(f"shock angle: {math.degrees(beta4)} deg")

    V4_dwn = math.sqrt(pt4_dwn.u**2 + pt4_dwn.v**2)
    a4_dwn = math.sqrt(gas.a0**2 - 0.5*(gas.gam-1)*V4_dwn**2)
    M4_dwn = V4_dwn/a4_dwn
    print(f"downstream mach: {M4_dwn}\n")

    print("interior solution:")
    pt1 = Point(u=500, v=0, x=1, y=0.5)
    [pt4_dwn_ref, pt4_ups, def4_upd, beta4, pt3p] = interior_shock_point(pt4_ups, pt4_dwn, beta4, def_4, pt1, pt3p, pcTOL, delta, gas, "neg")
    V4_ups = math.sqrt(pt4_ups.u**2 + pt4_ups.v**2)
    a4_ups = math.sqrt(gas.a0**2 - 0.5*(gas.gam-1)*V4_ups**2)
    M4_ups = V4_ups/a4_ups
    print(f"upstream mach: {M4_ups}")
    print(f"flow deflection: {math.degrees(def_4)}")
    print(f"shock angle: {math.degrees(beta4)} deg")
    V4_dwn = math.sqrt(pt4_dwn.u**2 + pt4_dwn.v**2)
    a4_dwn = math.sqrt(gas.a0**2 - 0.5*(gas.gam-1)*V4_dwn**2)
    M4_dwn = V4_dwn/a4_dwn
    print(f"downstream mach: {M4_dwn}\n")


    print("to wall solution:")
    [pt4_dwn, pt4_ups, def_4, beta4, delta_thet_w, pt3p] = to_wall_shock_point(pt4_ups, pt4_dwn_ref, beta4, def_4, pt1, pt3p, y_x_low, dydx_low, pcTOL, delta, gas, "neg")
    V4_ups = math.sqrt(pt4_ups.u**2 + pt4_ups.v**2)
    a4_ups = math.sqrt(gas.a0**2 - 0.5*(gas.gam-1)*V4_ups**2)
    M4_ups = V4_ups/a4_ups
    print(f"upstream mach: {M4_ups}")
    print(f"flow deflection: {math.degrees(def_4)}")
    print(f"shock angle: {math.degrees(beta4)} deg")
    V4_dwn = math.sqrt(pt4_dwn.u**2 + pt4_dwn.v**2)
    a4_dwn = math.sqrt(gas.a0**2 - 0.5*(gas.gam-1)*V4_dwn**2)
    M4_dwn = V4_dwn/a4_dwn
    print(f"downstream mach: {M4_dwn}\n")

    print("**Expected answer for cases: M4_dwn = 1.77255166\n")

    print("reflection from wall")
    [pt4_dwn, pt4_ups, def_4, beta4, ptw_dwn, pt3p] = wall_shock_point(pt4_dwn, y_x_low, dydx_low, pt4_dwn_ref, pcTOL, delta, gas, "pos")
    print(f"shock angle: {math.degrees(beta4)} deg")