import math 
import numpy as np 

"""
Development of unit processes for rotational, isentropic, steady, 2D/axisymmetric method of characteristics
From Zucrow and Hoffman Ch 17
"""
class Point:

    def __init__(self,x,y,u,v):
          self.x,self.y,self.u,self.v = x,y,u,v


class Funcs: 

    def __init__(self):
        self.V =        lambda u,v: math.sqrt(u**2 + v**2)
        self.thet =     lambda u,v: math.atan(v/u)
        self.t =        lambda T, V, cp: T - V**2/(2*cp)
        self.p =        lambda P, t, T, gam: P*(t/T)**(gam/(gam-1))
        self.rho =      lambda p, R, t: p/(R*t)
        self.a =        lambda gam, p, rho: math.sqrt(gam*p/rho)
        self.M =        lambda V, a: V/a
        self.alph =     lambda M: math.asin(1/M)
        self.Q =        lambda M, rho, V: math.sqrt(M**2 - 1)/(rho*V**2)
        self.S_plus =   lambda delta, thet, y, M, alph: delta*math.sin(thet)/(y*M*math.cos(thet + alph))
        self.S_min =    lambda delta, thet, y, M, alph: delta*math.sin(thet)/(y*M*math.cos(thet - alph))
        self.lam_plus = lambda thet, alph: math.tan(thet+alph)
        self.lam_min =  lambda thet, alph: math.tan(thet-alph)
        self.T_plus =   lambda S_p, x4, x2, Q_p, p2, thet2: -S_p*(x4-x2) + Q_p*p2 + thet2
        self.T_min =    lambda S_m, x4, x1, Q_m, p1, thet1: -S_m*(x4-x1) + Q_m*p1 - thet1
        self.lam_0 =    lambda thet1, thet2: 0.5*(thet1 + thet2)
        self.lam_12 =   lambda y1, x1, y2, x2: (y1 - y2)/(x1 - x2)
        self.R0 =       lambda rho, V: rho*V
        self.A0 =       lambda a: a**2
        self.T01 =      lambda R0, V3, p3: R0*V3 + p3
        self.T02 =      lambda A0, rho3, p3: p3 - A0*rho3
        self.linInt =   lambda x1,x2,x3,z1,z2: (z2-z1)/(x2-x1)*(x3-x1)+z1

        def solve_point_4(lam_p, x2, y2, lam_m, x1, y1):
            return np.linalg.solve(np.array([[-lam_p, 1],[-lam_m, 1]]), np.array([y2-lam_p*x2, y1-lam_m*x1])) 
        self.linSolvePt4 = solve_point_4

        def solve_point_3(lam_0, x4, y4, lam_12, x2, y2):
            return np.linalg.solve(np.array([[-lam_0, 1],[-lam_12, 1]]), np.array([y4-lam_0*x4, y2-lam_12*x2])) 
        self.linSolvePt3 = solve_point_3

        def solve_point_4_props(Q_p, Q_m, T_p, T_m, R0, A0, T01, T02):
            a = np.array([[Q_p,1,0,0],[Q_m,-1,0,0],[1,0,R0,0],[1,0,0,-A0]])
            b = np.array([T_p, T_m, T01, T02])
            return np.linalg.solve(a,b)
        self.linSolvePt4Props = solve_point_4_props
             

def solve_interior(pt1, pt2, gasProps, delta, pcTOL, funcs): 

    #unpack
    u1, v1, x1, y1 = pt1.u, pt1.v, pt1.x, pt1.y
    u2, v2, x2, y2 = pt2.u, pt2.v, pt2.x, pt2.y
    R, T, cp, P, gam = gasProps.R, gasProps.T, gasProps.cp, gasProps.P, gasProps.gam

    #step 0 - get initial values at point 1 and 2
    V1, V2 = funcs.V(u1,v1), funcs.V(u2,v2)
    t1, t2 = funcs.t(T, V1, cp), funcs.t(T, V2, cp)
    thet1, thet2 = funcs.thet(u1, v1), funcs.thet(u2, v2)
    p1, p2 = funcs.p(P, t1, T, gam), funcs.p(P, t2, T, gam)
    rho1, rho2 = funcs.rho(p1, R, t1), funcs.rho(p2, R, t2)

    #step 1 - solve for a, M, alpha at points 1 and 2 
    a1, a2 = funcs.a(gam, p1, rho1), funcs.a(gam, p2, rho2)
    M1, M2 = funcs.M(V1, a1), funcs.M(V2, a2)
    alph1, alph2 = funcs.alph(M1), funcs.alph(M2)

    #step 2 - calculate coefficients: Q, S, lambda for pos and neg characteristics 
    Q_p, Q_m = funcs.Q(M2, rho2, V2), funcs.Q(M1, rho1, V1)
    S_p, S_m = funcs.S_plus(delta, thet2, y2, M2, alph2), funcs.S_min(delta, thet1, y1, M1, alph1)
    lam_p, lam_m = funcs.lam_plus(thet2, alph2), funcs.lam_min(thet1, alph1)

    #step 3 - solve linear system to obtain position of point 4 
    [x4, y4] = funcs.linSolvePt4(lam_p, x2, y2, lam_m, x1, y1)
    T_p, T_m = funcs.T_plus(S_p, x4, x2, Q_p, p2), funcs.T_min(S_m, x4, x1, Q_m, p1, thet1) 

    #step 4 - iteratively locate point 3  
    thet3 = 0.5*(thet1 + thet2) #intialization value 
    lam12 = (y1 - y2)/(x1 - x2)
    pcChange = pcTOL 
    while pcChange >= pcTOL: 
        lam_0 = math.tan(thet3)
        x3_old, y3_old = x3, y3 
        [x3, y3] = funcs.linSolvePt3(lam_0, x4, y4, lam12, x2, y2)
        posChange = [x3-x3_old, y3-y3_old] 
        pcChange = math.sqrt(posChange[0]**2 + posChange[1]**2)/math.sqrt(x3_old**2 + y3_old**2)
        thet3 = thet2 + (y3-y2)/(y1-y2)*(thet1 - thet2)

    #linear interpolate to get properties at point 3
    p3 = funcs.linInt(x1,x2,x3,p1,p2)
    rho3 = funcs.linInt(x1,x2,x3,rho1,rho2)
    V3 = funcs.linInt(x1,x2,x3,V1,V2)
    a3 = funcs.a(gam, p3, rho3)

    #step 5 - calculate p, thet, V, rho at point 4 for the predictor
    R0 = funcs.R0(rho3, V3)
    A0 = funcs.A0(a3)
    T01, T02 = funcs.T01(R0, V3, p3), funcs.T02(A0, rho3, p3)

    [p4, thet4, V4, rho4] = funcs.linSolvePt4Props(Q_p, Q_m, T_p, T_m, R0, A0, T01, T02)

    pcChange = pcTOL
    while pcChange >= pcTOL:
        #step 6 - calculate average properties along characteristics 

        
        #step 7 - solve for point 4

        
        #step 8 - solve for point 3 iteratively


        #step 9 - solve for new point 4 


        #check if convergence has been reached

        pass 
