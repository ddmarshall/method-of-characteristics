import math
import scipy.optimize as sci_opt
"""
Python implementation of numerical integration on governing equations of Taylor-Maccoll (conical) flow 
Based on Zucrow & Hoffman - Gas Dynamics Vol II - pg.173
"""
def get_Cone_Shock(M, gam, N, delta_c):
    """
    Estimates attached oblique shock and flow properties of axisymmetric flow over an infinite, straight cone
    M: freestream mach number 
    gam: specific heat ratio 
    N: number of angular incremetns between shock wave and cone surface
    delta_c: cone semi-angle (rads)

    returns [eps, M2, ] 
        eps: shock wave angle (rads)
        M2: Mach number immediately behind shock 
        ...
    """
    #calculate trial value for shock wave angle (epsilon)
    eps1 = delta_c + 0.5*math.asin(1/M)

    #numerically determine shock wave angle which corresponds to cone geometry
    def Taylor_Maccoll_Flow(eps, M1, gam, N):
        
        #determine angular step size for numerical integration algorithm 
        del_psi = -(eps - delta_c)/N

        #determine flow properties behind shock wave 
        M1str = math.sqrt(((gam+1)*M1**2)/(2 + (gam-1)*M1**2))
        beta = math.atan(2*math.tan(eps)/(gam+1)*(1/(M1**2*math.sin(eps)**2) + (gam-1)/2))
        thet_s = eps1 - beta

        M2str = M1str*math.sin(eps)/math.sin(beta)*(2/((gam+1)*M1**2*math.sin(eps)**2) + (gam-1)/(gam+1))
        u_strbr = M2str*math.cos(beta)
        v_strbr = -M2str*math.sin(beta)

        #Numerically integrate to the surface of the cone

        def TMc_RK4(gam, u_strbr, v_strbr, psi, del_psi):
            #single iteration of Runge-Kutta 4th order numerical integration for Taylor-Maccoll flow

            a_astr2_1 = (gam+1)/2 - (gam-1)/2*(u_strbr**2 + v_strbr**2)
            l1 = (-u_strbr + (a_astr2_1*(u_strbr + v_strbr/math.tan(psi)))/(v_strbr**2 - a_astr2_1))*del_psi
            k1 = v_strbr*del_psi
            
            A2 = u_strbr + 0.5*k1
            B2 = v_strbr + 0.5*l1
            C2 = psi + 0.5*(del_psi)

            a_astr2_2 = (gam+1)/2 - ((gam-1)/2)*(A2**2 + B2**2)
            l2 = (-A2 + (a_astr2_2*(A2 + B2/math.tan(C2)))/(B2**2 - a_astr2_2))*del_psi
            k2 = B2*del_psi

            A3 = u_strbr + 0.5*k2
            B3 = v_strbr + 0.5*l2

            a_astr2_3 = (gam+1)/2 - ((gam-1)/2)*(A3**2 + B3**2)
            l3 = (-A3 + (a_astr2_3*(A3 + B3/math.tan(C2)))/(B3**2 - a_astr2_3))*del_psi
            k3 = B3*del_psi

            A4 = u_strbr + k3 
            B4 = v_strbr + l3
            psi_upd = psi + del_psi

            a_astr2_4 = (gam+1)/2 - ((gam-1)/2)*(A4**2 + B4**2)
            l4 = (-A4 + (a_astr2_4*(A4 + B4/math.tan(psi_upd)))/(B4**2 - a_astr2_4))*del_psi
            k4 = B4*del_psi

            u_strbr_upd = u_strbr + (k1 + 2*k2 + 2*k3 + k4)/6
            v_strbr_upd = v_strbr + (l1 + 2*l2 + 2*l3 + l4)/6

            return [u_strbr_upd, v_strbr_upd, psi_upd]
        
        [u_strbr_upd, v_strbr_upd, psi_upd] = TMc_RK4(gam, u_strbr, v_strbr, eps, del_psi)

        while v_strbr_upd*v_strbr > 0:
            #break when value of v_strbr switches signs from negative to positive
            u_strbr = u_strbr_upd
            v_strbr = v_strbr_upd
            psi = psi_upd
            [u_strbr_upd, v_strbr_upd, psi_upd] = TMc_RK4(gam, u_strbr, v_strbr, psi, del_psi)
            #print(f" psi = {math.degrees(psi_upd)} (deg) \tu* = {u_strbr_upd} \tv* = {v_strbr_upd}")

        def func(del_psi, gam, u_strbr, v_strbr, psi):
            [_,v_strbr_upd,_] = TMc_RK4(gam, u_strbr, v_strbr, psi, del_psi)
            return v_strbr_upd

        del_psi = sci_opt.fsolve(func, psi_upd - psi, args=(gam, u_strbr, v_strbr, psi)) #solving for del_psi which makes v component = 0
        [u_strbr_upd, v_strbr_upd, psi_upd] = TMc_RK4(gam, u_strbr, v_strbr, psi, del_psi) 

        return [u_strbr_upd, v_strbr_upd, psi_upd, M2str]

    def func(eps, M, gam, N, delta_c):
        [_, _, psi_upd,_] = Taylor_Maccoll_Flow(eps, M, gam, N)
        return abs(psi_upd - delta_c)

    res = sci_opt.minimize(func, eps1, args=(M1, gam, N, delta_c))
    eps = res.x

    #calculate flow using converged-upon shock angle
    [u_brstr, v_brstr, psi, M2_str] = Taylor_Maccoll_Flow(eps, M, gam, N)

    #determine final flow properties on cone surface (using 2d velocity components)
    u_str_c = u_brstr*math.cos(psi) - v_brstr*math.sin(psi)
    v_str_c = u_brstr*math.sin(psi) + v_brstr*math.cos(psi)
    M_str_c = math.sqrt(u_str_c**2 + v_str_c**2)
    thet_c = math.atan(v_str_c/u_str_c)
    
    pcPc = (1 - (gam-1)*M_str_c**2/(gam+1))**(gam/(gam-1))
    rhocrho0c = (1 - (gam-1)*M_str_c**2/(gam+1))**(1/(gam-1))

    #determine final flow property ratios immediately behind shock
    p2p1 = 2*gam/(gam+1)*(M**2*math.sin(eps)**2 - (gam-1)/(2*gam))
    rho1rho2 = 2/(gam+1)*(1/(M**2*math.sin(eps)**2) + (gam-1)/2)

    p2P2 = (1 - (gam-1)*M2_str**2/(gam+1))**(gam/(gam-1))
    rho2rho02 = (1 - (gam-1)*M2_str**2/(gam+1))**(1/(gam-1))

    #determine property ratios between cone surface and immediately behind shock 
    pcp1 = pcPc*p2p1/p2P2
    rhocrho1 = rhocrho0c/(rho2rho02*rho1rho2)

    return [eps, M2_str, pcp1, rhocrho1]
    
#Example 16.9 - pg. 177
delta_c = math.radians(30)
M1 = 3
gam = 1.4
N = 10 

get_Cone_Shock(M1, gam, N, delta_c)

pass