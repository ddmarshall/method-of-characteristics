import os 
import sys 
sys.path.append(os.getcwd() + "\\taylor_maccoll_cone") #add path to taylor maccoll module
import taylor_maccoll as tmc
import math 
import scipy.optimize as scp_opt
import scipy.interpolate as scp_int
import numpy as np 
import matplotlib.pyplot as plt

"""
Playing around with taylor maccoll function and generating an initial data line for moc computations
"""
class generate_tmc_initial_data_line: 
    def __init__(self, tmc_res, curve):
        """
        tmc_res: results object from taylor maccoll function computation
        curve: object containing attributes: 
            curve.y_x - y(x) curve of initial data line
                or 
            curve.coords - form ((x's),(y's)) - discrete form 

            curve.dist - array-like from 0 to 1 denoting distribution of data points along line starting from x_a
            curve.endpoints - (x_a,x_b) endpoints of data line (only needed if curve.y_x is specified)

        ultimately return initial data line as object which can be handed off to moc sequence 
        """
        self.curveParams = curve
        self.check_inputs(tmc_res)
        self.generate_idl(tmc_res)

    def check_inputs(self, tmc_res):
        """
        TODO: docstring
        TODO: check for geometry collision
        """
        #check if curve lies between cone surface and shock surface - throw error if so

        return

    def generate_idl(self, tmc_res):
        """
        TODO: docstring
        """
        if hasattr(self.curveParams, 'y_x'):
            y_x = self.curveParams.y_x
        elif hasattr(self.curveParams, 'coords'):
            y_x = scp_int.interpolate.interp1d(self.curveParams.coords[0], self.curveParams.coords[1], kind='linear')
        else: 
            raise ValueError("Missing Curve")

        #compute length of curve from endpoint to endpoint
        n = 1000 #!Placeholder
        x_pts = list(np.linspace(curve.endpoints[0], curve.endpoints[1], n))
        y_pts = [y_x(x) for x in x_pts]
        
        l = 0 
        dist = lambda x1, x2, y1, y2: math.sqrt((x2-x1)**2 + (y2-y1)**2)
        for i in range(1, len(x_pts)): 
            dl = dist(x_pts[i], x_pts[i-1], y_pts[i], y_pts[i-1]) #length of segment
            l += dl
        
        spacing = [l*x for x in curve.dist]

        #iterating through curve
        x_alpha = curve.endpoints[0]
        y_alpha = y_x(x_alpha)

        xlist, ylist, ulist, vlist = [], [], [], []
        for i in range(1,len(spacing)): 
            #calculate and store x,y,u,v at point
            xlist.append(x_alpha), ylist.append(y_alpha)
            u,v = tmc_res.f_veloc_uv(math.atan(y_alpha/x_alpha))
            ulist.append(u), vlist.append(v)
            
            #get next point along curve 
            dl = spacing[i]- spacing[i-1]
            func = lambda x_beta: (x_beta - x_alpha)**2 + (y_x(x_beta) - y_alpha)**2 - dl**2
            x_beta = float(scp_opt.fsolve(func, 1.05*x_alpha))
            
            #update position
            x_alpha = x_beta
            y_alpha = curve.y_x(x_alpha)
        
        #x&y and u&v discrete points on data line
        self.x_idl = xlist
        self.y_idl =  ylist
        self.u_idl = ulist
        self.v_idl = vlist

def plot_tmc_idl(idl,coneSol, xinterval, annotate=False):
    plt.figure(figsize=(16,9)), plt.title(f"M = {coneSol.M_inf}, \u03B3 = {coneSol.gam}, R = {coneSol.R} J/(kg*K), T_0 = {coneSol.T0} K")
    #cone surface
    xint = np.linspace(xinterval[0], xinterval[1], 2)

    plt.plot(xint, [x*math.tan(coneSol.cone_ang) for x in xint], label=f'cone = {round(math.degrees(coneSol.cone_ang),2)} deg', color='k')
    plt.hlines(0, 0, xinterval[1], color='k')
    plt.vlines(xinterval[1],0, math.tan(coneSol.cone_ang)*xinterval[1], color = 'k')
    plt.plot(xint, [x*math.tan(coneSol.shock_ang) for x in xint], label=f'shock = {round(math.degrees(coneSol.shock_ang),2)} deg', color='r')
    plt.plot(idl.x_idl, idl.y_idl, '-o', label="idl")

    if annotate: 
        for i,x in enumerate(idl.x_idl):
            text = f"V={round(idl.u_idl[i],1)}, {round(idl.v_idl[i],1)}"
            xy = (x,idl.y_idl[i])
            plt.annotate(text, xy)

    plt.xlabel('x'), plt.ylabel('y'), plt.legend(), plt.show()

if __name__ == "__main__":
    gam = 1.4
    cone_ang = math.radians(30)
    M_inf = 2
    R = 287.05
    T0 = 288.15
    cone = tmc.TaylorMaccoll_Cone(cone_ang, M_inf, gam, R, T0) 

    class make_curve:
        def __init__(self, y_x, dist, endpoints):
            self.y_x, self.dist, self.endpoints = y_x, dist, endpoints
    dist = [0, 0.1, 0.2, 0.25, 0.30, 0.6, 0.8, 1]
    curve =  make_curve(lambda x: 0.5*(x-5)**2+3, dist, (3.6,5))
    idl = generate_tmc_initial_data_line(cone, curve)
    plot_tmc_idl(idl, cone, (0,6), annotate=True)
    pass