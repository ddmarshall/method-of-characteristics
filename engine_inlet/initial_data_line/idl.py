import math 
import scipy.optimize as scp_opt
import scipy.interpolate as scp_int
import numpy as np 
"""
class for generating initial data line object using flow solution from taylor maccoll module
"""
class generate_tmc_initial_data_line: 
    def __init__(self, tmc_res, curve, gasProps):
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
        self.get_properties_on_idl(gasProps)

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
        x_pts = list(np.linspace(self.curveParams.endpoints[0], self.curveParams.endpoints[1], n))
        y_pts = [y_x(x) for x in x_pts]
        
        l = 0 
        dist = lambda x1, x2, y1, y2: math.sqrt((x2-x1)**2 + (y2-y1)**2)
        for i in range(1, len(x_pts)): 
            dl = dist(x_pts[i], x_pts[i-1], y_pts[i], y_pts[i-1]) #length of segment
            l += dl
        
        spacing = [l*x for x in self.curveParams.dist]

        #iterating through curve
        x_alpha = self.curveParams.endpoints[0]
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
            y_alpha = self.curveParams.y_x(x_alpha)
        
        #x&y and u&v discrete points on data line
        self.x, self.y, self.u, self.v = xlist, ylist, ulist, vlist

    def get_properties_on_idl(self, gasProps):
        #unpacking
        gam, a0, T0, p0 = gasProps.gam, gasProps.a0, gasProps.T0, gasProps.p0
        
        self.T, self.p, self.mach = [],[],[]
        for i,_ in enumerate(self.x):
            V = math.sqrt(self.u[i]**2 + self.v[i]**2)
            a = math.sqrt(a0**2 + 0.5*(gam-1)*V**2)
            self.mach.append(V/a)
            self.T.append(T0/(1+0.5*(gam-1)*(V/a)**2))
            self.p.append(p0*(T0/self.T[i])**(gam/(gam-1)))