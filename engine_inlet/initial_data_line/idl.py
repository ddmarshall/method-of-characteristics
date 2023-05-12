import math 
import method_of_characteristics.oblique_shock as shock
import numpy as np 
"""
class for generating initial data line object using flow solution from taylor maccoll module
"""
class Generate_TMC_Initial_Data_Line: 
    """
    generates an initial data line using taylor-maccoll equation solution. Only 
    use if centerbody is a straight cone from tip to idl. IDL not valid for 
    non-straight geometry upstream of cowl idl. 
    """
    def __init__(self, geom, tmc_res, gasProps, nPts, endPoints=None):
        """
        geom: geometry object
        tmc_res: results object from taylor maccoll function computation
        gasProps: gas properties object
        nPts: number of points spaced along initial data line
        endpoints: (for straight IDL) x coordinates of centerbody and cowl
            points denoting ends of data line 
        
        ultimately returns initial data line as object which can be handed off to moc sequence 
        """

        if endPoints is not None: 
            self.generate_2_point_idl(geom, tmc_res, nPts, endPoints)

        if endPoints is None: 
            self.generate_initial_char_from_cowl_lip(geom, tmc_res, gasProps, nPts)

        self.p02_p01_inc_shock = tmc_res.p02_p01 #total pressure change across incident shock
        self.get_properties_on_idl(gasProps)

    def generate_2_point_idl(self, geom, tmc_res, nPts, endPoints):
       
        if endPoints[0][1] == "cowl":
            endPoints[0][1] = geom.y_cowl(endPoints[0][0])
            self.cowlPoint = True

        if endPoints[1][1] == "centerbody":
            endPoints[1][1] =  geom.y_centerbody(endPoints[1][0])
            self.cbPoint = True 
            
        a = np.array(endPoints[0])
        b = np.array(endPoints[1])
        spacing = np.linspace(0,1,nPts)
        pts = [np.multiply(spac, b-a) + a for spac in spacing]

        self.x = [pt[0] for pt in pts]
        self.y = [pt[1] for pt in pts]
        
        self.u = []
        self.v = []
        for i,X in enumerate(self.x):
            U,V = tmc_res.f_veloc_uv(math.atan(self.y[i]/X))
            self.u.append(U)
            self.v.append(V)

    def generate_initial_char_from_cowl_lip(self, geom, tmc_res, gasProps, nPts):
        """
        generates an initial positive characteristic from the cowl lip to the 
        centerbody. Use for non-straight centerbody geometry. Only works if 
        point on centerbody is upstream of any non-straight geometry. 
        """
        import scipy.integrate as scp_int
        import scipy.interpolate as scp_interp

        a0, gam = gasProps.a0, gasProps.gam
        x_i, y_i = geom.x_cowl_lip, geom.y_cowl(geom.x_cowl_lip)
        thet_i = math.atan(y_i/x_i)
        thet_f = tmc_res.cone_ang
        y0 = [x_i, y_i]
        
        def lambda_char(thet,y):
            #gives slope of + & - mach lines through point x,y in the flow 
            u,v = tmc_res.f_veloc_uv(thet)
            a =  math.sqrt(a0**2 - 0.5*(gam-1)*(u**2 + v**2))
            dydx = (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
            return [1, dydx]

        #numerically integrate until mach line meets the centerbody
        sol = scp_int.solve_ivp(lambda_char, t_span=(thet_i, thet_f), y0=y0, dense_output=True)

        #convert to x and y coordinates 
        thet_interp = np.linspace(thet_i, thet_f, 5000)
        y_interp = sol.sol(thet_interp)
        x = y_interp[0]
        y = y_interp[1]

        #find total length of curve
        L = 0
        for i,xval in enumerate(x):
            if i == 0: continue
            dl = math.sqrt((xval - x[i-1])**2 + (y[i] - y[i-1])**2)
            L += dl

        #generate evenly-spaced intial data line points
        distances = np.zeros(len(x))
        for i in range(1, len(x)):
            distances[i] = distances[i-1] + np.sqrt((x[i] - x[i-1])**2 + (y[i] - y[i-1])**2)
    
        # Interpolate the curve using distances as the independent variable
        fx = scp_interp.interp1d(distances, x, kind='linear')
        fy = scp_interp.interp1d(distances, y, kind='linear')
    
        # Create an array of evenly spaced distances
        new_distances = np.linspace(0, L, nPts)
    
        # Evaluate the interpolated curve at the new distances
        self.x = fx(new_distances)
        self.y = fy(new_distances)

        self.u = []
        self.v = []
        for i,X in enumerate(self.x):
            U,V = tmc_res.f_veloc_uv(math.atan(self.y[i]/X))
            self.u.append(U)
            self.v.append(V)

    def get_properties_on_idl(self, gasProps):
        #unpacking
        gam, a0, T0 = gasProps.gam, gasProps.a0, gasProps.T0

        self.p0 = gasProps.p0*self.p02_p01_inc_shock
        self.T, self.p, self.mach = [],[],[]
        for i,_ in enumerate(self.x):
            V = math.sqrt(self.u[i]**2 + self.v[i]**2)
            a = math.sqrt(a0**2 - 0.5*(gam-1)*V**2)
            self.mach.append(V/a)
            self.T.append(T0/(1+0.5*(gam-1)*(V/a)**2))
            self.p.append(self.p0*(T0/self.T[i])**(gam/(gam-1)))

class Generate_2D_Initial_Data_Line: 
    """
    generates an initial data line using 2D oblique shock relations. Only 
    use if centerbody is a straight ramp from tip to idl. IDL not valid for 
    non-straight geometry upstream of cowl idl
    """
    def __init__(self, inputs, shockObj, gasProps, nPts, endPoints): 

        T1 = inputs.T0/(1 + 0.5*(gasProps.gam-1)*inputs.M_inf**2)
        shockObj.T2 = shockObj.T2_T1*T1
        self.generate_straight_cowl_lip_idl(inputs.geom, nPts, shockObj, endPoints, gasProps)
        self.p02_p01_inc_shock = shockObj.p02_p01
        self.get_properties_on_idl(gasProps) 

    def generate_straight_cowl_lip_idl(self, geom, nPts, shockObj, endPoints, gasProps): 
        """
        generates straight initial data line which originates at the cowl lip
        with equidistant points along it 
        """
        V2 = shockObj.M2*math.sqrt(gasProps.gam*gasProps.R*shockObj.T2)
        u2, v2 = V2*math.cos(shockObj.deflec), V2*math.sin(shockObj.deflec)
        shockObj.V2 = V2

        if endPoints[0][1] == "cowl":
            endPoints[0][1] = geom.y_cowl(endPoints[0][0])
            self.cowlPoint = True

        if endPoints[1][1] == "centerbody":
            endPoints[1][1] =  geom.y_centerbody(endPoints[1][0])
            self.cbPoint = True 
            
        a = np.array(endPoints[0])
        b = np.array(endPoints[1])
        spacing = np.linspace(0,1,nPts)
        pts = [np.multiply(spac, b-a) + a for spac in spacing]

        self.x = [pt[0] for pt in pts]
        self.y = [pt[1] for pt in pts]
        
        self.u = []
        self.v = []
        for _ in self.x:
            self.u.append(u2)
            self.v.append(v2) 

    def get_properties_on_idl(self, gasProps):
        #unpacking
        gam, a0, T0 = gasProps.gam, gasProps.a0, gasProps.T0

        self.p0 = gasProps.p0*self.p02_p01_inc_shock
        self.T, self.p, self.mach = [],[],[]
        for i,_ in enumerate(self.x):
            V = math.sqrt(self.u[i]**2 + self.v[i]**2)
            a = math.sqrt(a0**2 - 0.5*(gam-1)*V**2)
            self.mach.append(V/a)
            self.T.append(T0/(1+0.5*(gam-1)*(V/a)**2))
            self.p.append(self.p0*(T0/self.T[i])**(gam/(gam-1)))