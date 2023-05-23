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
    def __init__(self, geom, tmc_res, gasProps, nPts=None, endPoints=None, upstream_idl=False):
        """
        geom: geometry object
        tmc_res: results object from taylor maccoll function computation
        gasProps: gas properties object
        nPts: number of points spaced along initial data line
        endpoints: (for straight IDL) x coordinates of centerbody and cowl
            points denoting ends of data line 
        
        ultimately returns initial data line as object which can be handed off to moc sequence 
        """
        if upstream_idl: 
            #if using class for region of flow upstream of cowl point and positive charcteristic
            self.generate_char_from_cowl_lip_to_incident_shock(geom, tmc_res, gasProps)
            self.p02_p01_inc_shock = tmc_res.p02_p01 #total pressure change across incident shock
            self.get_properties_on_idl(gasProps)
            return 

        if endPoints is not None: 
            self.generate_2_point_idl(geom, tmc_res, nPts, endPoints)

        else: 
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
        import scipy.interpolate as scp_interp

        a0, gam = gasProps.a0, gasProps.gam
        x_i, y_i = geom.x_cowl_lip, geom.y_cowl(geom.x_cowl_lip)
        thet_i = math.atan(y_i/x_i)
        thet_f = tmc_res.cone_ang
        
        def lambda_char(thet):
            #gives slope of + & - mach lines through point x,y in the flow 
            u,v = tmc_res.f_veloc_uv(thet)
            a =  math.sqrt(a0**2 - 0.5*(gam-1)*(u**2 + v**2))
            dydx = (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
            return dydx

        #solving using Euler integration
        #! get this working with a more efficient solver (Euler is quick and dirty fix)
        x, y = [x_i],[y_i]
        thet = thet_i
        delta_x = 1e-4
        while thet >= thet_f:
            x_i, y_i = x[-1], y[-1]
            dydx = lambda_char(thet)
            x_n = x_i - delta_x
            y_n = y_i + dydx*(-delta_x)
            x.append(x_n), y.append(y_n)
            thet = math.atan(y_n/x_n)
        

        #numerically integrate until mach line meets the centerbody
        #sol = scp_int.solve_ivp(lambda_char, t_span=(thet_i, thet_f), y0=y0, dense_output=True)

        #convert to x and y coordinates 
        #thet_interp = np.linspace(thet_i, thet_f, 5000)
        #y_interp = sol.sol(thet_interp)
        #x = y_interp[0]
        #y = y_interp[1]

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
        self.x = list(fx(new_distances))
        self.y = list(fy(new_distances))

        self.u = []
        self.v = []
        for i,X in enumerate(self.x):
            U,V = tmc_res.f_veloc_uv(math.atan(self.y[i]/X))
            self.u.append(U)
            self.v.append(V)
    
    def generate_char_from_cowl_lip_to_incident_shock(self, geom, tmc_res, gasProps):
        """
        docstring
        """
        a0, gam = gasProps.a0, gasProps.gam
        x_i, y_i = geom.x_cowl_lip, geom.y_cowl(geom.x_cowl_lip)
        thet_i = math.atan(y_i/x_i)
        thet_f = tmc_res.shock_ang
        
        def lambda_char(thet):
            #gives slope of + & - mach lines through point x,y in the flow 
            u,v = tmc_res.f_veloc_uv(thet)
            a =  math.sqrt(a0**2 - 0.5*(gam-1)*(u**2 + v**2))
            dydx = (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
            return dydx

        #solving using Euler integration
        self.x,self.y = [x_i],[y_i]
        thet = thet_i
        delta_x = 1e-2
        while thet <= thet_f:
            x_i, y_i = self.x[-1], self.y[-1]
            dydx = lambda_char(thet)
            x_n = x_i + delta_x
            y_n = y_i + dydx*(delta_x)
            self.x.append(x_n), self.y.append(y_n)
            thet = math.atan(y_n/x_n)
        
        self.u, self.v = [],[]
        for i,X in enumerate(self.x):
            U,V = tmc_res.f_veloc_uv(math.atan(self.y[i]/X))
            self.u.append(U)
            self.v.append(V)

    def get_properties_on_idl(self, gasProps):
        #unpacking
        gam, a0, T0, R = gasProps.gam, gasProps.a0, gasProps.T0, gasProps.R

        self.p0 = gasProps.p0*self.p02_p01_inc_shock
        self.mach = []
        self.T, self.p, self.rho = [],[],[]
        self.T_T0, self.p_p0, self.rho_rho0 = [],[],[]
        for i,_ in enumerate(self.x):
            V = math.sqrt(self.u[i]**2 + self.v[i]**2)
            a = math.sqrt(a0**2 - 0.5*(gam-1)*V**2)
            mach = V/a
            self.mach.append(mach)
            self.T.append(T0/(1+0.5*(gam-1)*(V/a)**2))
            self.T_T0.append((1 + 0.5*(gam-1)*mach**2)**-1)
            self.p.append(self.p0*(T0/self.T[i])**(gam/(gam-1)))
            self.p_p0.append(((1 + 0.5*(gam-1)*mach**2)**(gam/(gam-1)))**-1)
            self.rho.append(self.p[-1]/(R*self.T[-1]))
            self.rho_rho0.append(((1 + 0.5*(gam-1)*mach**2)**(1/(gam-1)))**-1)            
            

class Generate_2D_Initial_Data_Line(Generate_TMC_Initial_Data_Line): 
    """
    generates an initial data line using 2D oblique shock relations. Only 
    use if centerbody is a straight ramp from tip to idl. IDL not valid for 
    non-straight geometry upstream of cowl idl
    """
    def __init__(self, inputs, shockObj, gasProps, nPts=None, endPoints=None, upstream_idl=False): 

        T1 = inputs.T0/(1 + 0.5*(gasProps.gam-1)*inputs.M_inf**2)
        shockObj.T2 = shockObj.T2_T1*T1

        if upstream_idl: 
            #if using class for region of flow upstream of cowl point and positive charcteristic
            self.generate_char_from_cowl_lip_to_incident_shock(inputs.geom, shockObj, gasProps)
            self.p02_p01_inc_shock = shockObj.p02_p01
            self.get_properties_on_idl(gasProps)
            return

        if endPoints is not None: 
            self.generate_straight_cowl_lip_idl(inputs.geom, nPts, shockObj, endPoints, gasProps)
        else: 
            self.generate_initial_char_from_cowl_lip(inputs.geom, shockObj, gasProps, nPts)
        
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
    
    def generate_initial_char_from_cowl_lip(self, geom, shockObj, gasProps, nPts): 
        """
        generates an initial positive characteristic from the cowl lip to the 
        centerbody. Use for non-straight centerbody geometry. Only works if 
        point on centerbody is upstream of any non-straight geometry. 
        """        
        
        x_i, y_i = geom.x_cowl_lip, geom.y_cowl(geom.x_cowl_lip)
        thet_f = shockObj.deflec

        a = math.sqrt(gasProps.gam*gasProps.R*shockObj.T2)
        V = shockObj.M2*a
        u,v = V*math.cos(shockObj.deflec), V*math.sin(shockObj.deflec)
        lam_pos = (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
        x = (lam_pos*x_i - y_i)/(lam_pos - math.tan(thet_f))
        y = math.tan(thet_f)*x

        self.x = list(np.linspace(x_i, x, nPts))
        self.y = list(np.linspace(y_i, y, nPts))
        self.u = list(np.multiply(u, np.ones(len(self.x))))
        self.v = list(np.multiply(v, np.ones(len(self.y)))) 

    def generate_char_from_cowl_lip_to_incident_shock(self, geom, shockObj, gasProps, ):
        """
        creates a data points along a positive characteristic extending 
        from the cowl lip, downstream to the incident shock wave
        """
        x_c, y_c = geom.x_cowl_lip, geom.y_cowl(geom.x_cowl_lip)
        thet_f = shockObj.deflec

        a = math.sqrt(gasProps.gam*gasProps.R*shockObj.T2)
        V = shockObj.M2*a
        u,v = V*math.cos(shockObj.deflec), V*math.sin(shockObj.deflec)
        lam_pos = (u*v + a*math.sqrt(u**2 + v**2 - a**2))/(u**2 - a**2)
        y = (y_c - lam_pos*x_c)/(1 - lam_pos/math.tan(shockObj.beta))
        x = y/math.tan(shockObj.beta)

        self.x = [x_c, x]
        self.y = [y_c, y]
        self.u = [u,u]
        self.v = [v,v]
