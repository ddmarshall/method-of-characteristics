from scipy.optimize import curve_fit 
from scipy.interpolate import PchipInterpolator
import math 
"""
Class responsible for parametricising discrete geometry data into polynomial 
equations which can be used with AIMCAT
"""
class Inlet_Geom:

    def __init__(self, inputDict):
        """
        creates AIMCAT geometry object from discrete data
        """
        #unpack input dictionary 
        name, init_turn_ang_deg, geom_type = None, None, None
        cowl_data_list, centerbody_data_list = None, None
        cowl_list, centerbody_list = None, None
        cowl_coord_list, centerbody_coord_list = None, None
        x_lip = 0 
        for key in inputDict.keys():
            if key == "name": name = inputDict[key]
            if key == "initial deflection angle (deg)": init_turn_ang_deg = inputDict[key]
            if key == "2D/Axi": geom_type = inputDict[key]
            if key == "cowl fit data": cowl_data_list = inputDict[key]
            if key == "centerbody fit data": centerbody_data_list = inputDict[key]
            if key == "cowl functions": cowl_list = inputDict[key]
            if key == "centerbody functions": centerbody_list = inputDict[key]
            if key == "cowl coordinates": cowl_coord_list = inputDict[key]
            if key == "centerbody coordinates": centerbody_coord_list = inputDict[key]
            if key == "cowl lip x coord": x_lip = inputDict[key]

        self.x_cowlLip = x_lip
        self.name = name
        self.init_turn_ang_deg = init_turn_ang_deg
        self.geom_type = geom_type

        #starting with centerbody:
        if centerbody_list is not None:
            y_centerbody, dydx_centerbody, centerbody_bounds = self.unpack_surface_function_list(centerbody_list)

        elif centerbody_coord_list is not None:
            y_centerbody, dydx_centerbody, centerbody_bounds = self.process_coords_list(centerbody_coord_list) 

        elif centerbody_data_list is not None: #if geometry includes centerbody data 
            if len(centerbody_data_list) > 1: 
                curve = Composite_Curve(centerbody_data_list)

            elif len(centerbody_data_list) == 1:
                [type_, startpoint, xdata, ydata, startpoint_dydx, x_end, \
                    endpoint] = self.read_sub_dict(centerbody_data_list[0]) 
                
                curve = Curve(type_, startpoint, xdata=xdata, ydata=ydata, \
                              startpoint_dydx=startpoint_dydx, x_end=x_end, \
                                endpoint=endpoint)

            y_centerbody = curve.y
            dydx_centerbody = curve.dydx
            centerbody_bounds = curve.x_bounds 

        if init_turn_ang_deg is None and centerbody_list is not None:
            #if not given, calculate it
            self.init_turn_ang_deg = math.degrees(math.atan(self.dydx_centerbody(curve.x_bounds[0]))) 

        #next cowl
        if cowl_list is not None:
            y_cowl, dydx_cowl, cowl_bounds = self.unpack_surface_function_list(cowl_list)

        elif cowl_coord_list is not None: 
            y_cowl, dydx_cowl, cowl_bounds = self.process_coords_list(cowl_coord_list) 

        elif cowl_data_list is not None: #if geometry includes cowl data
            
            if len(cowl_data_list) > 1: 
                curve = Composite_Curve(cowl_data_list)
            
            elif len(cowl_data_list) == 1: 
                [type_, startpoint, xdata, ydata, startpoint_dydx, x_end, \
                    endpoint] = self.read_sub_dict(cowl_data_list[0])
                
                curve = Curve(type_, startpoint, xdata=xdata, ydata=ydata, \
                              startpoint_dydx=startpoint_dydx, x_end=x_end, \
                                endpoint=endpoint)

            y_cowl = curve.y
            dydx_cowl = curve.dydx    
            cowl_bounds = curve.x_bounds

        #set centerbody functions 
        self.y_centerbody = y_centerbody
        self.dydx_centerbody = dydx_centerbody
        self.centerbody_bounds = centerbody_bounds

        #add cowl lip offset to cowl functions and bounds
        self.y_cowl = lambda x: y_cowl(x-self.x_cowlLip)
        self.dydx_cowl = lambda x: dydx_cowl(x-self.x_cowlLip)
        self.cowl_bounds = [x+self.x_cowlLip for x in cowl_bounds]

    def unpack_surface_function_list(self, surface_list):
        """
        unpacks a subdictionary which defines the surface using functions rather 
        data points. Functions should be in lambda form
        """
        y_funcs, dydx_funcs = [], []
        endpoints = []

        for subDict in surface_list:
            endpoints.append(subDict["x bounds"])

            y_func = eval(subDict["y function"])
            dydx_func = eval(subDict["dydx function"])

            y_funcs.append(y_func), dydx_funcs.append(dydx_func)

        def y_func_combined(x):
            y = None
            for i,func in enumerate(y_funcs):
                if endpoints[i][0] <= x <= endpoints[i][-1]:
                    y = func(x)
                else:continue
            return y

        def dydx_func_combined(x):
            dydx = None
            for i,func in enumerate(dydx_funcs):
                if endpoints[i][0] <= x <= endpoints[i][-1]:
                    dydx = func(x)
                else:continue
            return dydx

        endpoints_flatten = []
        for pt in endpoints: 
            endpoints_flatten += pt

        return y_func_combined, dydx_func_combined, [min(endpoints_flatten), max(endpoints_flatten)]
       
    def process_coords_list(self, surface_list):
        """
        uses pchip interpolation to define the surface
        """
        y_list, dydx_list = [],[]
        x_tot = []
        for curve in surface_list: 
            x_tot += curve["x coords"]
            interpolator = PchipInterpolator(curve["x coords"], curve["y coords"], extrapolate=False)
            def y(x):
                y = None 
                try: return float(interpolator(x)) #position
                except: return y
            def dydx(x):
                dydx = None
                try: return float(interpolator.derivative()(x))
                except: return dydx 

            y_list.append(y)
            dydx_list.append(dydx)

        if len(surface_list) > 1: 
            
            def y_comb(x):
                for func in y_list:
                    y = func(x)
                    if y is not None: 
                        return y
            def dydx_comb(x):
                for func in dydx_list:
                    dydx = func(x)
                    if dydx is not None: 
                        return dydx

            return y_comb, dydx_comb, [min(x_tot), max(x_tot)]

        return y_list[0], dydx_list[0], [min(x_tot), max(x_tot)] 

    def read_sub_dict(self, subDict):
        """
        takes in a curve subdictionary and unpacks it
        """
        keys = subDict.keys()
        type_ = subDict["type"]
        startpoint = subDict["startpoint"]
        startpoint_dydx,x_end,endpoint,xdata,ydata = None,None,None,None,None
        
        if "startpoint dydx" in keys:
            startpoint_dydx = subDict["startpoint dydx"]
        if "endpoint x" in keys:
            x_end = subDict["endpoint x"]
        if "endpoint" in keys: 
            endpoint = subDict["endpoint"]
        if "xdata" in keys: 
            xdata, ydata = subDict["xdata"], subDict["ydata"]    

        return [type_, startpoint, xdata, ydata, startpoint_dydx, x_end, endpoint]


class Curve:
    """
    generates functions defining the position and first derivative of a 
    single curve from input data
    """
    def __init__(self, type_, startpoint, xdata=None, ydata=None, \
                 startpoint_dydx=None, x_end=None, endpoint=None):
        
        if type_=="linear":
            y, dydx = self.linear_poly(startpoint, startPoint_dydx=startpoint_dydx, \
                                  x_end=x_end, endPoint=endpoint)
            if x_end is not None: self.x_bounds = (startpoint[0], x_end)
            elif endpoint is not None: self.x_bounds = (startpoint[0], endpoint[0])

        elif type_=="polynomial":
            y,dydx = self.least_squares_5th_order_poly(xdata, ydata, startpoint,\
                                                startPoint_dydx=startpoint_dydx)
            self.x_bounds =(startpoint[0], max(xdata))

        self.y = y
        self.dydx = dydx

    def least_squares_5th_order_poly(self, xs, ys, startPoint, startPoint_dydx=None):
        """
        returns polymial equation for y_x as well as dydx of a curve defined 
        """
        y_x, dydx_x = None, None
    
        x_a, y_a = startPoint    

        def func_free(x, A, B, C, D, E):
            return A*(x-x_a)**5 + B*(x-x_a)**4 + C*(x-x_a)**3 + D*(x-x_a)**2 + E*(x-x_a) + y_a

        def func_const_deriv(x, A, B, C, D):
            return A*(x-x_a)**5 + B*(x-x_a)**4 + C*(x-x_a)**3 + D*(x-x_a)**2 + startPoint_dydx*(x-x_a) + y_a

        if startPoint_dydx is None: 
            coeffs, _ = curve_fit(func_free, xs, ys)
            A,B,C,D,E = coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4]
            def y_x(x):
                if x_a <= x <= max(xs):
                    return A*(x-x_a)**5 + B*(x-x_a)**4 + C*(x-x_a)**3 + D*(x-x_a)**2 + E*(x-x_a) + y_a
                else: return None
            def dydx_x(x):
                if x_a <= x <= max(xs):
                    return 5*A*(x-x_a)**4 + 4*B*(x-x_a)**3 + 3*C*(x-x_a)**2 + 2*D*(x-x_a) + E
                else: return None 

        if startPoint_dydx is not None: 
            coeffs, _ = curve_fit(func_const_deriv, xs, ys)
            A,B,C,D = coeffs[0],coeffs[1],coeffs[2],coeffs[3]
            def y_x(x):
                if x_a <= x <= max(xs): 
                    return A*(x-x_a)**5 + B*(x-x_a)**4 + C*(x-x_a)**3 + D*(x-x_a)**2 + startPoint_dydx*(x-x_a) + y_a
                else: return None
            def dydx_x(x):
                if x_a <= x <= max(xs):
                    return 5*A*(x-x_a)**4 + 4*B*(x-x_a)**3 + 3*C*(x-x_a)**2 + 2*D*(x-x_a) + startPoint_dydx
                else: return None

        #compile error and print to console
        y_x_err = [abs(y_x(x) - ys[i]) for i,x in enumerate(xs)]
        print(f"\nmaximum least squares position error: {max(y_x_err)}\n\
              average error: {sum(y_x_err)/len(y_x_err)}")

        return y_x, dydx_x

    def linear_poly(self, startPoint, startPoint_dydx=None, x_end=None, endPoint=None):
        """
        returns linear equation for any straight or near straight segments. 
        """
        y_x, dydx_x = None, None 
        x_a, y_a = startPoint

        if startPoint_dydx is not None and x_end is None: 
            raise ValueError("startPoint_dydx and x_end must both be specified")

        if endPoint is not None: #if an endpoint is specified 
            x_b, y_b = endPoint
            m = (y_b-y_a)/(x_b-x_a)

        elif startPoint_dydx is not None: #if start point slope is given
            x_b = x_end
            m = startPoint_dydx

        def y_x(x):
            if x_a <= x <= x_b:
                return m*(x-x_a) + y_a
            else: return None
        def dydx_x(x):
            if x_a <= x <= x_b:
                return m
            else: return None

        return y_x, dydx_x


class Composite_Curve(Inlet_Geom): 
    """
    generates functions defining the position and first derivative of a curve
    made of multiple sub-curves (polynomial or linear) from input data 
    """
    def __init__(self, curves_list):
        
        y_funcs = []
        dydx_funcs = []
        endpoint_prev = None
        endpoint_dydx_prev = None
        x_endpoints = []

        for i,subDict in enumerate(curves_list):
            
            [type_, startpoint, xdata, ydata, startpoint_dydx, x_end, endpoint] = self.read_sub_dict(subDict)
            if i == 0: x_endpoints.append(startpoint[0])

            if startpoint == "prev endpoint":
                startpoint = endpoint_prev
            if startpoint_dydx == "prev endpoint":
                startpoint_dydx = endpoint_dydx_prev

            curve = Curve(type_, startpoint, xdata, ydata, startpoint_dydx, x_end, endpoint)
            x_endpoints.append(curve.x_bounds[-1])
            y_funcs.append(curve.y), dydx_funcs.append(curve.dydx)
            endpoint_prev = (curve.x_bounds[-1], curve.y(curve.x_bounds[-1]))
            endpoint_dydx_prev = curve.dydx(curve.x_bounds[-1])

        self.x_bounds = (min(x_endpoints), max(x_endpoints))
        self.merge_curve_functions(y_funcs, dydx_funcs)

    def merge_curve_functions(self, y_funcs, dydx_funcs):
        """
        returns functions defining y and dydx of the whole composite curve
        """
        self.y = lambda x: sum([func(x) for func in y_funcs if func(x) is not None])
        self.dydx = lambda x: sum([func(x) for func in dydx_funcs if func(x) is not None])
