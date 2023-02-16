import matplotlib.pyplot as plt
import numpy as np 
import math 
import os 
import sys 
sys.path.append(os.getcwd() + "\\taylor_maccoll_cone") #add path to taylor maccoll module
import taylor_maccoll as tmc
"""
responsible for generating plots and figures
"""

class create_slice_plot:
    def __init__(self, coneSol=None, inletGeom=None, idl=None, mesh=None):
        #create plot depending on which objects you give it
        plt.style.use('dark_background') #dark > light mode 
        fig = plt.figure(figsize=(16,9)) #create figure object
        ax = fig.add_subplot(1,1,1)
        ax.set_xlabel('x'), ax.set_ylabel('y'), ax.grid(linewidth=0.3, color='grey')

        if coneSol is not None:
            self.plot_coneSol(coneSol, ax, inletGeom) 
        if inletGeom is not None: 
            self.plot_inletGeom(inletGeom, ax)
        if idl is not None: 
            self.plot_idl(idl, ax)
        if mesh is not None: 
            self.plot_mesh(mesh, ax)

        ax.legend()
        self.fig = fig

    def plot_coneSol(self, cone, axes, inletGeom):
        axes.set_title(f"M = {cone.M_inf}, \u03B3 = {cone.gam}, R = {cone.R} J/(kg*K), T_0 = {cone.T0} K") 
        xint = np.array([0, 1])
        if inletGeom is not None:
            #plot incident only which conforms to inlet geometry 
            xint = np.array([min(inletGeom.centerbody_bounds), max(inletGeom.centerbody_bounds)])
        else: 
            axes.plot(xint, [x*math.tan(cone.cone_ang) for x in xint], label=f'cone = {round(math.degrees(cone.cone_ang),2)}', color='w', linewidth=1.3) #plot straight cone surface

        axes.plot(xint, [x*math.tan(cone.shock_ang) for x in xint], label=f'shock = {round(math.degrees(cone.shock_ang),2)} deg', color='r', linewidth=0.7) 

    def plot_inletGeom(self, inletGeom, axes):
        #plot inlet geometry: 
        x_cowl = np.linspace(inletGeom.cowl_bounds[0], inletGeom.cowl_bounds[1], 1000)
        axes.plot(x_cowl, [inletGeom.y_cowl(x) for x in x_cowl], '-w', linewidth=1.3)
        x_cb = np.linspace(inletGeom.centerbody_bounds[0], inletGeom.centerbody_bounds[1], 1000)
        axes.plot(x_cb, [inletGeom.y_centerbody(x) for x in x_cb], '-w', linewidth=1.3)
        axes.axhline(0, color='w', linestyle='dashdot', linewidth=1) 
         
    def plot_idl(self, idl, axes, annotate=None): 
        axes.plot(idl.x_idl, idl.y_idl, '-o', label="idl", linewidth=0.5, markersize=2, color='gold')
        if annotate: 
            for i,x in enumerate(idl.x_idl):
                text = f"V={round(idl.u_idl[i],1)}, {round(idl.v_idl[i],1)}"
                xy = (x,idl.y_idl[i])
                axes.annotate(text, xy)

    def plot_mesh(self, mesh, axes):
        #TODO write when mesh generator is functional
        pass  

if __name__ == "__main__":
    #Testing out class
    #LOAD IN INLET GEOMETRY 
    import example_geometry as geom
    inlet = geom.inletGeom() 

    #CONE SOLUTION 
    gam = 1.4
    cone_ang = math.radians(12.5)
    M_inf = 2.5
    R = 287.05
    T0 = 288.15
    cone = tmc.TaylorMaccoll_Cone(cone_ang, M_inf, gam, R, T0) 

    #IDL 
    import idl as IDL
    class make_curve:
        def __init__(self, y_x, dist, endpoints):
            self.y_x, self.dist, self.endpoints = y_x, dist, endpoints
    #dist = [0, 0.1, 0.2, 0.30, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    dist = [0, 0.2, 0.4, 0.6, 0.8, 1]
    curve =  make_curve(lambda x: 4*(x-2.5)**2, dist, (2.01,2.15))
    idlObj = IDL.generate_tmc_initial_data_line(cone, curve)

    #MAKE PLOT 
    plotObj = create_slice_plot(coneSol=cone, inletGeom=inlet, idl=idlObj)
    plotObj.fig.show()
    pass             