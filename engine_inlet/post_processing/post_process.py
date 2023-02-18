import matplotlib.pyplot as plt
import numpy as np 
import math 
"""
module responsible for generating plots and figures
"""
class create_slice_plot:
    def __init__(self, coneSol=None, inletGeom=None, idl=None, annotateIdl=None, mesh=None):
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
            self.plot_idl(idl, ax, annotate=annotateIdl)
        if mesh is not None: 
            self.plot_mesh(mesh, ax)

        ax.legend()
        self.fig = fig
        plt.show()

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
        axes.plot(idl.x, idl.y, '-o', linewidth=0.5, markersize=2, color='gold')
        if annotate: 
            for i,x in enumerate(idl.x):
                text = f"V={round(idl.u[i],1)}, {round(idl.v[i],1)}"
                xy = (x,idl.y[i])
                axes.annotate(text, xy)

    def plot_mesh(self, mesh, axes):
        axes.scatter([pt.x for pt in mesh.meshPts],[pt.y for pt in mesh.meshPts], color='gold', s=2)
        for tri in mesh.triangle:
            a,b,c = tri
            if b is not None: 
                plt.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[1]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[1]].y], color='gold', linewidth=0.5)
            if c is not None:
                plt.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[2]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[2]].y], color='gold', linewidth=0.5)