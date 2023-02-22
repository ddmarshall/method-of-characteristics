import matplotlib.pyplot as plt
import numpy as np 
import math 
"""
module responsible for generating plots and figures
"""
class create_slice_plot:
    def __init__(self, plotDict, defSettDict, mainObj):
        
        fig,ax = self.initialize_figure() #TODO take in default settings dictionary
        self.plot_from_plotDict(plotDict, ax, mainObj) #

    def initialize_figure(self, settingDict=None): 
        fig = plt.figure(figsize=(16,6)) #create figure object
        ax = fig.add_subplot(1,1,1) 
        ax.set_ylim(0,1.25)
        ax.set_xlabel('x'), ax.set_ylabel('y'), ax.grid(linewidth=0.3, color='grey')
        return fig, ax 

    def plot_from_plotDict(self, plotDict, axes, mainObj):
        """
        Loads input file and converts entries to object attributes
        TODO write this
        """
        for key in plotDict.keys(): 

            typ = plotDict[key]["type"]

            if typ == "geom":
                self.plot_inletGeom(axes, mainObj.inputs.geom)

            elif typ == "mesh":
                anno=False
                if plotDict[key]["annotate"] == True: anno=True
                self.plot_mesh(axes, mainObj.mesh, annotate=anno)
                self.plot_idl(axes, mainObj.idlObj)

            elif typ == "scalar":
                scalar = plotDict[key]["scalar"]
                self.plot_scalar_contours(axes, scalar, mainObj.idlObj, mainObj.mesh)

            else: 
                raise ValueError(f"invalid displayer type: {typ}")

    def plot_coneSol(self, cone, axes, inletGeom):
        axes.set_title(f"M = {cone.M_inf}, \u03B3 = {cone.gam}, R = {cone.R} J/(kg*K), T_0 = {cone.T0} K") 
        xint = np.array([0, 1])
        if inletGeom is not None:
            #plot interval only which conforms to inlet geometry 
            xint = np.array([min(inletGeom.centerbody_bounds), max(inletGeom.centerbody_bounds)])
        else: 
            axes.plot(xint, [x*math.tan(cone.cone_ang) for x in xint], label=f'cone = {round(math.degrees(cone.cone_ang),2)}', color='w', linewidth=1.3) #plot straight cone surface

        axes.plot(xint, [x*math.tan(cone.shock_ang) for x in xint], label=f'shock = {round(math.degrees(cone.shock_ang),2)} deg', color='red', linewidth=0.7) 

    def plot_inletGeom(self, axes, inletGeom):
        #plot inlet geometry: 
        x_cowl = np.linspace(inletGeom.cowl_bounds[0], inletGeom.cowl_bounds[1], 1000)
        axes.plot(x_cowl, [inletGeom.y_cowl(x) for x in x_cowl], '-w', linewidth=1.3)
        x_cb = np.linspace(inletGeom.centerbody_bounds[0], inletGeom.centerbody_bounds[1], 1000)
        axes.plot(x_cb, [inletGeom.y_centerbody(x) for x in x_cb], '-w', linewidth=1.3)
        axes.axhline(0, color='w', linestyle='dashed', linewidth=1) 
         
    def plot_idl(self, axes, idl, annotate=None): 
        axes.plot(idl.x, idl.y, '-o', linewidth=0.5, markersize=2, color='gold')
        for i,x in enumerate(idl.x): 
            axes.plot([0,x],[0,idl.y[i]],linewidth=0.5,color='gold')
        if annotate: 
            for i,x in enumerate(idl.x):
                text = f"V={round(idl.u[i],1)}, {round(idl.v[i],1)}"
                xy = (x,idl.y[i])
                axes.annotate(text, xy)

    def plot_mesh(self, axes, mesh, annotate=False):
        
        axes.scatter([pt.x for pt in mesh.meshPts],[pt.y for pt in mesh.meshPts], color='gold', s=2)
        if annotate: 
            [axes.annotate(f"{pt.i}", (pt.x,pt.y)) for pt in mesh.meshPts]
                
        for tri in mesh.triangle:
            _,b,c = tri
            if b is not None: 
                plt.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[1]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[1]].y], color='gold', linewidth=0.5)
            if c is not None:
                plt.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[2]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[2]].y], color='gold', linewidth=0.5)

    def plot_scalar_contours(self, axes, scalar, idl=None, mesh=None, barLabel=None):

        xList = []
        yList = []
        scalarList = []

        if idl is not None: 

            x_init,y_init = [],[]
            r = 0.00001
            for i,x in enumerate(idl.x):
                thet = math.atan(idl.y[i]/x)
                x_init.append(r*math.cos(thet))
                y_init.append(r*math.sin(thet))

            xList = xList + idl.x + x_init
            yList = yList + idl.y + y_init
            scalarList = scalarList + getattr(idl, scalar) + getattr(idl,scalar)

        if mesh is not None: 
            xList = xList + [pt.x for pt in mesh.meshPts]
            yList = yList + [pt.y for pt in mesh.meshPts]
            scalarList = scalarList + [getattr(pt, scalar) for pt in mesh.meshPts]  

        tcf = axes.tricontourf(xList, yList, scalarList, 100, cmap='jet')
        if barLabel is None: barLabel = scalar
        plt.colorbar(tcf, orientation='horizontal', shrink=0.5, label=barLabel)