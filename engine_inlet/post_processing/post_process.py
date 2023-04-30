import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.tri
import numpy as np 
import math 
"""
module responsible for generating plots and figures
"""
class create_slice_plot:
    
    def __init__(self, plotDict, mainObj, plotSettings):
        
        xlims, ylims = None, None 
        if "xlim" in plotSettings.keys():
            if plotSettings["xlim"] not in ["none", "None", "NONE"]: xlims = plotSettings["xlim"]
        if "ylim" in plotSettings.keys():
            if plotSettings["ylim"] not in ["none", "None", "NONE"]: ylims = plotSettings["ylim"]

        self.M_inf = mainObj.freestream.mach
        figsize = plotSettings["figure size"]
        fig, ax = self.initialize_figure(figsize, xlims=xlims, ylims=ylims)
        #pos = [0.05, 1] #placeholder 
        #self.display_sol_params(ax, pos, mainObj)
        self.plot_from_plotDict(plotDict, ax, mainObj) #

    def initialize_figure(self, figsize, xlims=None, ylims=None): 
        fig = plt.figure(figsize=figsize) #create figure object
        ax = fig.add_subplot(1,1,1) 
        if xlims is not None: ax.set_xlim(xlims[0], xlims[-1])
        if ylims is not None: ax.set_ylim(ylims[0], ylims[-1])
        ax.set_xlabel('x'), ax.set_ylabel('y', rotation='horizontal'), ax.grid(linewidth=0.3, color='grey')
        return fig, ax 

    def display_sol_params(self, axes, pos, mainObj):
        textStr = f"Freestream:\n M = {mainObj.freestream.mach}\n T = {round(mainObj.freestream.T,2)} K \np = {round(mainObj.freestream.p,1)} Pa"
        axes.text(pos[0], pos[1], textStr)

    def plot_from_plotDict(self, plotDict, axes, mainObj):
        """
        Loads input file and converts entries to object attributes
        """
        for key in plotDict.keys(): 

            typ = plotDict[key]["type"]

            if typ == "geom":
                self.plot_inletGeom(axes, mainObj.inputs.geom)

            elif typ == "mesh":

                anno, mFlow, wall_flow = False, False, None
                if plotDict[key]["annotate"] == True: anno=True
                if plotDict[key]["show mass flow"] == True: mFlow=True 
                if "wall flow scalar" in plotDict[key].keys(): 
                    if plotDict[key]["wall flow scalar"] not in ["None", "none", "NONE"]:
                        wall_flow = plotDict[key]["wall flow scalar"]

                self.plot_coneSol(axes, mainObj.coneSol, mainObj.inputs.geom)
                self.plot_mesh(axes, mainObj.mesh, annotate=anno, mass_flow_plot=mFlow, wall_flow_plot=wall_flow)
                self.plot_idl(axes, mainObj.idlObj)


            elif typ == "scalar":
                self.plot_coneSol(axes, mainObj.coneSol, mainObj.inputs.geom)
                scalar = plotDict[key]["scalar"]
                lims = plotDict[key]["limits"]
                #zs = [getattr(pt, scalar) for pt in mainObj.mesh.meshPts]
                #self.plot_scalar_contours(axes, scalar, lims, idl=mainObj.idlObj, mesh=mainObj.mesh, coneSol=mainObj.coneSol, freeStream = mainObj.freestream)
                self.plot_scalar_contours(axes, scalar, lims, idl=mainObj.idlObj, mesh=mainObj.mesh) #dumbed down to not include freestream 

            else: 
                raise ValueError(f"invalid displayer type: {typ}")

    def plot_coneSol(self, axes, cone, inletGeom):
        xint = np.array([0, 1])
        if inletGeom is not None:
            #plot interval only which conforms to inlet geometry 
            xint = np.array([min(inletGeom.centerbody_bounds), max(inletGeom.centerbody_bounds)])
        else: 
            axes.plot(xint, [x*math.tan(cone.cone_ang) for x in xint], label=f'cone = {round(math.degrees(cone.cone_ang),2)}', color='w', linewidth=1.3) #plot straight cone surface

        axes.plot(xint, [x*math.tan(cone.shock_ang) for x in xint], label=f'shock = {round(math.degrees(cone.shock_ang),2)} deg', color='crimson', linewidth=2, linestyle='dashdot') 

    def plot_inletGeom(self, axes, inletGeom):
        #plot inlet geometry: 
        x_cowl = np.linspace(inletGeom.cowl_bounds[0], inletGeom.cowl_bounds[1], 100)
        axes.plot(x_cowl, [inletGeom.y_cowl(x) for x in x_cowl], '-w', linewidth=2)
        x_cb = np.linspace(inletGeom.centerbody_bounds[0], inletGeom.centerbody_bounds[1], 100)
        #axes.plot(x_cb, [inletGeom.y_centerbody(x) for x in x_cb], '-w', linewidth=2)
        axes.axhline(0, color='w', linestyle='dashed', linewidth=1) 

        fill_x = np.array([max(x_cb)])
        fill_x = np.append(fill_x, x_cb)
        fill_y = np.array([0])
        fill_y = np.append(fill_y, [inletGeom.y_centerbody(x) for x in x_cb])
        axes.fill(fill_x, fill_y, facecolor="black", edgecolor="white", zorder=15, hatch="\\\\", linewidth=2) 
         
    def plot_idl(self, axes, idl, annotate=None): 
        axes.plot(idl.x, idl.y, '-o', linewidth=0.5, markersize=2, color='aquamarine')
        for i,x in enumerate(idl.x): 
            axes.plot([0,x],[0,idl.y[i]],linewidth=0.5,color='aquamarine')
        if annotate: 
            for i,x in enumerate(idl.x):
                text = f"V={round(idl.u[i],1)}, {round(idl.v[i],1)}"
                xy = (x,idl.y[i])
                axes.annotate(text, xy)

    def plot_mesh(self, axes, mesh, annotate=False, mass_flow_plot=False, wall_flow_plot=None):
        
        
        axes.scatter([pt.x for pt in mesh.meshPts],[pt.y for pt in mesh.meshPts], color='aquamarine', s=2)

        if annotate: 
            [axes.annotate(f"{pt.i}", (pt.x,pt.y)) for pt in mesh.meshPts]
                
        for tri in mesh.triangle:
            a,b,c = tri
            if a is None: 
                continue
            if b is not None: 
                plt.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[1]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[1]].y], color='aquamarine', linewidth=0.5)
            if c is not None:
                plt.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[2]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[2]].y], color='aquamarine', linewidth=0.5)

        if hasattr(mesh, "shock_segs"):
            for i,ind in enumerate(mesh.shock_segs):
                if i != 0: 
                    prevPt = mesh.meshPts[mesh.shock_segs[i-1]]
                    pt = mesh.meshPts[ind]
                    plt.plot([pt.x, prevPt.x],[pt.y, prevPt.y], color='crimson', linewidth=2, linestyle='dashdot')

        if mass_flow_plot and hasattr(mesh, 'char_mass_flow'): 
            fig2 = plt.figure()
            ax2 = fig2.add_subplot(1,1,1)
            ax2.plot(mesh.char_mass_flow[0], mesh.char_mass_flow[1], 'o-')
            ax2.set_xlabel("characteristic line no."), ax2.set_ylabel("mass flow rate (kg/sec)")
            ax2.grid(linewidth=0.3, color='grey')
            ax2.set_xlim(0, max(mesh.char_mass_flow[0]))

        if wall_flow_plot is not None: 

            self.plot_surface_properties(mesh, wall_flow_plot)

    def plot_surface_properties(self, mesh, scalar):
        """
        plots scalar properties on the wall. Currently works for pressure only
        """        
        fig = plt.figure(figsize=(16,8)) #create figure object
        ax1 = fig.add_subplot(1,1,1) 
        ax1.set_xlim(0,4.3)
        
        if scalar == "p/p_inf":

            p_inf = mesh.gasProps.p0/((1 +  0.5*(mesh.gasProps.gam-1)*self.M_inf**2)**(mesh.gasProps.gam/(mesh.gasProps.gam-1)))
            ax1.plot([pt.x for pt in mesh.wallPtsUpper],[pt.p/p_inf for pt in mesh.wallPtsUpper], color='r', label="cowl")
            ax1.plot([pt.x for pt in mesh.wallPtsLower],[pt.p/p_inf for pt in mesh.wallPtsLower], color='b', label="centerbody")
            ax1.set_xlabel('x'), ax1.set_ylabel('p/p_inf'), ax1.grid(linewidth=0.3, color='grey')
            ax1.legend()
            ax1.set_ylim(0, 10)
        else: raise ValueError(f"unknown wall scalar: {scalar}")

    def plot_scalar_contours(self, axes, scalar, lims, idl=None, coneSol=None, mesh=None, freeStream=None, barLabel=None,):
        """
        TODO: freestream not plotting right. Screws up the blending
        """
        xList = []
        yList = []
        scalarList = []

        if idl is not None: 
            #plot tmc region upto idl 
            x_init,y_init = [],[]
            r = 0.00001
            for i,x in enumerate(idl.x):
                thet = math.atan(idl.y[i]/x)
                x_init.append(r*math.cos(thet))
                y_init.append(r*math.sin(thet))

            xList += idl.x + x_init
            yList += idl.y + y_init
            scalarList = scalarList + getattr(idl, scalar) + getattr(idl,scalar)

        if mesh is not None: 
            #plot mesh region

            #sort mesh points by regions
            mesh_point_regions = [[]]
            for pt in mesh.meshPts: 
                reg_ind = pt.reg
                while reg_ind > len(mesh_point_regions)-1: 
                    mesh_point_regions.append([])
                mesh_point_regions[reg_ind].append(pt)
                
            for pts in mesh_point_regions: 
                xList += [pt.x for pt in pts]
                yList += [pt.y for pt in pts]
                scalarList = scalarList + [getattr(pt, scalar) for pt in pts] 
                mocReg = matplotlib.tri.Triangulation(xList,yList) 
                axes.tricontourf(mocReg, scalarList, 100, cmap='jet', vmin=lims[0], vmax=lims[1])
                scalarList = [] 
                xList = []
                yList = []

        #plot far field triangle
        if freeStream is not None and coneSol is not None: 
            #Testing this out
            xpts = [0, 2, 0]
            ypts = [0, 1, 1]
            z = getattr(freeStream, scalar)
            scalarList = [z,z,z]
            freestrReg = matplotlib.tri.Triangulation(xpts, ypts)
            axes.tricontourf(freestrReg, scalarList, 100, cmap='jet', vmin=lims[0], vmax=lims[1])

        map_ = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=lims[0], vmax=lims[1]), cmap='jet')
        if barLabel is None: barLabel = scalar
        plt.colorbar(mappable=map_, orientation='horizontal', shrink=0.5, label=barLabel)
