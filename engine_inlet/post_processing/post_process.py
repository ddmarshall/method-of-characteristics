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

        self.M_inf = mainObj.inputs.freeStream.mach
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
                
                if hasattr(mainObj, "coneSol"): 
                    self.plot_coneSol(axes, mainObj.coneSol, mainObj.inputs.geom)  
                    
                elif hasattr(mainObj, "rampSol"):
                    self.plot_rampSol(axes, mainObj.rampSol, mainObj.inputs.geom)

                self.plot_mesh(axes, mainObj.mesh, annotate=anno, mass_flow_plot=mFlow, wall_flow_plot=wall_flow)
                if mainObj.inputs.init_method in ["STRAIGHT IDL", "MACH LINE"]:
                    self.plot_idl(axes, mainObj.idlObj, mainObj.inputs.delta)


            elif typ == "scalar":
                if hasattr(mainObj, "coneSol"): self.plot_coneSol(axes, mainObj.coneSol, mainObj.inputs.geom)
                elif hasattr(mainObj, "rampSol"): self.plot_rampSol(axes, mainObj.rampSol, mainObj.inputs.geom)
                scalar = plotDict[key]["scalar"]
                lims = plotDict[key]["limits"]
                #zs = [getattr(pt, scalar) for pt in mainObj.mesh.meshPts]
                #self.plot_scalar_contours(axes, scalar, lims, idl=mainObj.idlObj, mesh=mainObj.mesh, coneSol=mainObj.coneSol, freeStream = mainObj.freestream)
                self.plot_scalar_contours(axes, scalar, lims, idl=mainObj.idlObj, mesh=mainObj.mesh) #dumbed down to not include freestream 

            else: 
                raise ValueError(f"invalid displayer type: {typ}")

    def plot_coneSol(self, axes, cone_flow, inletGeom):
        xint = np.array([0, 1])
        if inletGeom is not None:
            #plot interval only which conforms to inlet geometry 
            xint = np.array([min(inletGeom.centerbody_bounds), max(inletGeom.centerbody_bounds)])
        else: 
            axes.plot(xint, [x*math.tan(cone_flow.cone_ang) for x in xint],\
                label=f'cone = {round(math.degrees(cone_flow.cone_ang),2)}',\
                    color='k', linewidth=1.3) #plot straight cone surface

        axes.plot(xint, [x*math.tan(cone_flow.shock_ang) for x in xint],\
                label=f'shock = {round(math.degrees(cone_flow.shock_ang),2)} \
                    deg', color='crimson', linewidth=2, linestyle='dashdot') 

    def plot_rampSol(self, axes, ramp_flow, inletGeom):
        xint = np.array([0,1])
        if inletGeom is not None: 
            xint = np.array([min(inletGeom.centerbody_bounds), max(inletGeom.centerbody_bounds)])
        else: 
            axes.plot(xint, [x*math.tan(ramp_flow.deflec) for x in xint],\
                label=f'cone = {round(math.degrees(ramp_flow.deflec),2)}',\
                    color='k', linewidth=1.3) #plot straight cone surface

        axes.plot(xint, [x*math.tan(ramp_flow.beta) for x in xint],\
                label=f'shock = {round(math.degrees(ramp_flow.beta),2)} \
                    deg', color='crimson', linewidth=2, linestyle='dashdot')

    def plot_inletGeom(self, axes, inletGeom):
        #plot inlet geometry: 

        #line_color = "white"
        line_color = "white"
        face_color="black"

        x_cowl = np.linspace(inletGeom.cowl_bounds[0], inletGeom.cowl_bounds[1], 100)
        axes.plot(x_cowl, [inletGeom.y_cowl(x) for x in x_cowl], color=line_color, linewidth=2)
        x_cb = np.linspace(inletGeom.centerbody_bounds[0], inletGeom.centerbody_bounds[1], 100)
        axes.plot(x_cb, [inletGeom.y_centerbody(x) for x in x_cb], color=line_color, linewidth=2)
        axes.axhline(0, color=line_color, linestyle='dashed', linewidth=1) 
        fill_x = np.array([max(x_cb)])
        fill_x = np.append(fill_x, x_cb)
        fill_y = np.array([0])
        fill_y = np.append(fill_y, [inletGeom.y_centerbody(x) for x in x_cb])
        axes.fill(fill_x, fill_y, facecolor=face_color, edgecolor=line_color, zorder=15, hatch="\\\\", linewidth=2) 
         
    def plot_idl(self, axes, idl, delta, annotate=None): 

        line_color = "aquamarine"

        axes.plot(idl.x, idl.y, '-o', linewidth=0.5, markersize=2, color=line_color)

        if delta == 1: 
            for i,x in enumerate(idl.x): 
                axes.plot([0,x],[0,idl.y[i]],linewidth=0.5,color=line_color)

        if annotate: 
            for i,x in enumerate(idl.x):
                text = f"V={round(idl.u[i],1)}, {round(idl.v[i],1)}"
                xy = (x,idl.y[i])
                axes.annotate(text, xy)

    def plot_mesh(self, axes, mesh, annotate=False, mass_flow_plot=False, wall_flow_plot=None):
        
        mesh_color = "aquamarine"
        #mesh_color = "dimgrey"
        axes.scatter([pt.x for pt in mesh.meshPts],[pt.y for pt in mesh.meshPts], color=mesh_color, s=2)
        
        if annotate: 
            [axes.annotate(f"{pt.i}", (pt.x,pt.y)) for pt in mesh.meshPts]
                
        for tri in mesh.triangle:
            a,b,c = tri
            if a is None: 
                continue
            if b is not None: 
                plt.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[1]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[1]].y], color=mesh_color, linewidth=0.5)
            if c is not None:
                plt.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[2]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[2]].y], color=mesh_color, linewidth=0.5)

        if hasattr(mesh, "shock_segs"):
            for i,ind in enumerate(mesh.shock_segs):
                if i != 0: 
                    prevPt = mesh.meshPts[mesh.shock_segs[i-1]]
                    pt = mesh.meshPts[ind]
                    plt.plot([pt.x, prevPt.x],[pt.y, prevPt.y], color='crimson', linewidth=2, linestyle='dashdot')

        if mass_flow_plot and hasattr(mesh, "mesh_mass_flow"): 
            fig2 = plt.figure(figsize=(5,4))
            ax2 = fig2.add_subplot(1,1,1)
            ax2.plot(mesh.mesh_mass_flow[0][0], mesh.mesh_mass_flow[0][1], "o-", label="+ Characteristics", markerfacecolor="none")
            ax2.plot(mesh.mesh_mass_flow[1][0], mesh.mesh_mass_flow[1][1], "o-", label="- Characteristics", markerfacecolor="none")
            mflows_max = max([max(mesh.mesh_mass_flow[0][1]), max(mesh.mesh_mass_flow[1][1])])
            mflows_min = min([min(mesh.mesh_mass_flow[0][1]), min(mesh.mesh_mass_flow[1][1])])
            mflows_range = mflows_max - mflows_min

            ax2.set_ylim(0.96, 1.02)
            ax2.axhline(mflows_max, linestyle="--", color="k"), ax2.axhline(mflows_min, linestyle="--", color="k")
            ax2.set_xlabel("mach line #"), ax2.set_ylabel("local mass flow ratio")
            ax2.set_title(f"range = {round(mflows_range,5)}")
            ax2.grid(linewidth=0.3, color='grey')
            ax2.legend()
            

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

class create_figure:
    """
    creates a single figure depending on a figure dictionary
    """
    def __init__(self, fig_dict, mainObj, plotSettings):

        self.set_default_settings(plotSettings)
        self.read_fig_dict(fig_dict)
        ax1, ax2 = self.initialize_figure()
        self.plot_inlet_geom(ax1, mainObj.inputs.geom)  

        if self.type == "mesh":
            self.plot_incident_shock_solution(ax1, mainObj) 
            self.plot_initial_data_line(ax1, mainObj)
            self.plot_mesh(ax1, mainObj.mesh, annotate=self.mesh_annotate)
            if hasattr(self, "mflow_plot"):
                self.plot_mass_flow_ratio(ax2, mainObj.mesh) 
        
        elif self.type == "scalar":
            #interpret scalar parameter
            scalar_interpreter = {
                "mach":                 "mach",
                "pressure ratio":       "p_p0",
                "temperature ratio":    "T_T0",
                "density ratio":        "rho_rho0"
            }            
            self.get_shock_endpoint(mainObj)
            self.scalar = scalar_interpreter[self.scalar_param]
            self.plot_scalar(ax1, mainObj)

            if self.surf_props:
                self.plot_surface_properties(ax2, mainObj.mesh)
            
    def set_default_settings(self, plotSettings):
        """
        sets initial object attributes
        """
        self.figname = None
        self.figsize = None
        self.xlim, self.ylim = None, None
        self.mflow_xlim, self.mflow_ylim = None, None 
        self.surfPlot_xlim, self.surfPlot_ylim = None, None 

        if "figure size" in plotSettings.keys():
            self.figsize = plotSettings["figure size"]

        if "xlim" in plotSettings.keys():
            self.xlim = plotSettings["xlim"]
        
        if "ylim" in plotSettings.keys():
            self.ylim = plotSettings["ylim"]
        
        if "style" in plotSettings.keys():
            plt.style.use(plotSettings["style"])

        if "mass flow xlim" in plotSettings.keys():
            self.mflow_xlim = plotSettings["mass flow xlim"]

        if "mass flow ylim" in plotSettings.keys():
            self.mflow_ylim = plotSettings["mass flow ylim"]

        if "surf props xlim" in plotSettings.keys():
            self.surfPlot_xlim = plotSettings["surf props xlim"]

        if "surf props ylim" in plotSettings.keys():
            self.surfPlot_ylim = plotSettings["surf props ylim"]

    def read_fig_dict(self, fig_dict):

        translator = {
            "figure name":              "figname",
            "figure size":              "figsize",
            "type":                     "type",
            "parameter":                "scalar_param",
            "show surface plots":       "surf_props",
            "colorbar limits":          "cbar_lims",
            "colorbar label":           "cbar_label",
            "annotate":                 "mesh_annotate",
            "show mass flow":           "mflow_plot",
            "geom x limits":            "xlim",
            "geom y limits":            "ylim",
            "mass flow x limits":       "mflow_xlim",
            "mass flow y limits":       "mflow_ylim",
            "surface plot x limits":    "surfPlot_xlim",
            "surface plot y limits":    "surfPlot_ylim"
        }

        #set attributes 
        for key in fig_dict.keys():
            if key in translator.keys():
                setattr(self, translator[key], fig_dict[key])

    def initialize_figure(self):
        
        self.fig = plt.figure(figsize=self.figsize)
        if hasattr(self, "figname"):
            self.fig.canvas.manager.set_window_title(self.figname)
        
        if hasattr(self, "mflow_plot") or self.surf_props:
            ax1 = self.fig.add_subplot(212)
            ax2 = self.fig.add_subplot(211)
            ax2.grid(linewidth=0.3, color='grey')
            #self.fig.tight_layout()

        else: 
            ax1 = self.fig.add_subplot(111)
            ax2 = None

        ax1.grid(linewidth=0.3, color='grey')

        if self.xlim is not None: 
            ax1.set_xlim(self.xlim[0], self.xlim[-1])
        if self.ylim is not None: 
            ax1.set_ylim(self.ylim[0], self.ylim[-1])

        return ax1, ax2 

    def plot_inlet_geom(self,ax, inletGeom):
        """
        plots the upper and lower cowl surfaces of the inlet
        """
        line_color = "white"
        face_color= "black"

        x_cowl = np.linspace(inletGeom.cowl_bounds[0], inletGeom.cowl_bounds[1], 100)
        ax.plot(x_cowl, [inletGeom.y_cowl(x) for x in x_cowl], color=line_color, linewidth=2)
        x_cb = np.linspace(inletGeom.centerbody_bounds[0], inletGeom.centerbody_bounds[1], 100)
        ax.plot(x_cb, [inletGeom.y_centerbody(x) for x in x_cb], color=line_color, linewidth=2)
        ax.axhline(0, color=line_color, linestyle='dashed', linewidth=1) 
        fill_x = np.array([max(x_cb)])
        fill_x = np.append(fill_x, x_cb)
        fill_y = np.array([0])
        fill_y = np.append(fill_y, [inletGeom.y_centerbody(x) for x in x_cb])
        ax.fill(fill_x, fill_y, facecolor=face_color, edgecolor=line_color, zorder=15, hatch="\\\\", linewidth=2) 

    def get_shock_endpoint(self, mainObj):
        """
        
        """
        xint = [min(mainObj.inputs.geom.centerbody_bounds), max(mainObj.inputs.geom.centerbody_bounds)]
        if hasattr(mainObj, "coneSol"):
            self.shock_endpoint = [xint[-1], math.tan(mainObj.coneSol.shock_ang)*xint[-1]]
        elif hasattr(mainObj, "rampSol"):
            self.shock_endpoint = [xint[-1], math.tan(mainObj.rampSol.beta)*xint[-1]]

    def plot_incident_shock_solution(self, ax, mainObj):
        """
        plots the incident shock wave on the ramp or cone incident point
        """
        inletGeom = mainObj.inputs.geom
        xint = np.array([min(inletGeom.centerbody_bounds), max(inletGeom.centerbody_bounds)])

        if hasattr(mainObj, "coneSol"):
            cone_flow = mainObj.coneSol
            
            #ax.plot(xint, [x*math.tan(cone_flow.cone_ang) for x in xint],\
            #    label=f'cone = {round(math.degrees(cone_flow.cone_ang),2)}',\
            #        color='k', linewidth=1.3) #plot straight cone surface

            ax.plot(xint, [x*math.tan(cone_flow.shock_ang) for x in xint],\
                label=f'shock = {round(math.degrees(cone_flow.shock_ang),2)} \
                    deg', color='crimson', linewidth=2, linestyle='dashdot') 

        elif hasattr(mainObj, "rampSol"):
            ramp_flow = mainObj.rampSol
            #ax.plot(xint, [x*math.tan(ramp_flow.deflec) for x in xint],\
            #    label=f'cone = {round(math.degrees(ramp_flow.deflec),2)}',\
            #        color='k', linewidth=1.3) #plot straight cone surface

            ax.plot(xint, [x*math.tan(ramp_flow.beta) for x in xint],\
                label=f'shock = {round(math.degrees(ramp_flow.beta),2)} \
                    deg', color='crimson', linewidth=2, linestyle='dashdot')

    def plot_initial_data_line(self, ax, mainObj):
        """
        docstring
        """
        idl = mainObj.idlObj
        line_color = "aquamarine"
        ax.plot(idl.x, idl.y, '-o', linewidth=0.5, markersize=2, color=line_color)
        
        if hasattr(mainObj, "coneSol") == 1: 
            for i,x in enumerate(idl.x): 
                ax.plot([0,x],[0,idl.y[i]],linewidth=0.5,color=line_color)
        
        elif hasattr(mainObj, "rampSol"):
            ramp = mainObj.rampSol
            for i,x_i in enumerate(idl.x):
                y_i = idl.y[i]
                y = (y_i - math.tan(ramp.deflec)*x_i)/(1 - (math.tan(ramp.deflec)/math.tan(ramp.beta)))
                x = y/math.tan(ramp.beta)
                ax.plot([x,x_i],[y,y_i], linewidth=0.5, color=line_color)

    def plot_mesh(self, ax, mesh, annotate=False):
        """

        """
        mesh_color = "aquamarine"
        #mesh_color = "dimgrey"
        ax.scatter([pt.x for pt in mesh.meshPts],[pt.y for pt in mesh.meshPts], color=mesh_color, s=2)
        
        if annotate: 
            [ax.annotate(f"{pt.i}", (pt.x,pt.y)) for pt in mesh.meshPts]
                
        for tri in mesh.triangle:
            a,b,c = tri
            if a is None: 
                continue
            if b is not None: 
                ax.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[1]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[1]].y], color=mesh_color, linewidth=0.5)
            if c is not None:
                ax.plot([mesh.meshPts[tri[0]].x, mesh.meshPts[tri[2]].x],[mesh.meshPts[tri[0]].y, mesh.meshPts[tri[2]].y], color=mesh_color, linewidth=0.5)

        if hasattr(mesh, "shock_segs"):
            for i,ind in enumerate(mesh.shock_segs):
                if i != 0: 
                    prevPt = mesh.meshPts[mesh.shock_segs[i-1]]
                    pt = mesh.meshPts[ind]
                    ax.plot([pt.x, prevPt.x],[pt.y, prevPt.y], color='crimson', linewidth=2, linestyle='dashdot')

    def plot_scalar(self, ax, mainObj):
        """
        
        """
        if hasattr(self, "cbar_lims"):
            lims = self.cbar_lims
        else:
            if self.scalar in ["p_p0", "T_T0", "rho_rho0"]:
                lims = [0,1]
            elif self.scalar == "mach":
                lims = [1, mainObj.inputs.M_inf]

        #plot region past incident shock but before idl region 
        data_ups = mainObj.upstream_data
        x_reg, y_reg, scal_reg = [],[],[]
        r = 1e-5
        for i,x in enumerate(data_ups.x):
            thet = math.atan(data_ups.y[i]/x)
            x_reg.append(r*math.cos(thet))
            y_reg.append(r*math.sin(thet))
            scal_reg.append(getattr(data_ups, self.scalar)[i])
        x_reg += data_ups.x
        y_reg += data_ups.y
        scal_reg += getattr(data_ups, self.scalar)
        tri_reg = matplotlib.tri.Triangulation(x_reg, y_reg)
        ax.tricontourf(tri_reg, scal_reg, 100, cmap='jet', vmin=lims[0],vmax=lims[1])

        #plot region from tip of cone to initial data line 
        idl = mainObj.idlObj
        x_reg, y_reg, scal_reg = [],[],[]
        for i,x in enumerate(idl.x):
            thet = math.atan(idl.y[i]/x)
            x_reg.append(r*math.cos(thet))
            y_reg.append(r*math.sin(thet))
            scal_reg.append(getattr(idl, self.scalar)[i])
        x_reg += idl.x 
        y_reg += idl.y
        scal_reg += getattr(idl, self.scalar)
        tri_reg = matplotlib.tri.Triangulation(x_reg, y_reg)
        ax.tricontourf(tri_reg, scal_reg, 100, cmap='jet', vmin=lims[0],vmax=lims[1])
            
        #plot characteristic mesh region
        xList, yList, scalarList = [], [], []
        mesh_point_regions = [[]]
        for pt in mainObj.mesh.meshPts: 
            reg_ind = pt.reg
            while reg_ind > len(mesh_point_regions)-1: 
                mesh_point_regions.append([])
            mesh_point_regions[reg_ind].append(pt)
            
        for pts in mesh_point_regions: 
            xList += [pt.x for pt in pts]
            yList += [pt.y for pt in pts]
            for pt in pts: 
                obj = pt
                scalarList = scalarList + [getattr(pt, self.scalar)]

            mocReg = matplotlib.tri.Triangulation(xList,yList) 
            ax.tricontourf(mocReg, scalarList, 100, cmap='jet', vmin=lims[0], \
                           vmax=lims[1])
            scalarList = [] 
            xList = []
            yList = []

        #plot far field triangle
        xpts = [0, self.shock_endpoint[0], 0]
        ypts = [0, self.shock_endpoint[-1], self.shock_endpoint[-1]]
        z = getattr(mainObj.inputs.freeStream, self.scalar)
        scalarList = [z,z,z]
        freestrReg = matplotlib.tri.Triangulation(xpts, ypts)
        ax.tricontourf(freestrReg, scalarList, 100, cmap='jet', vmin=lims[0], \
                       vmax=lims[1])

        map_ = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(\
            vmin=lims[0], vmax=lims[1]), cmap='jet')
        
        if hasattr(self, "cbar_label"): barLabel = self.cbar_label
        else: barLabel = self.scalar
        plt.colorbar(mappable=map_, orientation='horizontal', shrink=0.5, \
                     label=barLabel)

    def plot_surface_properties(self, ax, mesh):
        """
        
        """
        ax.plot([pt.x for pt in mesh.wallPtsUpper],[getattr(pt,self.scalar) \
                                    for pt in mesh.wallPtsUpper], label="cowl")
        ax.plot([pt.x for pt in mesh.wallPtsLower],[getattr(pt,self.scalar) \
                                    for pt in mesh.wallPtsLower], label="centerbody")
        ax.set_xlabel('x'), ax.set_ylabel(f'{self.scalar}'), \
                                            ax.grid(linewidth=0.3, color='grey')
        ax.set_xlim(self.xlim[0], self.xlim[-1])
        ax.legend()
        
    def plot_mass_flow_ratio(self, ax, mesh):
        """
        plots local mass flow ratio throughout the mesh. Used for checking order
        of accuracy of solution
        """
        ax.plot(mesh.mesh_mass_flow[0][0], mesh.mesh_mass_flow[0][1], "o-", label="+ Characteristics", markerfacecolor="none")
        ax.plot(mesh.mesh_mass_flow[1][0], mesh.mesh_mass_flow[1][1], "o-", label="- Characteristics", markerfacecolor="none")
        mflows_max = max([max(mesh.mesh_mass_flow[0][1]), max(mesh.mesh_mass_flow[1][1])])
        mflows_min = min([min(mesh.mesh_mass_flow[0][1]), min(mesh.mesh_mass_flow[1][1])])
        mflows_range = mflows_max - mflows_min
        #ax.set_ylim(0.96, 1.02)
        ax.axhline(mflows_max, linestyle="--", color="white"), ax.axhline(mflows_min, linestyle="--", color="white")
        ax.set_xlabel("mach line #"), ax.set_ylabel("local mass flow ratio")
        ax.set_title(f"range = {round(mflows_range,5)}")
        ax.grid(linewidth=0.3, color='grey')
        ax.legend()