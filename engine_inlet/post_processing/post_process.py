import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.tri
import numpy as np 
import math 
"""
module responsible for generating plots and figures
"""
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
            if self.mflow_plot:
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

        ax1, ax2 = None,None
        if hasattr(self, "mflow_plot"):
            if self.mflow_plot: 
                ax1 = self.fig.add_subplot(212) 
                ax2 = self.fig.add_subplot(211)
                ax2.grid(linewidth=0.3, color='grey')
        elif hasattr(self, "surf_props"):
            if self.surf_props: 
                ax1 = self.fig.add_subplot(212)
                ax2 = self.fig.add_subplot(211)
                ax2.grid(linewidth=0.3, color='grey')            
        
        if ax1 is None: 
            ax1 = self.fig.add_subplot(111)

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
        plt.colorbar(mappable=map_, orientation='horizontal', shrink=0.5, ax=ax, \
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