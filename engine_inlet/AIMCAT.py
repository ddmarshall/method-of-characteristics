"""
***MAIN RUN FILE***
Author: Shay Takei 
"""
class Main:
    """
    This class controls all modules of AIMCAT.
    """
    def __init__(self, inputFile:str, geomFile:str, export:bool=False, \
                 plotFile:str=None, preview_geom:bool=False) -> None:
       
        self.load_inputs(inputFile, geomFile) #generate input object 
        if preview_geom: 
            import post_processing.post_process as post_process
            import json
            try: plotDict = json.load(open(plotFile, 'r'))
            except: plotDict = json.load(open("post_processing/"+plotFile, 'r'))
            plotSettings = plotDict["global plot settings"]
            post_process.Preview_Geom(self, plotSettings)

        self.run_solution() #run solution 
        self.print_details() #prints important details to the console

        if export: 
            self.export_results()

        if plotFile is not None: #if plotfile is provided, run it 
            self.plot_solution(plotFile)

    def load_inputs(self, inpFile:str, geomFile:str) -> None:
        import Input.input as inp
        import math 
        inpObj = inp.Input(inpFile, geomFile)
        self.inputs = inpObj

        #Kind of a dirty solution, but adds freestream object for plotting purposes
        class freeStream:
            def __init__(frst, inputObj):
                M, T0, gam, R = inputObj.M_inf, inputObj.T0, inputObj.gam, inputObj.R
                frst.mach = M
                frst.p_p0f = (1 + 0.5*(gam - 1)*M**2)**(-gam/(gam-1))
                frst.T = T0/(1 + 0.5*(gam - 1)*M**2)
                frst.T_T0 = frst.T/T0
                frst.rho_rho0f = (1 + 0.5*(gam - 1)*M**2)**(-1/(gam-1))
                a = math.sqrt(gam*R*frst.T)
                frst.u = M*a
                frst.v = 0 #!change if want angled flow for 2d case
                frst.V = math.sqrt(frst.u**2 + frst.v**2)  

        self.inputs.freeStream = freeStream(inpObj)

    def run_solution(self) -> None:
        print("\nrunning solution...")
        import method_of_characteristics.moc_mesh_engine as moc
        import initial_data_line.idl as idl
        import math
        import time

        t0 = time.perf_counter() #intial time
        inp = self.inputs
        class gasProps:
            def __init__(self, gam, R, T0): 
                self.gam, self.R, self.T0 = gam, R, T0
                self.a0 = math.sqrt(gam*R*T0)
        inp.gasProps = gasProps(inp.gam, inp.R, inp.T0) #create gas properties object 

        #RUNNING TAYLOR-MACCOLL or 2D WEDGE SOLUTION: 
        if inp.delta==1: #cone
            import taylor_maccoll_cone.taylor_maccoll as tmc 
            self.coneSol = tmc.TaylorMaccoll_Cone(math.radians(inp.geom.init_turn_ang_deg), inp.M_inf, inp.gasProps)
            #check if incident shock crosses centerbody geometry:
            if inp.geom.y_cowl(inp.geom.x_cowl_lip) > math.tan(self.coneSol.shock_ang)*inp.geom.x_cowl_lip:
                #!assumes straight incident shock up to cowl lip  
                raise ValueError("Incident Shock Crosses Cowl Geometry. Solution Cannot Proceed.")

        elif inp.delta==0: #wedge
            import method_of_characteristics.oblique_shock as shock 
            deflec = math.radians(inp.geom.init_turn_ang_deg)
            self.rampSol = shock.Oblique_Shock(inp.M_inf, inp.gasProps.gam, inp.gasProps.R, deflec=deflec)
            #check if incident shock crosses centerbody geometry:
            if inp.geom.y_cowl(inp.geom.x_cowl_lip) > math.tan(self.rampSol.beta)*inp.geom.x_cowl_lip:
                #!assumes straight incident shock up to cowl lip  
                raise ValueError("Incident Shock Crosses Cowl Geometry. Solution Cannot Proceed.")
            
        #GENERATING INITIAL DATA LINE
        if inp.init_method == "STRAIGHT IDL":
            if inp.delta == 1: #axisymmetric case 
                self.idlObj = idl.Generate_TMC_Initial_Data_Line(inp.geom, self.coneSol, inp.gasProps, nPts=inp.nIdlPts, endpoints=inp.idlEndPts)
            elif inp.delta == 0: #2D case
                self.idlObj = idl.Generate_2D_Initial_Data_Line(inp, self.rampSol, inp.gasProps, nPts=inp.nIdlPts, endpoints=inp.idlEndPts)

        elif inp.init_method == "MACH LINE":
            if inp.delta==1: #axisymmetric
                self.idlObj = idl.Generate_TMC_Initial_Data_Line(inp.geom, self.coneSol, inp.gasProps, nPts=inp.nIdlPts)
            elif inp.delta==0: #2D case 
                self.idlObj = idl.Generate_2D_Initial_Data_Line(inp, self.rampSol, inp.gasProps, nPts=inp.nIdlPts)
        else: 
            raise ValueError(f"Invalid Mesh Initialization Method Specified: {inp.init_method}")
        
        #GENERATING POINTS UPSTREAM OF IDL + CHAR 
        if inp.delta == 1:
            self.upstream_data = idl.Generate_TMC_Initial_Data_Line(inp.geom, self.coneSol, inp.gasProps, upstream_idl=True)
        elif inp.delta == 0: 
            self.upstream_data = idl.Generate_2D_Initial_Data_Line(inp, self.rampSol, inp.gasProps,upstream_idl=True)
        
        #GENERATING MESH 
        if inp.compute_shocks:        
            self.mesh = moc.Mesh(inp, idl=self.idlObj, explicit_shocks=True) #shocked mesh
        else:  
            self.mesh = moc.Mesh(inp, idl=self.idlObj) #shockless mesh 
        
        self.solution_runtime = time.perf_counter()-t0 #storing runtime

    def plot_solution(self, plotFile:str) -> None:
        """
        load in plot file and create plots from save file
        """
        import matplotlib.pyplot as plt
        import post_processing.post_process as post_process 
        import json 

        print("\ngenerating figures...\n")
        try: plotDict = json.load(open(plotFile, 'r'))
        except: plotDict = json.load(open("post_processing/"+plotFile, 'r'))
        plotSettings = plotDict["global plot settings"]
        del plotDict["global plot settings"]
        
        for key in plotDict.keys():
            subDict = plotDict[key] 
            subDict["figure name"] = key
            #hand off subDict to the post processing module
            #post_process.create_slice_plot(subDict, self, plotSettings)
            post_process.create_figure(subDict, self, plotSettings)

        plt.show() 

    def print_details(self) -> None:
        """
        prints relevant solution information to console
        to be called after run_solution() is called 
        run time
        total pressure loss
        number of mesh points
        taylor maccoll solution
        TODO
        """
        print("\nRESULTS:")
        print(f"\tRun Time: {round(self.solution_runtime,3)} seconds")
        print(f"\tNumber of Mesh Points: {len(self.mesh.meshPts)}")   
        print(f"\tNumber of Regions: {len(self.mesh.p0_ratio_by_region)}") 
        print(f"\tTotal Pressure Ratio By Region:\
               {[round(p,4) for p in self.mesh.p0_ratio_by_region]}")  

    def export_results(self) -> None:
        import os 
        import pandas as pd
        import numpy as np 
        import math
        print("\nexporting solution...")
        export_name = "save_" + self.inputs.geom.name + f"_M{self.inputs.M_inf}.csv"

        if os.path.isfile(export_name):
            abort = True
            inp = input(f"\t{export_name} already exists. Overwrite? [y/n]: ")
            if inp in ["Y", "y"]:
                abort = False 

            if abort: 
                print("\texport aborted\n")
                return 

        geomtype = ["AXI" if self.inputs.delta == 1 else "2D"][0]
        basic_info_list = [
            f"Geom Type:,           {geomtype}\n",
            f"M_inf:,               {self.inputs.M_inf}\n"
            f"gamma:,               {self.inputs.gam}\n"
            f"p0 region ratios:,    {[round(p, 4) for p in self.mesh.p0_ratio_by_region]}\n"
        ]

        if hasattr(self, "coneSol"):
            method = "Taylor-Maccoll Equation"
            def_ = math.degrees(self.coneSol.cone_ang)
            beta = math.degrees(self.coneSol.shock_ang)
            p02_p01 = self.coneSol.p02_p01

        elif hasattr(self, "rampSol"):
            method = "2D Oblique Shock"
            def_ = math.degrees(self.rampSol.deflec)
            beta = math.degrees(self.rampSol.beta)
            p02_p01 = self.rampSol.p02_p01

        init_solution_list = [
            f"solution method:,             {method}\n"
            f"initial deflection (deg):,    {def_}\n",
            f"wave angle (deg):,            {beta}\n"
            f"p02/p01:,                     {p02_p01}\n"
        ]

        x_cb = np.linspace(self.inputs.geom.centerbody_bounds[0], \
                           self.inputs.geom.centerbody_bounds[-1], 100)
        x_cowl = np.linspace(self.inputs.geom.cowl_bounds[0], \
                             self.inputs.geom.cowl_bounds[-1], 100)
        geometry_dict = {
            "centerbody x":         x_cb, 
            "centerbody y":         [self.inputs.geom.y_centerbody(x) for x in x_cb], 
            "centerbody dydx":      [self.inputs.geom.dydx_centerbody(x) for x in x_cb],
            "cowl x":               x_cowl, 
            "cowl y":               [self.inputs.geom.y_cowl(x) for x in x_cowl], 
            "cowl dydx":            [self.inputs.geom.dydx_cowl(x) for x in x_cowl] 
            
        } 
        df_geom = pd.DataFrame(geometry_dict)
        df_geom.insert(3, None, None)

        cb_mesh_points = {"index": [],"x": [],"y": [],"u": [],"v": [],"Mach": [],\
                           "Region": [], "p/p0f": [], "T/T0f": [], "rho/rho0f": []}
        for pt in self.mesh.wallPtsLower:
            cb_mesh_points["index"].append(pt.i)
            cb_mesh_points["x"].append(pt.x)
            cb_mesh_points["y"].append(pt.y)
            cb_mesh_points["u"].append(pt.u)
            cb_mesh_points["v"].append(pt.v)
            cb_mesh_points["Mach"].append(pt.mach)
            cb_mesh_points["Region"].append(pt.reg)
            cb_mesh_points["p/p0f"].append(pt.p_p0f)
            cb_mesh_points["T/T0f"].append(pt.T_T0)
            cb_mesh_points["rho/rho0f"].append(pt.rho_rho0f)
        df_cb_pts = pd.DataFrame(cb_mesh_points)

        cowl_mesh_points = {"index": [],"x": [],"y": [],"u": [],"v": [],"Mach": [],\
                           "Region": [], "p/p0f": [], "T/T0f": [], "rho/rho0f": []}
        for pt in self.mesh.wallPtsUpper:
            cowl_mesh_points["index"].append(pt.i)
            cowl_mesh_points["x"].append(pt.x)
            cowl_mesh_points["y"].append(pt.y)
            cowl_mesh_points["u"].append(pt.u)
            cowl_mesh_points["v"].append(pt.v)
            cowl_mesh_points["Mach"].append(pt.mach)
            cowl_mesh_points["Region"].append(pt.reg)
            cowl_mesh_points["p/p0f"].append(pt.p_p0f)
            cowl_mesh_points["T/T0f"].append(pt.T_T0)
            cowl_mesh_points["rho/rho0f"].append(pt.rho_rho0f)
        df_cowl_pts = pd.DataFrame(cowl_mesh_points)

        with open(export_name, 'w') as file: 
            file.write("BASIC INFO\n")
            [file.write(line) for line in basic_info_list]
            file.write("\nINITIALIZATION SOLUTION\n")
            [file.write(line) for line in init_solution_list]
            file.write("\nCENTERBODY & COWL GEOMETRY\n")
            df_geom.to_csv(file, index=False, lineterminator="\n")
            file.write("\nCENTERBODY MESH POINTS\n")
            df_cb_pts.to_csv(file, index=False, lineterminator="\n")
            file.write("\nCOWL MESH POINTS\n")
            df_cowl_pts.to_csv(file, index=False, lineterminator="\n")

        print(f"\texported: {export_name} to {os.getcwd()}\n")
