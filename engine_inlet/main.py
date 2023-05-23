"""
***MAIN RUN FILE***
Author: Shay Takei 
"""
class Main:
    """
    This class controls all modules of AIMCAT.
    """
    def __init__(self, inputFile=None, geomFile=None, plotFile=None):
        print("="*160)
        ans = None
       
        if inputFile is not None and geomFile is not None: 
            self.load_inputs(inputFile, geomFile) #generate input object 
            self.run_solution() #run solution 
            self.print_details() #prints important details to the console

        if plotFile is not None: #if plotfile is provided, run it 
            self.plot_solution(plotFile)

    def load_inputs(self, inpFile, geomFile):
        import input as inp
        import math 
        inpObj = inp.inputObj(inpFile, geomFile)
        self.inputs = inpObj

        #Kind of a dirty solution, but adds freestream object for plotting purposes
        class freeStream:
            def __init__(frst, inputObj):
                M, p0, T0, gam, R = inputObj.M_inf, inputObj.p0, inputObj.T0, inputObj.gam, inputObj.R
                frst.mach = M
                frst.p = p0*(1 + 0.5*(gam - 1)*M**2)**(-gam/(gam-1))
                frst.p_p0 = frst.p/p0
                frst.T = T0/(1 + 0.5*(gam - 1)*M**2)
                frst.T_T0 = frst.T/T0
                frst.rho = frst.p/(R*frst.T)
                frst.rho_rho0 = frst.rho/(p0/(R*T0))
                a = math.sqrt(gam*R*frst.T)
                frst.u = M*a
                frst.v = 0 #!change if want angled flow for 2d case  

        self.inputs.freeStream = freeStream(inpObj)

    def run_solution(self):
        print("\nrunning solution...")
        import method_of_characteristics.moc_mesh_engine as moc
        import initial_data_line.idl as idl
        import math
        import time

        t0 = time.perf_counter() #intial time
        inp = self.inputs
        class gasProps:
            def __init__(self, gam, R, T0, p0): 
                self.gam, self.R, self.T0, self.p0 = gam, R, T0, p0
                self.a0 = math.sqrt(gam*R*T0)
        inp.gasProps = gasProps(inp.gam, inp.R, inp.T0, inp.p0) #create gas properties object 

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
            self.mesh = moc.Mesh(inp, eval(inp.kill), idl=self.idlObj, explicit_shocks=True) #shocked mesh
        else:  
            self.mesh = moc.Mesh(inp, eval(inp.kill), idl=self.idlObj) #shockless mesh 
        
        self.solution_runtime = time.perf_counter()-t0 #storing runtime

    def plot_solution(self, plotFile):
        """
        load in plot file and create plots from save file
        """
        import matplotlib.pyplot as plt
        import post_processing.post_process as post_process 
        import json 

        print("\ngenerating figures...\n")
        plotDict = json.load(open(plotFile, 'r'))
        plotSettings = plotDict["global plot settings"]
        del plotDict["global plot settings"]
        
        for key in plotDict.keys():
            subDict = plotDict[key] 
            subDict["figure name"] = key
            #hand off subDict to the post processing module
            #post_process.create_slice_plot(subDict, self, plotSettings)
            post_process.create_figure(subDict, self, plotSettings)

        plt.show() 


    def print_details(self):
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
        print(f"\tNumber of Regions: {len(self.mesh.tot_press_by_region)}") 
        print(f"\tTotal Pressure Ratio By Region: {[round(p/self.inputs.p0, 4) for p in self.mesh.tot_press_by_region]}")  

if __name__ == "__main__":

    inletFile = "geometry/single_cone_12_5deg.json"
    #inletFile = "geometry/2D_isentropic_ramp_5deg.json"
    #plotfile = "plot_profile_mesh_only.json"
    #plotfile = "plot_profile_test.json"
    plotfile = "plot_settings_test.json"
    #inputFile = 'test_idl_straight_inputs.json'
    inputFile = 'test_mach_line_idl_straight_inputs.json'
    sol = Main(inputFile=inputFile, geomFile=inletFile, plotFile=plotfile) #run solution then plot results