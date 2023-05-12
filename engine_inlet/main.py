"""
***MAIN RUN FILE***
Author: Shay Takei 
"""
class Main:
    
    def __init__(self, inputFile=None, geomObj=None, saveFile=None, plotFile=None):
        ans = None
       
        if inputFile is not None and geomObj is not None: 
            self.load_inputs(inputFile, geomObj) #generate input object 

            if saveFile is not None:#if save file provided
                from os import path 
                if path.exists(saveFile): #check if save file is going to overwritten and warn user 
                    ans = input(f"\nWarning!: {saveFile} already exists. Proceed anyways? [y/n]: ")
                    if ans in ["n,N"]: return

            self.run_solution() #run solution 

            if saveFile is not None: 
                self.store_solution(saveFile) #save file

        if saveFile is not None and inputFile is None: #if only a save file is provided, load it
            self.load_solution(saveFile)

        if plotFile is not None: 
            self.plot_solution(plotFile)

    def load_inputs(self, inpFile, geomObj):
        print(f"\nloading input file: {inpFile}")
        import input as inp
        import math 
        inpObj = inp.inputObj(inpFile, geomObj)
        self.inputs = inpObj

        #Kind of a dirty solution, but adds freestream object for plotting purposes
        class freeStream:
            def __init__(frst, inputObj):
                M, p0, T0, gam, R = inputObj.M_inf, inputObj.p0, inputObj.T0, inputObj.gam, inputObj.R
                frst.mach = M
                frst.p = p0*(1 + 0.5*(gam - 1)*M**2)**(-gam/(gam-1))
                frst.T = T0/(1 + 0.5*(gam - 1)*M**2)
                a = math.sqrt(gam*R*frst.T)
                frst.u = M*a
                frst.v = 0 #!change if want angled flow for 2d case  

        self.inputs.freeStream = freeStream(inpObj)

    def run_solution(self):
        print("\nrunning solution...\n")
        import method_of_characteristics.moc_mesh_engine as moc
        import initial_data_line.idl as idl
        import math

        inp = self.inputs
        class gasProps:
            def __init__(self, gam, R, T0, p0): 
                self.gam, self.R, self.T0, self.p0 = gam, R, T0, p0
                self.a0 = math.sqrt(gam*R*T0)
        inp.gasProps = gasProps(inp.gam, inp.R, inp.T0, inp.p0) #create gas properties object 

        #RUNNING TAYLOR-MACCOLL or 2D RAMP SOLUTION: 
        if inp.delta==1:
            import taylor_maccoll_cone.taylor_maccoll as tmc 
            self.coneSol = tmc.TaylorMaccoll_Cone(math.radians(inp.geom.cone_ang_deg), inp.M_inf, inp.gasProps)
        elif inp.delta==0:
            import method_of_characteristics.oblique_shock as shock 
            deflec = math.radians(inp.geom.cone_ang_deg)
            self.rampSol = shock.Oblique_Shock(inp.M_inf, inp.gasProps.gam, inp.gasProps.R, deflec=deflec)

        #GENERATING INITIAL DATA LINE
        if inp.init_method == "STRAIGHT IDL":
            if inp.delta == 1: #axisymmetric case 
                self.idlObj = idl.Generate_TMC_Initial_Data_Line(inp.geom, self.coneSol, inp.gasProps, inp.nIdlPts, inp.idlEndPts)
            elif inp.delta == 0: #2D case
                self.idlObj = idl.Generate_2D_Initial_Data_Line(inp, self.rampSol, inp.gasProps, inp.nIdlPts, inp.idlEndPts)

        elif inp.init_method == "MACH LINE":
            if inp.delta==1: #axisymmetric
                self.idlObj = idl.Generate_TMC_Initial_Data_Line(inp.geom, self.coneSol, inp.gasProps, inp.nIdlPts)
            elif inp.delta==0: #2D case 
                self.idlObj = idl.Generate_2D_Initial_Data_Line(inp.geom, self.coneSol, inp.gasProps, inp.nIdlPts)
        
        else: 
            raise ValueError(f"Invalid Mesh Initialization Method Specified: {inp.init_method}")
        
        #GENERATING MESH        
        self.mesh = moc.Mesh(inp, eval(inp.kill), idl=self.idlObj, explicit_shocks=True) #shocked mesh 
        #self.mesh = moc.Mesh(inp, eval(inp.kill), idl=self.idlObj) #shockless mesh 

    def store_solution(self, saveFile):
        #calling this function will overwrite existing files
        #! Currently broken (pickling doesn't work will stored functions I think...)
        print(f"\npickling solution results to {saveFile}")
        import pickle 
        file = open(saveFile, 'ab')
        pickle.dump(self.mesh, file)
        file.close()

    def load_solution(self, saveFile):
        print(f"\nloading solution file: {saveFile}")
        import pickle
        file = open(saveFile, 'rb')
        try: 
            res = pickle.load(file)
            self.mesh = res.mesh 
            self.idlObj = res.idlObj
            self.inputs = res.inputs
            self.coneSol = res.coneSol
        except: pass 
        file.close()

    def plot_solution(self, plotFile):
        """
        load in plot file and create plots from save file
        """
        import matplotlib.pyplot as plt
        import post_processing.post_process as post_process 
        import json 

        plt.style.use('dark_background') #!temporary location
        plotDict = json.load(open(plotFile, 'r'))
        plotSettings = plotDict["default plot settings"]
        del plotDict["default plot settings"]
        
        for key in plotDict.keys():
            subDict = plotDict[key] 
            #hand off subDict to the post processing module
            post_process.create_slice_plot(subDict, self, plotSettings)

        plt.show() 

    def print_details(self):
        """
        prints all relevant solution information to console
        TODO
        """
        pass 


if __name__ == "__main__":

    import example_geometry as geom
    inlet = geom.Geom()
    plotfile = "plot_profile_mesh_only.json"
    #plotfile = "plot_profile_test.json"
    #inputFile = 'test_idl_straight_inputs.json'
    inputFile = 'test_mach_line_idl_straight_inputs.json'
    sol = Main(inputFile=inputFile, geomObj=inlet, plotFile=plotfile) #run solution then plot results