"""
Strings everything together and runs an inlet solution
! saving to pickle is current not working
"""
class main:
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
        inpObj = inp.inputObj(inpFile, geomObj)
        self.inputs = inpObj

        #Kind of a dirty solution, but adds freestream object for plotting purposes
        #TODO consider rewriting
        class freeStream:
            def __init__(frst, inputObj):
                M, p0, T0, gam = inputObj.M_inf, inputObj.p0, inputObj.T0, inputObj.gam
                frst.mach = M
                frst.p = p0*(1 + 0.5*(gam - 1)*M**2)**(-gam/(gam-1))
                frst.T = T0/(1 + 0.5*(gam - 1)*M**2)

        self.freestream = freeStream(inpObj)

    def run_solution(self):
        print("\nrunning solution...\n")
        import method_of_characteristics.moc_mesh_generator as moc
        import taylor_maccoll_cone.taylor_maccoll as tmc
        import initial_data_line.idl as idl
        import post_processing.post_process as post_process
        import math
        import numpy as np

        inp = self.inputs
        
        class gasProps:
            def __init__(self, gam, R, T0, p0): 
                self.gam, self.R, self.T0, self.p0 = gam, R, T0, p0
                self.a0 = math.sqrt(gam*R*T0)
        gas = gasProps(inp.gam, inp.R, inp.T0, inp.p0) #create gas properties object 

        #run taylor maccoll
        self.coneSol = tmc.TaylorMaccoll_Cone(math.radians(inp.geom.cone_ang_deg), inp.M_inf, gas) 

        #generate idl 
        #class make_curve:
        #    def __init__(self, y_x, dist, endpoints):
        #        self.y_x, self.dist, self.endpoints = y_x, dist, endpoints
        
        #if inp.idlDist[0] == "linear": 
        #    Dist = np.linspace(0,1+(1/inp.idlDist[1]),inp.idlDist[1]+1)

        #curve = make_curve(eval(inp.idlFuncStr), Dist, inp.idlEndPts)
        self.idlObj = idl.generate_tmc_initial_data_line(self.inputs.geom, self.coneSol, gas, self.inputs.nIdlPts, self.inputs.idlEndPts)

        #generate mesh
        mesh = moc.mesh(self.idlObj, inp.geom, gas, inp.delta, inp.pcTOL) #create mesh object
        mesh.generate_mesh(eval(inp.kill)) #generate mesh
        del mesh.funcs
        self.mesh = mesh 

    def store_solution(self, saveFile):
        #calling this function will overwrite existing files
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
            figname = key #TODO do something with this?
            #hand off subDict to the post processing module
            post_process.create_slice_plot(subDict, plotSettings, self)
        plt.show() 

    def print_details(self):
        """
        prints all relevant solution information to console
        TODO
        """
        pass 

if __name__ == "__main__":
    import post_processing.post_process as post_process
    import example_geometry as geom
    inlet = geom.Geom()
    #sol = main(inputFile='user_inputs.json', geomObj=inlet) #run solution then plot results
    sol = main(inputFile='user_inputs.json', geomObj=inlet, plotFile="plot_profile_test.json") #run solution then plot results
    #post_process.create_slice_plot(coneSol=sol.coneSol, inletGeom=sol.inputs.geom, idl=sol.idlObj, mesh=sol.mesh)
pass        