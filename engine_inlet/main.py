"""
***MAIN RUN FILE***
Author: Shay Takei 
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
        class freeStream:
            def __init__(frst, inputObj):
                M, p0, T0, gam = inputObj.M_inf, inputObj.p0, inputObj.T0, inputObj.gam
                frst.mach = M
                frst.p = p0*(1 + 0.5*(gam - 1)*M**2)**(-gam/(gam-1))
                frst.T = T0/(1 + 0.5*(gam - 1)*M**2)

        self.freestream = freeStream(inpObj)

    def run_solution(self):
        print("\nrunning solution...\n")
        import method_of_characteristics.moc_mesh_engine as moc
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
        if self.inputs.delta == 1: #axisymmetric
            self.coneSol = tmc.TaylorMaccoll_Cone(math.radians(inp.geom.cone_ang_deg), inp.M_inf, gas) 

        #generate IDL
        self.idlObj = idl.generate_tmc_initial_data_line(self.inputs.geom, self.coneSol, gas, self.inputs.nIdlPts, self.inputs.idlEndPts)

        #generate mesh
        #mesh = moc.mesh(self.idlObj, inp.geom, gas, inp.delta, inp.pcTOL) #create mesh object
        #mesh.generate_mesh(eval(inp.kill)) #generate mesh

        mesh = moc.Mesh(self.idlObj, inp.geom, gas, inp.delta, inp.pcTOL, eval(inp.kill), explicit_shocks=True) #shocked mesh 
        #mesh = moc.Mesh(self.idlObj, inp.geom, gas, inp.delta, inp.pcTOL, eval(inp.kill)) #shockless mesh 

        self.mesh = mesh 

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
    #plotfile = "plot_profile_mesh_only.json"
    plotfile = "plot_profile_test.json"
    sol = main(inputFile='user_inputs.json', geomObj=inlet, plotFile=plotfile) #run solution then plot results