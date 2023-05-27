import json
"""
This class generates a universal input object from a set of user input files
"""
class Input: 
    
    def __init__(self, inpFile, geomFile):
        #check file name extension for filetype
        fileName = inpFile.split("/")[-1]
        ext = fileName.split(".")[-1]
        print(f"loading input file: {inpFile}")
        if ext == "json":
            self.import_json(inpFile) #if json is supplied
        elif ext == "toml":
            self.import_toml(inpFile) #if toml is supplied 
        else: 
            raise NameError(f'Invalid Filetype: .{ext}')

        self.import_geom(geomFile, preview_geom=True) #store geometry object as attribute
        self.print_inputs()

    def import_json(self,file):
        
        try: 
            json_data = json.load(open(file))
        except:
            json_data = json.load(open("Input/" + file))
        json_translator = { #converts dictionary keys to object attributes 
            #Gas Properties
            "Freestream Stag Temp (K)":     "T0",
            "Freestream Stag Pres (Pa)":    "p0",
            "Spec Heat Ratio":              "gam",
            "Ideal Gas Constant (J/kgK)":   "R",
            #Flow Properties
            "Freestream Mach":              "M_inf",
            #MOC Settings
            #"Delta":                        "delta",
            "Initiation":                   "init_method",
            "Compute Shocks":               "compute_shocks",
            "Convergence Tolerance":        "pcTOL",
            "Kill Function":                "kill",
            #Initial Data
            "Endpoints":                    "idlEndPts",
            "Num Points":                   "nIdlPts",
        }
        #Set attributes according to translator. Good luck comprehending this comprehension:      
        [setattr(self, json_translator[key_j], json_data[key_i][key_j]) for \
         key_i in list(json_data.keys()) for key_j in json_data[key_i].keys()\
              if key_j in json_translator.keys()]

    def import_geom(self, geomFile, preview_geom=True):
        import geometry.geom_pre_processor as gproc
        """
        Accesses geometry file and stores in input object
        """
        try: 
            geomDict = json.load(open(geomFile))
        except: 
            geomDict = json.load(open("geometry/" + geomFile))
            
        geomObj = gproc.Inlet_Geom(geomDict)
        self.geom = geomObj
        self.geom.x_cowl_lip = geomObj.cowl_bounds[0]
        if self.geom.geom_type in ["2D", "2d"]:
            self.delta = 0
        elif self.geom.geom_type in ["Axi", "AXI", "axi"]:
            self.delta = 1

        if preview_geom:
            #generate a preview of the loaded geometry before continuing 
            import matplotlib.pyplot as plt
            import numpy as np 
            x_range_cb = np.linspace(geomObj.centerbody_bounds[0], geomObj.centerbody_bounds[-1], 100)
            x_range_cowl = np.linspace(geomObj.cowl_bounds[0], geomObj.cowl_bounds[-1], 100)
            plt.figure(figsize=(10,4))
            plt.plot(x_range_cb, [geomObj.y_centerbody(x) for x in x_range_cb], label="centerbody")
            plt.plot(x_range_cowl, [geomObj.y_cowl(x) for x in x_range_cowl], label="cowl")
            plt.grid(), plt.legend(), plt.title("Inlet Geometry Preview")
            plt.show()

    def print_inputs(self):
        """
        prints important inputs to console
        """
        print("\nINPUTS:")
        print(f"\tMach Number: {self.M_inf}")
        print(f"\tSpecific Heat Ratio: {self.gam}")
        if self.delta == 0: 
            print("\tGeometry: 2-Dimensional")
        elif self.delta == 1:
            print("\tGeometry: Axisymmetric")
        print(f"\tMesh Initiation: {self.init_method}")
