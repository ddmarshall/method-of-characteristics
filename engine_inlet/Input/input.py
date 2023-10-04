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

        self.import_geom(geomFile) #store geometry object as attribute
        self.print_inputs()

    def import_json(self,file):
        
        try: 
            json_data = json.load(open(file))
        except:
            json_data = json.load(open("Input/" + file))
        json_translator = { #converts dictionary keys to object attributes 
            #Gas Properties
            "Freestream Stag Temp (K)":     "T0",
            "Spec Heat Ratio":              "gam",
            "Ideal Gas Constant (J/kgK)":   "R",
            #Flow Properties
            "Freestream Mach":              "M_inf",
            #MOC Settings
            "Initiation":                   "init_method",
            "Compute Shocks":               "compute_shocks",
            "Convergence Tolerance":        "pcTOL",
            #Initial Data
            "Endpoints":                    "idlEndPts",
            "Num Points":                   "nIdlPts",
        }
        #Set attributes according to translator. Good luck comprehending this comprehension:      
        [setattr(self, json_translator[key_j], json_data[key_i][key_j]) for \
         key_i in list(json_data.keys()) for key_j in json_data[key_i].keys()\
              if key_j in json_translator.keys()]

    def import_geom(self, geomFile):
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
