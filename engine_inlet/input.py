"""
This class generates a universal input object from a set of user input files
"""
class inputObj: 
    
    def __init__(self, inpFile, geomObj):
        
        #check file name extension for filetype
        fileName = inpFile.split("/")[-1]
        ext = fileName.split(".")[-1]

        if ext == "json":
            self.import_json(inpFile) #if json is supplied
        elif ext == "toml":
            self.import_toml(inpFile) #if toml is supplied 
        else: 
            raise NameError(f'Invalid Filetype: .{ext}')

        self.import_geom(geomObj) #store geometry object as attribute

    def import_json(self,file):

        import json
        json_data = json.load(open(file))
        json_translator = { #converts dictionary keys to object attributes 
            #Gas Properties
            "Freestream Stag Temp (K)":     "T0",
            "Freestream Stag Pres (Pa)":    "p0",
            "Spec Heat Ratio":              "gam",
            "Ideal Gas Constant (J/kgK)":   "R",
            #Flow Properties
            "Freestream Mach":              "M_inf",
            #MOC Settings
            "Delta":                        "delta",
            "Unit Process Converge TOL":    "pcTOL",
            "kill function":                "kill",
            #Initial Data Line
            #"function":                     "idlFuncStr",
            #"distribution":                 "idlDist",
            "endpoints":                    "idlEndPts",
            "num points":                   "nIdlPts"
        }
        #Set attributes according to translator. Good luck comprehending this comprehension:      
        [setattr(self, json_translator[key_j], json_data[key_i][key_j]) for key_i in list(json_data.keys()) for key_j in json_data[key_i].keys() if key_j in json_translator.keys()]
        
    def import_toml(self, file):
        #TODO write this part if using .toml files for inputs
        import toml 

    def import_geom(self, geomObj):
        """
        Accesses geometry file and stores in input object
        TODO write this
        """
        self.geom = geomObj