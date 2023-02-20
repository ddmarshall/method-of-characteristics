"""
This class generates a universal input object from a set of user input files
"""
class inputObj: 
    def __init__(self, inpFile, geomFile):
        
        #check file name extension for filetype
        fileName = inpFile.split("/")[-1]
        ext = fileName.split(".")[-1]

        if ext == "json":
            self.import_json(inpFile) #if json is supplied
        elif ext == "toml":
            self.import_toml(inpFile) #if toml is supplied 
        else: 
            raise NameError(f'Invalid Filetype: .{ext}')

        self.import_geom(geomFile) #get geometry

    def import_json(self,file):

        import json
        json_data = json.load(open(file))
        json_translator = { #converts dictionary keys to object attributes 
            #Gas Properties
            "Freestream Stag Temp (K)":     "T0",
            "Spec Heat Ratio":              "gam",
            "Ideal Gas Constant (J/kgK)":   "R",
            #Flow Properties
            "Freestream Mach":              "M_inf",
            #Geometry 
            "file":                         "geomFile",
            #MOC Settings
            "Delta":                        "delta",
            "Unit Process Converge TOL":    "pcTOL",
            #Initial Data Line
            "function":                     "funcStr",
            "distribution":                 "dist",
            "endpoints":                    "endPts"
        }
        #Set attributes according to translator. Good luck comprehending this comprehension:      
        [setattr(self, json_translator[key_j], json_data[key_i][key_j]) for key_i in list(json_data.keys()) for key_j in json_data[key_i].keys() if key_j in json_translator.keys()]
        
    def import_toml(self, file):
        #TODO write this part if using .toml files for inputs
        import toml 

    def import_geom(self):
        """
        Accesses geometry file and stores in input object
        TODO write this
        """
        import self.geomFile
        geom = inletGeom()
        self.geom = geom