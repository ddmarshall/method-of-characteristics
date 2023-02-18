"""
This class generates an input object from a user defined input file
"""
class inputObj: 
    def __init__(self, file):
        fileName = file.split("/")[-1]
        ext = fileName.split(".")[-1]

        if ext == "json":
            self.import_json(file)
        elif ext == "toml":
            self.import_toml(file)
        else: 
            raise NameError(f'Invalid Filetype: .{ext}')

    def import_json(self,file):

        import json
        json_data = json.load(open(file))
        json_translator = {
            "Freestream Stag Temp (K)":     "T0",
            "Spec Heat Ratio":              "gam",
            "Ideal Gas Constant (J/kgK)":   "R",
            "Freestream Mach":              "M_inf",
            "Delta":                        "delta",
            "Velocity Tol (m/s)":           "vel_TOL"
        }
        #Set attributes according to translator. Good luck comprehending this comprehension:      
        [setattr(self, json_translator[key_j], json_data[key_i][key_j]) for key_i in list(json_data.keys()) for key_j in json_data[key_i].keys() if key_j in json_translator.keys()]
        
    def import_toml(self, file):
        import toml 
        #TODO write this part if using .toml files for inputs    