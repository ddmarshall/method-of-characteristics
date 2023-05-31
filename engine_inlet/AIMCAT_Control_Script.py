import AIMCAT as AIMCAT
"""
Control script for AIMCAT module. Specify input file locations and run AIMCAT by
 calling class "Main" with the appropriate inputs filled out. 

AIMCAT.Main(inputFile:str, geomFile:str, plotFile:str, export:bool=False, 
    preview_geom:bool=False)
    
    Inputs: 
        inputFile:      file path of main input .json file (relative to cwd)
        geomFile:       file path of geometry .json file (relative to cwd)
        plotFile:       file path of plot .json file (relative to cwd)
        export:         export results to csv when solution finished 
        preview_geom:   preview of loaded geometry before running 
    Returns: 
        None
"""
#SPECIFY INLET GEOMETRY FILE####################################################
inletFile = "single_cone_12_5deg.json"      
#inletFile = "2D_isentropic_ramp_5deg.json"
#inletFile = "NASA_D6078_Inlet.json"

#SPECIFY PLOTTING FILE##########################################################
plotfile = "plot_settings_test.json"
#plotfile = "plot_mesh.json"

#SPECIFY USER INPUT FILE########################################################
#inputFile = 'test_idl_straight_inputs.json'
inputFile = 'test_mach_line_idl_straight_inputs.json'

#RUN SOLUTION###################################################################
AIMCAT.Main(inputFile=inputFile, geomFile=inletFile, plotFile=None, \
            export=True, preview_geom=False)