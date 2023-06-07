import aimcat as aimcat
"""
Control script for aimcat module. Specify input file locations and run AIMCAT by
calling class "Main" with the appropriate inputs filled out. 

aimcat.Main(inputFile:str, geomFile:str, plotFile:str or None, export:bool, 
    preview_geom:bool)
    
    Inputs: 
        inputFile:      file name of main input file
        geomFile:       file name of geometry file 
        plotFile:       file name of plot file 
        export:         export results to csv when solution finished 
        preview_geom:   preview of loaded geometry before running 
    Returns: 
        None
"""
#SPECIFY INLET GEOMETRY FILE####################################################
inletFile = "single_cone_12_5deg.json"      
#inletFile = "2D_isentropic_ramp_5deg.json"
#inletFile = "NASA_D6078_Inlet_Interpolated.json"
#inletFile = "NASA_D6078_Inlet_least_squares.json"

#SPECIFY PLOTTING FILE##########################################################
#plotfile = "plot_settings_test.json"
plotfile = "plot_mesh.json"

#SPECIFY USER INPUT FILE########################################################
#inputFile = 'test_idl_straight_inputs.json'
inputFile = 'test_mach_line_idl_straight_inputs.json'

#RUN SOLUTION###################################################################
aimcat.Main(inputFile=inputFile, geomFile=inletFile, plotFile=plotfile, \
            export=False, preview_geom=False)