"""
Strings everything together and runs an inlet solution
TODO: take in a run file with everything known ahead of time 
"""

if __name__ == "__main__":
    #Importing Stuff
    import sys
    import os
    import math
    sys.path.append(os.getcwd() + "\\method_of_characteristics")
    sys.path.append(os.getcwd() + "\\taylor_maccoll_cone")
    sys.path.append(os.getcwd() + "\\initial_data_line")
    sys.path.append(os.getcwd() + "\\post_processing")
    import moc_mesh_generator as moc
    import taylor_maccoll as tmc
    import idl 
    import post_process

    #Geometry 
    import example_geometry as geom
    inlet = geom.inletGeom()

    #Flow Conditions
    gam = 1.4
    M_inf = 2.5
    R = 287.05
    T0 = 288.15
    class gasProps:
        def __init__(self, gam, R, T0): 
            self.gam, self.R, self.T0 = gam, R, T0
    gas = gasProps(gam, R, T0)

    #Run Taylor Maccoll Solution 
    cone = tmc.TaylorMaccoll_Cone(math.radians(inlet.cone_ang_deg), M_inf, gam, R, T0) 

    #Generate initial data line 
    class make_curve:
        def __init__(self, y_x, dist, endpoints):
            self.y_x, self.dist, self.endpoints = y_x, dist, endpoints
    dist = [0, 0.2, 0.4, 0.6, 0.8, 1]
    curve =  make_curve(lambda x: 4*(x-2.5)**2, dist, (2.01,2.15))
    idlObj = idl.generate_tmc_initial_data_line(cone, curve)
    
    #Generate Mesh
    delta=1 #axisymmetric flow
    velTOL = 0.0001 #velocity delta
    masterMesh = moc.mesh(idlObj, inlet, gas, delta, velTOL) #create mesh object
    masterMesh.generate_mesh(lambda masterMesh: len(masterMesh.gens) > 10) #generate mesh

    #Plot results
    plotObj = post_process.create_slice_plot(coneSol=cone, inletGeom=inlet, idl=idlObj, mesh=masterMesh)

    pass