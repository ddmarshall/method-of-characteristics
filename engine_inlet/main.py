"""
Strings everything together and runs an inlet solution
TODO: take in a run file (.toml or .json) 
"""
if __name__ == "__main__":
    #Importing Stuff
    import math
    import numpy as np
    import method_of_characteristics.moc_mesh_generator as moc
    import taylor_maccoll_cone.taylor_maccoll as tmc
    import initial_data_line.idl as idl
    import post_processing.post_process as post_process

    #Geometry 
    import example_geometry as geom
    inlet = geom.inletGeom()

    #Flow Conditions
    gam, M_inf, R, T0, p0 = 1.4, 2.5, 287.05, 288.15, 101325
    class gasProps:
        def __init__(self, gam, R, T0, p0): 
            self.gam, self.R, self.T0, self.p0 = gam, R, T0, p0
            self.a0 = math.sqrt(gam*R*T0)
    gas = gasProps(gam, R, T0, p0)

    #Run Taylor Maccoll Solution 
    cone = tmc.TaylorMaccoll_Cone(math.radians(inlet.cone_ang_deg), M_inf, gas) 

    #Generate initial data line 
    class make_curve:
        def __init__(self, y_x, dist, endpoints):
            self.y_x, self.dist, self.endpoints = y_x, dist, endpoints
    nPoints = 10
    dist = np.linspace(0,1+(1/nPoints),nPoints+1)
    curve =  make_curve(lambda x: 4*(x-2.5)**2, dist, (2.01,2.15))
    idlObj = idl.generate_tmc_initial_data_line(cone, curve, gas)
    
    #Generate Mesh
    delta=1 #axisymmetric flow
    pcTOL=0.0001 #percent change convergence tolerance
    masterMesh = moc.mesh(idlObj, inlet, gas, delta, pcTOL) #create mesh object
    masterMesh.generate_mesh(lambda masterMesh: masterMesh.numGens > 10) #generate mesh

    #Plot results
    plotObj = post_process.create_slice_plot(coneSol=cone, inletGeom=inlet, idl=idlObj, mesh=masterMesh)