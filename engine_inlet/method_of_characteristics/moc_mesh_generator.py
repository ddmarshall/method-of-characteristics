import unit_processes as moc_op
import math
"""
class for generating mesh and mesh point objects
"""
class mesh: 
    def __init__(self, idl, Geom):

        #create first generation from idl object
        firstGen = []
        for i,x in enumerate(idl.x): 
            firstGen.append(mesh_point(x, idl.y[i], idl.u[i], idl.v[i]))
        self.gens = [firstGen] #create generation list
        
        funcs = moc_op.operator_functions() #operator functions

    def next_generation(self, Geom, funcs):
        #creates next generation of mesh points
        #!WORK IN PROGRESS...
        newGen = []
        for i,pt in enumerate(self.gens[-1]): #iterate thru mesh points in latest generation
            
            #run interior point solution for interior pairs
            pt1 = self.gens[-1][i]
            pt2 = self.gens[-1][i+1]
            x3, y3, u3, v3 = moc_op.interior_point(pt1, pt2, gasProps, delta, velTOL, funcs)
            newGen.append(mesh_point(x3,y3,u3,v3))
        
            #run wall solution for edges
            if i == 0 or i == len(self.gens[-1]-1):
                #run wall or direct wall 
                pass
        
        self.gens.append(newGen) #add next generation of points to

class mesh_point: 
    def __init__(self,x,y,u,v):
         self.x,self.y,self.u,self.v = x,y,u,v 

if __name__ == "__main__":
    #TESTING MODULE
    import os 
    import sys
    sys.path.append(os.getcwd())
    #LOAD IN INLET GEOMETRY 
    import example_geometry as geom
    inlet = geom.inletGeom()

    #CONE SOLUTION 
    sys.path.append(os.getcwd()+ "\\taylor_maccoll_cone")
    import taylor_maccoll as tmc
    gam = 1.4
    cone_ang = math.radians(12.5)
    M_inf = 2.5
    R = 287.05
    T0 = 288.15
    cone = tmc.TaylorMaccoll_Cone(cone_ang, M_inf, gam, R, T0) 

    #IDL 
    import idl as IDL
    class make_curve:
        def __init__(self, y_x, dist, endpoints):
            self.y_x, self.dist, self.endpoints = y_x, dist, endpoints
    dist = [0, 0.2, 0.4, 0.6, 0.8, 1]
    curve =  make_curve(lambda x: 4*(x-2.5)**2, dist, (2.01,2.15))
    idlObj = IDL.generate_tmc_initial_data_line(cone, curve)

    masterMesh = mesh(idlObj, inlet)
    pass