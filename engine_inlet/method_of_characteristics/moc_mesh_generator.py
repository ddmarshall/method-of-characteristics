import unit_processes as moc_op
import math
"""
class for generating mesh and mesh point objects
"""
class mesh: 
    def __init__(self, idl, Geom, gasProps, delta, velTOL):

        #create first generation from idl object
        firstGen = []
        for i,x in enumerate(idl.x): 
            firstGen.append(mesh_point(x, idl.y[i], idl.u[i], idl.v[i]))
        self.gens = [firstGen] #create generation list
        
        self.funcs = moc_op.operator_functions() #operator functions
        self.gasProps = gasProps
        self.delta = delta
        self.velTOL = velTOL
        self.geom = Geom

    def next_generation(self):
        #creates next generation of mesh points
        #!WORK IN PROGRESS...
        newGen = []
        for i,pt2 in enumerate(self.gens[-1]): #iterate thru mesh points in most recent generation
            
            #run interior point solution for interior pairs
            pt1 = self.gens[-1][i-1]
            x3, y3, u3, v3 = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.velTOL, self.funcs)
            newGen.append(mesh_point(x3,y3,u3,v3))

            #run wall solution for edges
            if i == 0 or i == len(self.gens[-1])-1:
                #run wall or direct wall 
                pass
        
        self.gens.append(newGen) #add next generation of points to

    def generate_mesh(self, kill_func):
        #computes subsequent generations until kill_func(mesh object) evaluates as true
        while kill_func(self) != True: 
            self.next_generation()
        pass 

class mesh_point: 
    def __init__(self,x,y,u,v):
         self.x,self.y,self.u,self.v = x,y,u,v 

    def get_point_flow_properties(self):
        #calculates flow properties at piont
        #Mach, Static/Stagnation ratio, Entropy, etc. 
        #TODO
        pass 

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
    sys.path.append(os.getcwd()+"\\initial_data_line")
    import idl as IDL
    class make_curve:
        def __init__(self, y_x, dist, endpoints):
            self.y_x, self.dist, self.endpoints = y_x, dist, endpoints
    dist = [0, 0.2, 0.4, 0.6, 0.8, 1]
    curve =  make_curve(lambda x: 4*(x-2.5)**2, dist, (2.01,2.15))
    idlObj = IDL.generate_tmc_initial_data_line(cone, curve)

    #OPERATOR INPUTS
    class gasProps:
        def __init__(self, gam, R, T0): 
            self.gam, self.R, self.T0 = gam, R, T0
    gas = gasProps(1.4, 287.05, 288.15)
    delta=1 #axisymmetric flow
    velTOL = 0.0001 #velocity delta 

    #GENERATING
    masterMesh = mesh(idlObj, inlet, gas, delta, velTOL)
    while True:
        masterMesh.next_generation()
        print(masterMesh.gens[-1])
        input() #press return to move to next generation

    pass