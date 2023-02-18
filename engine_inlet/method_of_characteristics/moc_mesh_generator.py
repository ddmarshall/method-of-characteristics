#import os
#import sys
#sys.path.append(os.getcwd())
import method_of_characteristics.unit_processes as moc_op
import math
"""
class for generating mesh and mesh point objects
"""
class mesh: 
    def __init__(self, idl, Geom, gasProps, delta, velTOL):

        #create first generation from idl object
        
        self.meshPts = []
        self.triangle = []
        self.idlLen = len(idl.x) #number of points in idl
        self.numGens = 1
        for i,x in enumerate(idl.x): 
            self.meshPts.append(mesh_point(x, idl.y[i], idl.u[i], idl.v[i], i))
        self.currGen = self.meshPts.copy()

        self.funcs = moc_op.operator_functions() #operator functions
        self.gasProps = gasProps
        self.delta = delta
        self.velTOL = velTOL
        self.geom = Geom

    def next_generation(self):
        #creates next generation of mesh points
        #!WORK IN PROGRESS...
        #TODO add check to see if interior point will end up crossing a boundary 
        nextGen = []
        for i,pt in enumerate(self.currGen): 
            #interior solution:
            if i == 0 and pt.isWall is not True:
                #upper wall solution
                pt1 = pt
                x3, y3, u3, v3 = moc_op.direct_wall_abv(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.velTOL, self.funcs) 
                ind = len(self.meshPts)
                pt3 = mesh_point(x3, y3, u3, v3, ind,isWall=True)
                self.meshPts.append(pt3)
                self.triangle.append([pt3.i, None, pt1.i])
                nextGen.append(pt3)

            if i < len(self.currGen)-1: 
                #interior point
                pt2 = pt #y(pt2) > y(pt1)
                pt1 = self.currGen[i+1]
                x3, y3, u3, v3 = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.velTOL, self.funcs)
                if self.check_boundary_breach(x3,y3) == False: #check if new point breaches upper or lower boundary
                    ind = len(self.meshPts)
                    pt3 = mesh_point(x3, y3, u3, v3, ind)
                    self.meshPts.append(pt3)
                    self.triangle.append([pt3.i, pt2.i, pt1.i])
                    nextGen.append(pt3)

            elif i == len(self.currGen)-1 and pt.isWall is not True:
                #lower wall solution
                pt2 = pt
                x3, y3, u3, v3 = moc_op.direct_wall_bel(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.velTOL, self.funcs)
                ind = len(self.meshPts)
                pt3 = mesh_point(x3, y3, u3, v3, ind, isWall=True)
                self.meshPts.append(pt3)
                self.triangle.append([pt3.i, pt2.i, None])
                nextGen.append(pt3)     

        self.numGens += 1
        self.currGen = nextGen

    def generate_mesh(self, kill_func):
        #computes subsequent generations until kill_func(mesh object) evaluates as true
        while kill_func(self) != True: 
            self.next_generation() 

    def check_boundary_breach(self,x,y):
        #check if a point object has breached the boundary
        if y >= self.geom.y_cowl(x) or y <= self.geom.y_centerbody(x):
            return True
        return False

class mesh_point: 
    def __init__(self,x,y,u,v,ind,isWall=False):
         self.x,self.y,self.u,self.v = x,y,u,v 
         self.i = ind
         self.isWall = isWall #is the point on the boundary? 

    def get_temp_pressure(self, a0, T0, p0, gam): 
        #temperature 
        V = math.sqrt(self.u**2 + self.v**2)
        a = math.sqrt(a0**2 + 0.5*(gam-1)*V**2)
        self.T = T0/(1+0.5*(gam-1)*(V/a)**2)
        #pressure
        self.p = p0*(T0/self.T)**(gam/(gam-1))

"""
if __name__ == "__main__":
    #TESTING MODULE
    
    #LOAD IN INLET GEOMETRY 
    import example_geometry as geom
    inlet = geom.inletGeom()

    #CONE SOLUTION 
    import taylor_maccoll_cone.taylor_maccoll as tmc
    cone_ang = math.radians(12.5)
    M_inf = 2.5
    class gasProps:
        def __init__(self, gam, R, T0): 
            self.gam, self.R, self.T0 = gam, R, T0
    gas = gasProps(1.4, 287.05, 288.15)
    cone = tmc.TaylorMaccoll_Cone(cone_ang, M_inf, gas) 

    #IDL 
    import initial_data_line.idl as IDL
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

"""