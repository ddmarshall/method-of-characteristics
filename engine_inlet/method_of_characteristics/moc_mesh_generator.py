import method_of_characteristics.unit_processes as moc_op
import math
"""
class for generating mesh and mesh point objects
"""
class mesh: 
    def __init__(self, idl, Geom, gasProps, delta, pcTOL):

        #create first generation from idl object
        
        self.meshPts = []
        self.triangle = []
        self.idlLen = len(idl.x) #number of points in idl
        self.numGens = 0
        for i,x in enumerate(idl.x): 
            self.meshPts.append(mesh_point(x, idl.y[i], idl.u[i], idl.v[i], i))
        self.currGen = self.meshPts.copy()

        self.funcs = moc_op.operator_functions() #operator functions
        self.gasProps = gasProps
        self.delta = delta
        self.pcTOL = pcTOL
        self.geom = Geom

    def next_generation(self):
        #creates next generation of mesh points
        #!WORK IN PROGRESS...
        nextGen = []
        for i,pt in enumerate(self.currGen): 
            #interior solution:
            if i == 0 and pt.isWall is not True:
                #upper wall solution
                pt1 = pt
                x3, y3, u3, v3 = moc_op.direct_wall_abv(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs) 
                ind = len(self.meshPts)
                pt3 = mesh_point(x3, y3, u3, v3, ind,isWall=True)
                self.meshPts.append(pt3)
                self.triangle.append([pt3.i, None, pt1.i])
                nextGen.append(pt3)

            if i < len(self.currGen)-1: 
                #interior point
                pt2 = pt #y(pt2) > y(pt1)
                pt1 = self.currGen[i+1]
                x3, y3, u3, v3 = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.pcTOL, self.funcs)
                if self.check_boundary_breach(x3,y3) == False: #check if new point breaches upper or lower boundary
                    ind = len(self.meshPts)
                    pt3 = mesh_point(x3, y3, u3, v3, ind)
                    self.triangle.append([pt3.i, pt2.i, pt1.i])
                    self.meshPts.append(pt3)
                    nextGen.append(pt3)

            elif i == len(self.currGen)-1 and pt.isWall is not True:
                #lower wall solution
                pt2 = pt
                x3, y3, u3, v3 = moc_op.direct_wall_bel(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs)
                ind = len(self.meshPts)
                pt3 = mesh_point(x3, y3, u3, v3, ind, isWall=True)
                self.meshPts.append(pt3)
                self.triangle.append([pt3.i, pt2.i, None])
                nextGen.append(pt3)     

        self.numGens += 1
        print(f"{len(nextGen)} points added in generation: {self.numGens}")
        self.currGen = nextGen

    def generate_mesh(self, kill_func):
        #continuously computes subsequent generations of characteristic mesh until kill_func(mesh object) evaluates as true
        while kill_func(self) != True: 
            self.next_generation() 

        for meshPt in self.meshPts: 
            meshPt.get_temp_pressure(self.gasProps)

    def check_boundary_breach(self,x,y):
        #check if a point object has breached the boundary
        if y >= self.geom.y_cowl(x) or y <= self.geom.y_centerbody(x):
            return True
        return False

    def check_curve_intersect(self):
        #!Work-In-Progress
        #checks if a new point creates a line which crosses the same family of characteristic
        #call when a new point is generated
        #pChar is of form 
        def ccw(A,B,C):
            return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0])
        # Return true if line segments AB and CD intersect
        def intersect(A,B,C,D):
            return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

class mesh_point: 
    def __init__(self,x,y,u,v,ind,isWall=False):
         self.x,self.y,self.u,self.v = x,y,u,v 
         self.i = ind
         self.isWall = isWall #is the point on the boundary? 

    def get_temp_pressure(self, gasProps): 
        #unpacking
        gam, a0, T0, p0 = gasProps.gam, gasProps.a0, gasProps.T0, gasProps.p0
        #temperature 
        V = math.sqrt(self.u**2 + self.v**2)
        a = math.sqrt(a0**2 + 0.5*(gam-1)*V**2)
        self.T = T0/(1+0.5*(gam-1)*(V/a)**2)
        #pressure
        self.p = p0*(T0/self.T)**(gam/(gam-1))