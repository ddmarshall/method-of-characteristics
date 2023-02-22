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
            #TODO add check for endpoints on the wall and set booleans to true
        self.currGen = self.meshPts.copy()

        self.funcs = moc_op.operator_funcs() #operator functions
        self.gasProps = gasProps
        self.delta = delta
        self.pcTOL = pcTOL
        self.geom = Geom

    def next_generation(self):
        #creates next generation of mesh points
        #!WORK IN PROGRESS...
        #Add clause for symmetry boundary 
        nextGen = []
        for i,pt in enumerate(self.currGen): 
            #interior solution:
            if i == 0 and pt.isWall is not True:
                #upper wall solution (direct)
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
                    intersecPts = self.check_curve_intersect(pt3, pt2, pt1)
                    self.triangle.append([pt3.i, pt2.i, pt1.i])
                    self.meshPts.append(pt3)
                    if intersecPts is not None: 
                        self.correct_curve_intersect(intersecPts) #apply correction to crossed characteristics
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

        print("\nkill function triggered\n")
        for meshPt in self.meshPts: 
            meshPt.get_point_properties(self.gasProps)

    def check_boundary_breach(self,x,y):
        #check if a point object has breached the boundary
        if y >= self.geom.y_cowl(x) or y <= self.geom.y_centerbody(x):
            return True
        return False

    def check_curve_intersect(self, pt3, pt2, pt1):
        #!Work-In-Progress
        #TODO This function suffers from egregious nesting (keep the tech debt flowing!)
        #checks if a new point creates a line which crosses the same family of characteristic
        #call when a new point is generated
        #pChar is of form 
        def ccw(A,B,C):
            return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])
        # Return true if line segments AB and CD intersect
        def intersect(A,B,C,D):
            #check if points intersect at the ends (guard clause)
            if A in [C,D] or B in [C,D]:
                return False
            return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

        for pt in self.currGen: 
            #search through triangle matrix to find relevant segments
            #recall triangle[i] = [node point i, upstream above, upstream below]
            posSeg_new = [[pt3.x, pt3.y],[pt1.x, pt1.y]] #proposed positive line segment 
            negSeg_new = [[pt3.x, pt3.y],[pt2.x, pt2.y]] #proposed negative line segment 

            for tri in [x for x in self.triangle if x[0] == pt.i]:
                if tri[2] is not None: 
                    posPt = [pt for pt in self.meshPts if pt.i == tri[2]][0] #should only be one point
                    posSeg = [[pt.x,pt.y],[posPt.x, posPt.y]]
                    #check for intersection
                    if intersect(posSeg_new[0], posSeg_new[1], posSeg[0], posSeg[1]): 
                        print(f"\nintersection of mesh segment [{pt3.i}-{pt1.i}]")
                        return [[pt, posPt],[pt3, pt1]]
                if tri[1] is not None:  
                    negPt = [pt for pt in self.meshPts if pt.i == tri[1]][0] #should only be one point 
                    negSeg = [[pt.x, pt.y],[negPt.x, negPt.y]]
                    if intersect(negSeg_new[0], negSeg_new[1], negSeg[0], negSeg[1]):
                        print(f"\nintersection of mesh segment [{pt3.i}-{pt2.i}]")
                        return [[pt, negPt],[pt3, pt2]]

        return None

    def correct_curve_intersect(self, points): 
        """
        applies corrective measures to deal with characteristic curve intersection
        see Z&H Vol2 Figure 16.19(c) for method
        TODO combine with check function (functions will always be called together if check intersect returns true)
        
        !Careful with organizing the points list 
        points = [[A,B],[C,D]]
        intersectionpoint <-> B is the segment to be deleted
        C<->D is segment to be kept and interpolated along 
        A<->B will be trimmed to A<->intersectionpoint
        """
        #unpacking
        D,C = points[1][0], points[1][1] #segment running to new point
        B,A = points[0][0], points[0][1] #segment running to previous gen point
            
        #find coordinates of intersection point:
        m_ab = (A.y - B.y)/(A.x - B.x)
        m_cd = (C.y - D.y)/(C.x - D.x)
        x = (m_ab*A.x - m_cd*C.x - A.y + C.y)/(m_ab - m_cd) 
        y = m_cd*(x - C.x) + C.y

        #linear interpolate to get velocity components at intersection point:
        l_cd = math.sqrt((C.x - D.x)**2 + (C.y - D.y)**2)
        l_cpt = math.sqrt((C.x - x)**2 + (C.y - y)**2)
        u = (D.u - C.u)*l_cpt/l_cd + C.u
        v = (D.v - C.v)*l_cpt/l_cd + C.v

        #create mesh point and triangle segments
        ind = B.i
        pt = mesh_point(x,y,u,v,ind)
        self.meshPts.remove(B)
        self.meshPts.insert(ind, pt)
        [self.triangle.remove(x) for x in self.triangle if x[0] == B.i] #delete segments connecting B to A and other
        #delete segment connecting C and D
        for i,tri in enumerate(self.triangle):
            if tri[0] == D.i: 
                self.triangle[i] = [D.i,pt.i,None]
                break
            if ind in tri[1:2]: 
                self.triangle.remove(tri)

        self.triangle.append([ind, C.i, A.i]) #adding segments connecting intersectpoint and A and C
        self.currGen.remove(B)
        self.currGen.insert(B.i,D)
class mesh_point: 
    def __init__(self,x,y,u,v,ind,isWall=False):
         self.x,self.y,self.u,self.v = x,y,u,v 
         self.i = ind
         self.isWall = isWall #is the point on the boundary? 

    def get_point_properties(self, gasProps): 
        #unpacking
        gam, a0, T0, p0 = gasProps.gam, gasProps.a0, gasProps.T0, gasProps.p0
        V = math.sqrt(self.u**2 + self.v**2)
        a = math.sqrt(a0**2 + 0.5*(gam-1)*V**2)
        self.mach = V/a
        self.T = T0/(1+0.5*(gam-1)*(V/a)**2)
        self.p = p0*(T0/self.T)**(gam/(gam-1))