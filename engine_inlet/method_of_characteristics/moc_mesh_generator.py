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
        if hasattr(idl, 'cowlPoint'):
            self.meshPts[0].isWall = True 
        if hasattr(idl, 'cbPoint'):
            self.meshPts[-1].isWall = True
        
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
        if self.geom.y_cowl(x) is None or self.geom.y_centerbody(x) is None: 
            return True

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

class mesh2:
    def __init__(self, idl, Geom, gasProps, delta, pcTOL, kill_func):
        self.C_pos, self.C_neg = [],[] #containers for characteristics lines
        self.triangle_obj = [] #container for point line segments
        self.idl = []
        self.funcs = moc_op.operator_funcs() #operator functions
        self.gasProps = gasProps
        self.delta = delta
        self.pcTOL = pcTOL
        self.geom = Geom
        self.f_kill = [kill_func, False]

        for i,x in enumerate(idl.x):
            self.idl.append(mesh_point(x, idl.y[i], idl.u[i], idl.v[i], None, isIdl=True))
            #TODO add check for endpoints on the wall and set booleans to true
        if hasattr(idl, 'cowlPoint'):
            self.idl[0].isWall = True 
        if hasattr(idl, 'cbPoint'):
            self.idl[-1].isWall = True
        
        self.generate_mesh()
        self.compile_mesh_points()
        [pt.get_point_properties(self.gasProps) for pt in self.meshPts]

    def generate_mesh(self):
        """
        generates the characteristic mesh until the kill function is triggered 
        """
        charDir = "neg" #!hard coded for now... starting direction
        self.generate_initial_mesh(charDir)
        while self.f_kill[1] == False: 

            if charDir == "neg":
                #generate initial above-wall point 
                pt1 = self.C_neg[-1][1] #2nd point in most recent characteristic
                [x3,y3,u3,v3] = moc_op.direct_wall_abv(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs)
                pt3 = mesh_point(x3,y3,u3,v3, None, isWall=True)
                self.triangle_obj.append([pt3, None, pt1]) 
                #self.C_neg.append([pt3])
                i,_ = self.find_mesh_point(pt1, self.C_pos)
                self.C_pos[i].append(pt3) #add to positive
                prev_n_char = self.C_neg[-1][1:]
                self.compute_next_neg_char(pt3, prev_n_char)

            elif charDir == "pos":

                #generate initial below-wall point 
                pt2 = self.C_pos[-1][1] #2nd point in most recent characteristic
                [x3,y3,u3,v3] = moc_op.direct_wall_bel(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs)
                pt3 = mesh_point(x3,y3,u3,v3, None, isWall=True)
                self.triangle_obj.append([pt3, pt2, None]) 
                #self.C_pos.append([pt3])
                i,_ = self.find_mesh_point(pt2, self.C_neg)
                self.C_neg[i].append(pt3) #add to positive
                
                prev_p_char = self.C_pos[-1][1:]
                self.compute_next_pos_char(pt3, prev_p_char)
            
        pass
 
    def generate_initial_mesh(self, charDir):
        """
        computes the initial mesh from the idl. For a vertical IDL this should form a triangle with either a leading + or - characterstic spanning wall to wall 
        """
        idl = self.idl 
        if charDir == "pos": 
            self.C_pos.append([idl[0]])
            [self.C_neg.append([pt]) for pt in idl]
            self.C_neg.reverse() #first neg line at bottom 
            
        elif charDir == "neg":
            self.C_neg.append([idl[-1]])
            [self.C_pos.append([pt]) for pt in idl] #first pos line at the top
            idl.reverse() #reverse idl to start at bottom

        #generating intial mesh: 
        for i,pt in enumerate(idl):
            if i == 0: continue #first point has no - char passing through it 
            if charDir == "neg":
                self.compute_next_neg_char(pt, self.C_neg[i-1])
            elif charDir == "pos":
                self.compute_next_pos_char(pt, self.C_pos[i-1])

    def compute_next_neg_char(self, init_point, prev_n_char, continueChar = False):
        """
        Generates the next leading negative characteristic by advancing the mesh along the previous
        init_point could be wall or idl point
        !NOTE TO SELF: CHANGES HERE NEED TO BE REFLECTED IN TWIN FUNCTION
        """
             
        if continueChar is False:
            self.C_neg.append([init_point])
        for pt in prev_n_char:
            pt2 = init_point #above
            pt1 = pt #below 
            [x3, y3, u3, v3] = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.pcTOL, self.funcs)
            pt3 = mesh_point(x3, y3, u3, v3, ind = None)
            if self.check_for_int_intersect(pt3, pt2, pt1, "neg"):
                self.trim_mesh_after_intersect(pt2, pt1, "neg")
                init_point = self.C_neg[-1][-1]
                i,j = self.find_mesh_point(init_point, self.C_pos)
                pt0 = self.C_pos[i][j-2]
                i,j = self.find_mesh_point(pt0, self.C_neg)
                prev_n_char = self.C_neg[i][j+1:]
                self.compute_next_neg_char(init_point, prev_n_char, continueChar=True) #recursion ooooh spooky 
                return 

            self.C_neg[-1].append(pt3)
            self.triangle_obj.append([pt3, pt2, pt1])

            #find point 1 and append the new point 
            i,_ = self.find_mesh_point(pt1, self.C_pos)
            self.C_pos[i].append(pt3)

            init_point = pt3 #update initial point

            if self.f_kill[0](self) == True:
                self.f_kill[1] = True
                return

        #terminating wall point
        pt2 = self.C_neg[-1][-1]
        [x3,y3,u3,v3] = moc_op.direct_wall_bel(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs)
        pt3 = mesh_point(x3, y3, u3, v3, None, isWall=True)
        self.C_neg[-1].append(pt3)
        self.triangle_obj.append([pt3, pt2, None])
        self.C_pos.append([pt3])

    def compute_next_pos_char(self, init_point, prev_p_char):
        """
        Generates the next leading positive characteristic by advancing the mesh along the previous
        init_point could be wall or idl point
        !OUTDATED !UPDATE THIS
        """
        self.C_pos.append([init_point])
        for pt in prev_p_char:
            pt1 = init_point #below #!
            pt2 = pt #above #!
            [x3, y3, u3, v3] = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.pcTOL, self.funcs)
            pt3 = mesh_point(x3, y3, u3, v3, ind = None)
            if self.check_for_int_intersect(pt3, pt2, pt1, "pos"):
                self.trim_mesh_after_intersect(pt2, pt1, "pos")
                return
            self.C_pos[-1].append(pt3)
            self.triangle_obj.append([pt3, pt2, pt1])

            #find point 2 and append the new point 
            i,_ = self.find_mesh_point(pt2, self.C_neg)
            self.C_neg[i].append(pt3)

            init_point = pt3

            if self.f_kill[0](self) == True:
                self.f_kill[1] = True
                return

        #terminating wall point
        pt1 = self.C_pos[-1][-1]
        [x3,y3,u3,v3] = moc_op.direct_wall_abv(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs)
        pt3 = mesh_point(x3, y3, u3, v3, None, isWall=True)
        self.C_pos[-1].append(pt3)
        self.triangle_obj.append([pt3, None, pt1])
        self.C_neg.append([pt3])

    def check_for_int_intersect(self, pt3, pt2, pt1, charDir):
        """
        checks for a same-family characteristic intersection for an interior point solution
        pt3 = new point from interior solution 
        pt2 = above parent point 
        pt1 = below parent point
        """
        hasIntersected = False 
        def ccw(A,B,C):
            return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])
        # Return true if line segments AB and CD intersect
        def intersect(A,B,C,D):
            #return true is segments A-B and C-D intersect 
            #check if points intersect at the ends (guard clause)
            if A in [C,D] or B in [C,D]:
                return False
            return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D) 

        if charDir == "neg":
            
            #find pt0: 
            i,j = self.find_mesh_point(pt2, self.C_pos)
            pt0 = self.C_pos[i][j-1]
            #check for intersection:
            A,B = [pt2.x, pt2.y], [pt3.x, pt3.y]
            C,D = [pt0.x, pt0.y], [pt1.x, pt1.y]

        if charDir == "pos":
            #find pt0: 
            i,j = self.find_mesh_point(pt1, self.C_neg)
            pt0 = self.C_neg[i][j-1]
            #check for intersection:
            A,B = [pt1.x, pt1.y], [pt3.x, pt3.y]
            C,D = [pt1.x, pt1.y], [pt0.x, pt0.y]

        if intersect(A,B,C,D): hasIntersected = True
        return hasIntersected
        
    def compile_mesh_points(self):
        """
        dumps all mesh points into one bucket, assigns them all indices, and makes a new triangle list of point indices
        need to run this before plotting the mesh
        """
        Clist = self.C_neg #both C_neg and C_pos should contain the same points so this should be fine. Maybe put a check here? 

        self.meshPts = []
        i = 0 
        for char in Clist:
            for pt in char: 
                pt.ind = i
                self.meshPts.append(pt)
                i += 1 

        self.triangle = []
        for i,tri in enumerate(self.triangle_obj): #!TEMPORAY IMPROVE LATER
            self.triangle.append([pt.ind if pt is not None else None for pt in tri])
            
    def find_mesh_point(self, pt, C_posneg):
        """
        gets the index of point in C_posneg (self.C_pos or self.C_neg)
        returns [i,j]
            i = line index 
            j = point index
        """
        for i,char in enumerate(C_posneg):
            for j,p in enumerate(char): 
                if p == pt:
                    return [i,j]

    def trim_mesh_after_intersect(self, pt2, pt1, charDir):
        """
        deletes downstream portion of crossed characteristic
        """
        if charDir == "neg":
            #find pt1 in neg families 
            i,j = self.find_mesh_point(pt1, self.C_neg)
            #delete points
            delPts = [pt for ii,pt in enumerate(self.C_neg[i]) if ii >= j]
        
        elif charDir == "pos":
            #find pt2 in pos families
            i,j = self.find_mesh_point(pt2, self.C_pos)
            #delete points
            delPts = [pt for ii,pt in enumerate(self.C_pos[i]) if ii >= j]

        self.delete_mesh_points(delPts)

    def delete_mesh_points(self, delPts):
        """
        #!for some reason, the C_neg and C_pos lists are not equal in size
        given a list of mesh points, this function will delete them from the mesh, line segments included
        """
        #delete from pos & neg char list 
        for i,char in enumerate(self.C_pos):
            newchar = []
            for pt in char: 
                if pt not in delPts: 
                    newchar.append(pt)
            self.C_pos[i] = newchar        
        
        for i,char in enumerate(self.C_neg):
            newchar = []
            for pt in char: 
                if pt not in delPts: 
                    newchar.append(pt)
            self.C_neg[i] = newchar   

        #delete from triangle list
        removeTriInd = []
        for i,tri in enumerate(self.triangle_obj): 
            for pt in tri:
                if pt in delPts and pt is not None: 
                    removeTriInd.append(i)
                    break
        self.triangle_obj = [tri for i,tri in enumerate(self.triangle_obj) if i not in removeTriInd]

class mesh_point: 
    def __init__(self,x,y,u,v,ind,isWall=False, isIdl=False):
         self.x,self.y,self.u,self.v = x,y,u,v 
         self.i = ind
         self.isWall = isWall #is the point on the boundary? 
         self.isIdl = isIdl #is the point on the initial data line? 

    def get_point_properties(self, gasProps): 
        #unpacking
        gam, a0, T0, p0 = gasProps.gam, gasProps.a0, gasProps.T0, gasProps.p0
        V = math.sqrt(self.u**2 + self.v**2)
        a = math.sqrt(a0**2 + 0.5*(gam-1)*V**2)
        self.mach = V/a
        self.T = T0/(1+0.5*(gam-1)*(V/a)**2)
        self.p = p0*(T0/self.T)**(gam/(gam-1))