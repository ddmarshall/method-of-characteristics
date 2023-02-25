import method_of_characteristics.unit_processes as moc_op
import math
"""
class for generating mesh and mesh point objects
"""
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
        self.alternateChar = False

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
        self.generate_initial_mesh(self.idl, charDir) #generate mesh from initial data line
        self.compile_mesh_points()
        while self.f_kill[1] == False: 
            #!following sections were written for initiating with neg lines. Untested for starting with positive
            if charDir == "neg":
                #generate initial above-wall point 
                pt1 = self.C_neg[-1][1] #2nd point in most recent characteristic
                [x3,y3,u3,v3] = moc_op.direct_wall_abv(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs)
                pt3 = mesh_point(x3,y3,u3,v3, None, isWall=True)
                self.triangle_obj.append([pt3, None, pt1]) 
                i,_ = self.find_mesh_point(pt1, self.C_pos)
                self.C_pos[i].append(pt3) #add to positive
                prev_n_char = self.C_neg[-1][2:]
                self.compute_next_neg_char(pt3, prev_n_char)
                if self.alternateChar: #if intersection occured, need to iterate over familiy to capture reflection 
                    charDir = "pos"

            elif charDir == "pos":

                pass

            if charDir == "pos" and self.alternateChar:
                #execute after finishing a crossed characteristic
                self.alternateChar = False
                self.generate_mesh_from_line(self.C_neg[-1], "pos")
                break 

            if charDir == "neg" and self.alternateChar: 
                pass 

    def generate_initial_mesh(self, dl, charDir):
        """
        computes the initial mesh from the idl. For a vertical IDL this should form a triangle with either a leading + or - characterstic spanning wall to wall 
        """
        
        if charDir == "pos": 
            self.C_pos.append([dl[0]])
            [self.C_neg.append([pt]) for pt in dl]
            self.C_neg.reverse() #first neg line at bottom 
        elif charDir == "neg":
            self.C_neg.append([dl[-1]])
            [self.C_pos.append([pt]) for pt in dl] #first pos line at the top
            dl.reverse() #reverse idl to start at bottom                
        #generating intial mesh: 
        for i,pt in enumerate(dl):
            if i == 0: continue #first point has no - char passing through it 
            if charDir == "neg":
                self.compute_next_neg_char(pt, self.C_neg[i-1])
            elif charDir == "pos":
                self.compute_next_pos_char(pt, self.C_pos[i-1])

    def generate_mesh_from_line(self, dl, charDir):
    
        if charDir == "neg": 
            dl.reverse() #reverse idl to start at bottom                
        
        if charDir == "pos":
            #top wall solution
            pt1 = dl[1]
            [x3,y3,u3,v3] = moc_op.direct_wall_abv(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs)
            pt3 = mesh_point(x3,y3,u3,v3,None,isWall=True)
            self.triangle_obj.append([pt3, None, pt1])
            self.C_neg.append([pt3])
            i,_ = self.find_mesh_point(pt1, self.C_pos)
            self.C_pos[i].append(pt3)
        
        elif charDir == "neg":
            pass 

        for i,initPt in enumerate(dl):
            if i >= 2:
                if charDir == "neg":
                    pass 
                elif charDir == "pos":
                    i,j = self.find_mesh_point(initPt, self.C_neg) 
                    pt0 = self.C_neg[i][j-1]
                    i,j = self.find_mesh_point(pt0, self.C_pos)
                    prev_p_char = self.C_pos[i][j+1:]
                    self.compute_next_pos_char(initPt, prev_p_char, continueChar=True);  

    def compute_next_neg_char(self, init_point, prev_n_char, continueChar=False):
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
                self.alternateChar = True
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

    def compute_next_pos_char(self, init_point, prev_p_char, continueChar=False):
        """
        Generates the next leading positive characteristic by advancing the mesh along the previous
        init_point could be wall or idl point
        !OUTDATED !UPDATE THIS
        """
        if continueChar is False:
            self.C_pos.append([init_point])
        for pt in prev_p_char:
            pt2 = pt #above
            pt1 = init_point #below 
            [x3, y3, u3, v3] = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.pcTOL, self.funcs)
            pt3 = mesh_point(x3, y3, u3, v3, ind = None)
            if self.check_for_int_intersect(pt3, pt2, pt1, "pos"):
                self.trim_mesh_after_intersect(pt2, pt1, "pos")
                init_point = self.C_pos[-1][-1]
                i,j = self.find_mesh_point(init_point, self.C_neg)
                pt0 = self.C_neg[i][j-2]
                i,j = self.find_mesh_point(pt0, self.C_pos)
                prev_p_char = self.C_pos[i][j+1:]
                self.compute_next_pos_char(init_point, prev_p_char, continueChar=True) #recursion ooooh spooky 
                self.alternateChar = True
                return

            if continueChar:
                i,j = self.find_mesh_point(pt1, self.C_pos) 
                self.C_pos[i].append(pt3)
            else: 
                self.C_pos[-1].append(pt3)

            self.triangle_obj.append([pt3, pt2, pt1])

            #find point 2 and append the new point 
            i,_ = self.find_mesh_point(pt2, self.C_neg)
            self.C_neg[i].append(pt3)

            init_point = pt3 #update initial point

            if self.f_kill[0](self) == True:
                self.f_kill[1] = True
                return

        #terminating wall point
        if continueChar:
            pt1 = self.C_neg[-1][-1]
            [x3,y3,u3,v3] = moc_op.direct_wall_abv(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs)
            pt3 = mesh_point(x3, y3, u3, v3, None, isWall=True)
            i,_ = self.find_mesh_point(pt1, self.C_pos)
            self.C_pos[i].append(pt3)
            
        else: 
            pt1 = self.C_pos[-1][-1]
            [x3,y3,u3,v3] = moc_op.direct_wall_abv(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs)
            pt3 = mesh_point(x3, y3, u3, v3, None, isWall=True)
            self.C_pos[-1].append(pt3)

        self.C_neg.append([pt3])
        self.triangle_obj.append([pt3, None, pt1])

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
            C,D = [pt0.x, pt0.y], [pt2.x, pt2.y]

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
                pt.i = i
                self.meshPts.append(pt)
                i += 1 

        self.triangle = []
        for i,tri in enumerate(self.triangle_obj): #!TEMPORAY IMPROVE LATER
            self.triangle.append([pt.i if pt is not None else None for pt in tri])
            
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