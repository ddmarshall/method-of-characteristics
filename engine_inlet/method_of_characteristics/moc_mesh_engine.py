import method_of_characteristics.unit_processes as moc_op
import method_of_characteristics.oblique_shock_point_irrot as shock
import math
import numpy as np
"""
Module responsible for generating method-of-characteristics mesh and mesh points
TODO: for shock mesh, make function to perform inverse wall operations in compression corners and add to shockPoint list
"""

class Mesh:

    def __init__(self, inputObj, kill_func, idl, explicit_shocks=False):
        self.C_pos, self.C_neg = [],[] #containers for characteristics lines
        self.triangle_obj = [] #container for point line segments
        self.shock_upst_segments_obj = [] #container for shock line segments
        self.shockPts_frontside = [] 
        self.shockPts_backside = []
        self.funcs = moc_op.operator_funcs() #operator functions
        self.gasProps = inputObj.gasProps
        self.delta = inputObj.delta
        self.pcTOL = inputObj.pcTOL
        self.geom = inputObj.geom
        self.f_kill = [kill_func, False]
        self.alternateChar = False
        self.numPointsGen = 0
        self.hasIntersected = False
        self.shockMesh = explicit_shocks
        self.init_method = inputObj.init_method
        self.working_region = 0 #iterator to denote multiple flowfield regions (bounded by shocks and walls)
        self.idl_p0_p0f = idl.p0_p0f
        self.idl = [Mesh_Point(x, idl.y[i], idl.u[i], idl.v[i], self.working_region, isIdl=True) for i,x in enumerate(idl.x)]

        #if hasattr(idl, 'cowlPoint'):
        self.idl[0].isWall = True 
        #if hasattr(idl, 'cbPoint'):
        self.idl[-1].isWall = True  
        self.cowl_point =  self.idl[0]

        charDir = "neg" #!hardcoded 
        if self.init_method == "STRAIGHT IDL":
            self.generate_initial_mesh_from_straight_idl(charDir) #generate initial mesh from idl
        elif self.init_method == "MACH LINE":
            self.generate_initial_mesh_from_mach_line_idl(charDir) #generate intial mesh from idl 
            
            #check if characteristic intersection with centerbody is acceptable 
            #if not: use backup function
          
        if explicit_shocks: 
            shockPt_init = self.idl[-1]
            shockPt_init.isShock = True
            self.shock_object_list = [[]]
            self.shockPts_frontside.append(shockPt_init)
            self.shock_upst_segments_obj.append(shockPt_init)
            self.generate_mesh_with_shocks_from_cowl() #generate mesh with shocks
                    
        else: 
            self.generate_mesh_from_cowl()
        
        self.compile_mesh_points()
        self.get_point_properties_and_minmax()
        self.compile_wall_points(self.geom.y_cowl, self.geom.y_centerbody) 
        self.compute_local_mass_flow()

    def generate_initial_mesh_from_straight_idl(self, charDir):
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
                if self.f_kill[0](self) == True:
                        self.f_kill[1] = True
                        return 
            elif charDir == "pos":
                self.compute_next_pos_char(pt, self.C_pos[i-1])
                if self.f_kill[0](self) == True:
                        self.f_kill[1] = True
                        return
    
    def generate_initial_mesh_from_mach_line_idl(self, charDir):
        """
        creates initial portion of the mesh from an initial data line defined 
        by a mach line extending from the cowl lip to the centerbody
        """
        idl = self.idl
        idl.reverse()
        if charDir == "pos":
            pass 
        elif charDir == "neg":
            self.C_pos.append([pt for pt in idl])
            [self.C_neg.append([pt]) for pt in idl]

        #generating initial mesh
        for i,pt in enumerate(idl):
            if i==0: continue
            if i==1:
                [x3,y3,u3,v3] = moc_op.direct_wall(pt, self.geom.y_centerbody,\
                                self.geom.dydx_centerbody, self.gasProps, \
                                    self.delta, self.pcTOL, self.funcs, "neg")
                pt3 = Mesh_Point(x3,y3,u3,v3, self.working_region, isWall=True)
                self.triangle_obj.append([pt3, None, pt])
                i,_ = self.find_mesh_point(pt, self.C_neg)
                self.C_neg[i].append(pt3)
                self.C_pos.append([pt3])
            else: 
                self.compute_next_neg_char(pt, self.C_neg[i-1][1:], continueChar=True) 
    
    def generate_initial_char_from_cb_point(self):

        """
        generates an initial positive characteristic from the centerbody. Uses 
        inverse wall to connect mesh to cowl point. Use this when 
        generate_initial_char_from_cowl_lip fails to capture non-straight 
        geometry. 
        """

    def generate_mesh_from_cowl(self):
        """
        generates shock-less characteristic mesh until the kill function is triggered 
        Handles shock waves implicitly i.e. trims mesh after same-family intersection and does NOT include oblique shock calculations
        """
        charDir = "neg" #!hard coded for now... starting direction
        self.compile_mesh_points()
        while self.f_kill[1] == False: 
            #!following sections were written for initiating with neg lines. Untested for starting with positive
            if charDir == "neg" and self.alternateChar == False:
                #TODO use generate_mesh_from_line here instead 
                
                #generate initial above-wall point 
                pt1 = self.C_neg[-1][1] #2nd point in most recent characteristic
                [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
                pt3 = Mesh_Point(x3,y3,u3,v3, self.working_region, isWall=True)
                self.triangle_obj.append([pt3, None, pt1]) 
                i,_ = self.find_mesh_point(pt1, self.C_pos)
                self.C_pos[i].append(pt3) #add to positive
                prev_n_char = self.C_neg[-1][2:]
                self.compute_next_neg_char(pt3, prev_n_char)
                if self.alternateChar: #if intersection occured, need to iterate over family to capture reflection 
                    charDir = "pos"

            elif charDir == "pos" and self.alternateChar == False:
                #TODO add if starting with positive characteristics
                pass

            if charDir == "pos" and self.alternateChar:
                #execute after finishing a crossed characteristic
                self.alternateChar = False
                if self.hasIntersected == True: 
                    self.generate_until_non_intersect(self.C_neg[-1], "neg")  
                self.generate_mesh_from_line(self.C_neg[-1], "pos")    
                if self.alternateChar: 
                    charDir = "neg"  
         
            if charDir == "neg" and self.alternateChar: 
                self.alternateChar = False
                if self.hasIntersected == True: 
                    self.generate_until_non_intersect(self.C_pos[-1], "pos") 
                self.generate_mesh_from_line(self.C_pos[-1],"neg")   
                if self.alternateChar:
                    charDir = "pos" 
                break  

    def generate_mesh_with_shocks_from_cowl(self):
        """
        generates characteristic mesh using oblique shock points
        """
        charDir = "neg"
        if charDir == "neg": charList = self.C_neg
        elif charDir == "pos": charList = self.C_pos
        if self.check_for_passed_shock_point(charList) == False:
            print("Shock Point Not Within Mesh")
            return  
        try: self.compute_wall_to_wall_shock(charDir, self.shockPts_frontside[-1]) 
        except: return

        if self.impending_shock_reflec:
            while True: 
                
                if charDir == "pos": charDir = "neg"
                elif charDir == "neg": charDir = "pos"

                try: self.generate_mesh_for_shock_reflec(charDir)  
                except: return 
                
                try: self.compute_wall_to_wall_shock(charDir, self.shockPts_backside[-1])
                except: return 

    def generate_mesh_from_line(self, dl, charDir):
        """
        Generates the characteristic mesh from an existing characteristic (either pos or neg) terminates when all opposite characteristics originating 
        from dl have been generated
        dl = list of mesh points from positive/negative characteristic spanning wall to wall 
        """
             
        if charDir == "pos":
            #top wall solution
            pt1 = dl[1]
            [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
            if self.f_kill[0](self) == True:
                self.f_kill[1] = True
                return 
            pt3 = Mesh_Point(x3,y3,u3,v3,self.working_region,isWall=True)
            self.triangle_obj.append([pt3, None, pt1])
            self.C_neg.append([pt3])
            i,_ = self.find_mesh_point(pt1, self.C_pos)
            self.C_pos[i].append(pt3)
        
        elif charDir == "neg":
            #bottom wall solution
            pt2 = dl[1]
            [x3,y3,u3,v3] = moc_op.direct_wall(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs, "neg")
            if self.f_kill[0](self) == True:
                self.f_kill[1] = True
                return 
            pt3 = Mesh_Point(x3,y3,u3,v3,self.working_region,isWall=True)
            self.triangle_obj.append([pt3, pt2, None])
            self.C_pos.append([pt3])
            i,_ = self.find_mesh_point(pt2, self.C_neg)
            self.C_neg[i].append(pt3) 

        for ind,initPt in enumerate(dl):
            
            if ind >= 2:
                self.hasIntersected = False  
                if charDir == "neg":
                    i,j = self.find_mesh_point(initPt, self.C_pos) 
                    pt0 = self.C_pos[i][j-1]
                    i,j = self.find_mesh_point(pt0, self.C_neg)
                    prev_n_char = self.C_neg[i][j+1:]
                    try: self.compute_next_neg_char(initPt, prev_n_char, continueChar=True);
                    except: #if math error occurs (due to subsonic) kill mesher 
                        self.f_kill[1] = True
                        return 
                    
                    if self.f_kill[0](self) == True:
                        self.f_kill[1] = True
                        return 
 
                elif charDir == "pos":
                    i,j = self.find_mesh_point(initPt, self.C_neg)
                    pt0 = self.C_neg[i][j-1]
                    i,j = self.find_mesh_point(pt0, self.C_pos)
                    prev_p_char = self.C_pos[i][j+1:]
                    try:self.compute_next_pos_char(initPt, prev_p_char, continueChar=True);  
                    except: #if math error occurs (due to subsonic) kill mesher 
                        self.f_kill[1] = True
                        return

                    if self.f_kill[0](self) == True:
                        self.f_kill[1] = True
                        return 

    def generate_mesh_for_shock_reflec(self, charDir):
        
        if charDir == "pos":
            C_on = self.C_neg
            C_off = self.C_pos
            y_x, dydx = self.geom.y_cowl, self.geom.dydx_cowl
        elif charDir == "neg":
            C_on = self.C_pos
            C_off = self.C_neg
            y_x, dydx = self.geom.y_centerbody, self.geom.dydx_centerbody

        pointList = C_on[-1]
        count = 0
        while len(pointList) > 2:

            pt = pointList[1]
            ii,_ = self.find_mesh_point(pt, C_off)
            [x3, y3, u3, v3] = moc_op.direct_wall(pt, y_x, dydx, self.gasProps, self.delta, self.pcTOL, self.funcs, charDir)
            init_point = Mesh_Point(x3, y3, u3, v3, self.working_region, isWall=True)
            C_on.append([init_point]), C_off[ii].append(init_point), self.triangle_obj.append([init_point, pt, None])
            dataLine = pointList[2:]


            if charDir=="pos": 
                C_on, C_off = self.compute_next_neg_char(init_point, dataLine, onChars=C_on, offChars=C_off, continueChar=True, terminate_wall=False, check_for_intersect=False)
            elif charDir=="neg":   
                C_on, C_off = self.compute_next_pos_char(init_point, dataLine, onChars=C_on, offChars=C_off, continueChar=True, terminate_wall=False, check_for_intersect=False)

            pointList = C_on[-1]
            count += 1
        
        pt = pointList[1]
        ii,_ = self.find_mesh_point(pt, C_off)
        [x3, y3, u3, v3] = moc_op.direct_wall(pt, y_x, dydx, self.gasProps, self.delta, self.pcTOL, self.funcs, charDir)
        init_point = Mesh_Point(x3, y3, u3, v3, self.working_region, isWall=True)
        C_on.append([init_point]), C_off[ii].append(init_point), self.triangle_obj.append([init_point, pt, None])

        if charDir == "pos":
            self.C_neg = C_on
            self.C_pos = C_off
        elif charDir == "neg":
            self.C_pos = C_on
            self.C_neg = C_off

    def compute_next_neg_char(self, init_point, prev_n_char, onChars=None, offChars=None, continueChar=False, terminate_wall=True, check_for_intersect=True):
        """
        Generates the next leading negative characteristic by advancing the mesh along the previous one
        init_point could be wall or idl point
        !NOTE TO SELF: CHANGES HERE NEED TO BE REFLECTED IN TWIN FUNCTION
        """
        if onChars is not None and offChars is not None: 
            #if supplying existing characteristic list that is not self.C_neg or self.C_pos
            C_on, C_off = onChars, offChars
        else: 
            C_on, C_off = self.C_neg, self.C_pos

        if continueChar is False:
            C_on.append([init_point])
        for pt in prev_n_char:
            pt2 = init_point #above
            pt1 = pt #below 
            [x3, y3, u3, v3] = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.pcTOL, self.funcs)
            pt3 = Mesh_Point(x3, y3, u3, v3, self.working_region)
            if check_for_intersect: 
                if self.check_for_int_intersect(pt3, pt2, pt1, "neg"):
                    self.hasIntersected = True #indicate that an intersection has occurred 
                    self.trim_mesh_after_intersect(pt2, pt1, "neg")
                    init_point = pt2
                    i,j = self.find_mesh_point(init_point, self.C_pos)
                    pt0 = self.C_pos[i][j-2]
                    i,j = self.find_mesh_point(pt0, self.C_neg)
                    prev_n_char = self.C_neg[i][j+1:]
                    self.compute_next_neg_char(init_point, prev_n_char, continueChar=True) #recursion ooooh spooky 
                    self.alternateChar = True
                    return

            if continueChar: 
                i,j = self.find_mesh_point(pt2, C_on)
                C_on[i].append(pt3)
            else: 
                C_on[-1].append(pt3)

            self.triangle_obj.append([pt3, pt2, pt1])

            #find point 1 and append the new point 
            i,_ = self.find_mesh_point(pt1, C_off)
            C_off[i].append(pt3)

            init_point = pt3 #update initial point

        #terminating wall point
        if terminate_wall: 
            if continueChar:
                pt2 = C_off[-1][-1]
                [x3,y3,u3,v3] = moc_op.direct_wall(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs, "neg")
                pt3 = Mesh_Point(x3, y3, u3, v3, self.working_region, isWall=True)
                i,_ = self.find_mesh_point(pt2, C_on)
                C_on[i].append(pt3)

            else:
                pt2 = C_on[-1][-1]
                [x3,y3,u3,v3] = moc_op.direct_wall(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs, "neg")
                pt3 = Mesh_Point(x3, y3, u3, v3, self.working_region, isWall=True)
                C_on[-1].append(pt3)

            C_off.append([pt3])
            self.triangle_obj.append([pt3, pt2, None])

        if onChars is not None and offChars is not None: return C_on, C_off
        else: self.C_neg, self.C_pos = C_on, C_off

    def compute_next_pos_char(self, init_point, prev_p_char, onChars=None, offChars=None, continueChar=False, terminate_wall=True, check_for_intersect=True):
        """
        Generates the next leading positive characteristic by advancing the mesh along the previous one
        init_point could be wall or idl point
        """
        if onChars is not None and offChars is not None: 
            #if supplying existing characteristic list that is not self.C_neg or self.C_pos
            C_on, C_off = onChars, offChars
        else: 
            C_on, C_off = self.C_pos, self.C_neg

        if continueChar is False:
            C_on.append([init_point])
        for pt in prev_p_char:
            pt2 = pt #above
            pt1 = init_point #below 
            [x3, y3, u3, v3] = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.pcTOL, self.funcs)
            pt3 = Mesh_Point(x3, y3, u3, v3, self.working_region)
            if check_for_intersect: 
                if self.check_for_int_intersect(pt3, pt2, pt1, "pos"):
                    self.hasIntersected = True 
                    self.trim_mesh_after_intersect(pt2, pt1, "pos")
                    init_point = pt1
                    i,j = self.find_mesh_point(init_point, C_off)
                    pt0 = C_off[i][j-2]
                    i,j = self.find_mesh_point(pt0, C_on)
                    prev_p_char = C_on[i][j+1:]
                    self.compute_next_pos_char(init_point, prev_p_char, continueChar=True) #recursion ooooh spooky 
                    self.alternateChar = True
                    return

            if continueChar:
                i,j = self.find_mesh_point(pt1, C_on)
                C_on[i].append(pt3)
            else: 
                C_on[-1].append(pt3)

            self.triangle_obj.append([pt3, pt2, pt1])

            #find point 2 and append the new point 
            i,_ = self.find_mesh_point(pt2, C_off)
            C_off[i].append(pt3)

            init_point = pt3 #update initial point

        #terminating wall point
        if terminate_wall:
            if continueChar:
                pt1 = C_off[-1][-1]
                [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
                pt3 = Mesh_Point(x3, y3, u3, v3, self.working_region, isWall=True)
                i,_ = self.find_mesh_point(pt1, C_on)
                C_on[i].append(pt3)

            else: 
                pt1 = C_on[-1][-1]
                [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
                pt3 = Mesh_Point(x3, y3, u3, v3, self.working_region, isWall=True)
                C_on[-1].append(pt3)

            C_off.append([pt3])
            self.triangle_obj.append([pt3, None, pt1])

        if onChars is not None and offChars is not None: return C_on, C_off
        else: self.C_pos, self.C_neg = C_on, C_off

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

        elif charDir == "pos":
            #find pt0: 
            i,j = self.find_mesh_point(pt1, self.C_neg)
            pt0 = self.C_neg[i][j-1]
            #check for intersection:
            A,B = [pt1.x, pt1.y], [pt3.x, pt3.y]
            C,D = [pt0.x, pt0.y], [pt2.x, pt2.y]

        if intersect(A,B,C,D): hasIntersected = True
        return hasIntersected

    def check_for_wall_intersect(self, pt3, pt12, charDir):
        """
        checks for a same family intersection at the generation of a new wall point and an existing wall point
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
            i,j = self.find_mesh_point(pt12, self.C_pos)
            pt0 = self.C_pos[i][j-1]
            i,j = self.find_mesh_point(pt0, self.C_neg)
            pt1 = self.C_neg[i][j]

            A,B = [pt12.x, pt12.y], [pt3.x, pt3.y]
            C,D = [pt0.x, pt0.y], [pt1.x, pt1.y]

        elif charDir == "pos":
            #TODO code this 
            pass 

        if intersect(A,B,C,D): hasIntersected = True 
        return hasIntersected

    def generate_until_non_intersect(self, dl, charDir):
        """
        will generate successive + or - chars until one makes it wall to wall without a same-family intersection
        """
        while self.hasIntersected == True: 
            self.hasIntersected = False #set to false
            if charDir == "neg":
                pt1 = self.C_neg[-1][1] #2nd point in most recent characteristic
                [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
                pt3 = Mesh_Point(x3,y3,u3,v3,self.working_region,isWall=True)
                self.triangle_obj.append([pt3, None, pt1]) 
                i,_ = self.find_mesh_point(pt1, self.C_pos)
                self.C_pos[i].append(pt3) #add to positive
                init_point = pt3
                prev_n_char = self.C_neg[-1][2:]
                self.compute_next_neg_char(init_point, prev_n_char)

            elif charDir == "pos":
                
                pt2 = self.C_pos[-1][1] #2nd point in most recent characteristic
                [x3,y3,u3,v3] = moc_op.direct_wall(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs, "neg")
                pt3 = Mesh_Point(x3,y3,u3,v3,self.working_region,isWall=True)
                self.triangle_obj.append([pt3, pt2, None]) 
                i,_ = self.find_mesh_point(pt2, self.C_neg)
                self.C_neg[i].append(pt3) #add to positive
                init_point = pt3
                prev_p_char = self.C_pos[-1][2:] 
                self.compute_next_pos_char(init_point, prev_p_char) 

    def check_for_passed_shock_point(self, C_posneg):
        """
        checks if the mesh has passed a known wall shock point. To be used as 
        preliminary check for wall-to-wall shock function
        """
        ret = False
        leading_char = C_posneg[-1]
        wallPts = [pt for pt in leading_char if pt.isWall==True]
        
        if wallPts[0].x >= self.shockPts_frontside[-1].x and wallPts[-1].x >= self.shockPts_frontside[-1].x: 
            #!logic will not work if flow has negative x-velocity component
            ret = True
                 
        return ret 

    def compute_wall_to_wall_shock(self, shockDir, init_shock_point):
        """
        computes a shock wave from one wall to another, deleting and appending mesh points at it moves along
        """
        if shockDir == "neg":
            C_on, C_off = self.C_neg, self.C_pos
            y_x_i, dydx_i = self.geom.y_cowl, self.geom.dydx_cowl
            y_x_f, dydx_f = self.geom.y_centerbody, self.geom.dydx_centerbody
        elif shockDir == "pos": 
            C_on, C_off = self.C_pos, self.C_neg 
            y_x_i, dydx_i = self.geom.y_centerbody, self.geom.dydx_centerbody
            y_x_f, dydx_f = self.geom.y_cowl, self.geom.dydx_cowl

        #GENERATE FIRST SHOCK POINT#############################################
        pt_w_ups = init_shock_point
        [i,j] = self.find_mesh_point(pt_w_ups, C_on)
        pt = C_on[i][j+1]
        [ii,jj] = self.find_mesh_point(pt, C_off)
        pt1 = C_off[ii][jj-1]
        [i_,j_] = self.find_mesh_point(pt1, C_on)
        pt0 = C_on[i_][j_-1]

        [pt4_dwn, pt4_ups, def4, beta4, ptw_dwn, pt3p, shockObj_wall, shockObj] = shock.wall_shock_point(pt_w_ups, y_x_i, dydx_i, pt1, self.pcTOL, self.delta, self.gasProps, shockDir)        
        
        #check if new shock point intersects same-family characteristic:
        delPts = None
        if self.check_for_shock_char_intersect(pt0, pt1, pt_w_ups, shockPt2=pt4_ups):
            pt1 = C_off[ii][jj-2]
            [pt4_dwn, pt4_ups, def4, beta4, ptw_dwn, pt3p, shockObj_wall, shockObj] = shock.wall_shock_point(pt_w_ups, y_x_i, dydx_i, pt1, self.pcTOL, self.delta, self.gasProps, shockDir)
            delPts = C_on[i_][j_:]
        
        pt4_dwn = Mesh_Point(pt4_dwn.x, pt4_dwn.y, pt4_dwn.u, pt4_dwn.v, self.working_region+1, isShock=True)
        pt4_ups = Mesh_Point(pt4_ups.x, pt4_ups.y, pt4_ups.u, pt4_ups.v, self.working_region, isShock=True)
        ptw_dwn = Mesh_Point(ptw_dwn.x, ptw_dwn.y, ptw_dwn.u, ptw_dwn.v, self.working_region+1, isShock=True, isWall=True)
        pt3p = Mesh_Point(pt3p.x, pt3p.y, pt3p.u, pt3p.v, self.working_region+1, isWall=True)
        self.shock_object_list[-1].append(shockObj_wall), self.shock_object_list[-1].append(shockObj)
        
        self.shock_upst_segments_obj.append(pt4_ups) #add shock pt4_ups to segment list 
        [self.shockPts_backside.append(pt) for pt in [ptw_dwn, pt4_dwn]]
        self.shockPts_frontside.append(pt4_ups)

        self.triangle_obj.append([pt4_ups, pt1, pt_w_ups]) #add mesh segments to triangle list 
        self.triangle_obj.append([pt3p, pt4_dwn, None]) 

        C_off[ii-1].append(ptw_dwn)
        C_off[ii][jj] = pt4_ups #replace old point with shock point
        C_off[ii].append(pt4_dwn)
        C_off[ii].append(pt3p) #add new point 

        C_on[i][j+1] = pt4_ups #replace old point with shock point
        C_on.append([ptw_dwn])
        C_on[-1].append(pt4_dwn)
        C_on.append([pt3p])
        
        pt_s_ups, pt_s_dwn = pt4_ups, pt4_dwn 

        if delPts is not None: 
            self.delete_mesh_points(delPts, C_on=C_on, C_off=C_off)
        
        #INTERIOR POINT LOOP####################################################
        count = 0 
        while True: 
            
            [i,j] = self.find_mesh_point(pt_s_ups, C_on)
            pt = C_on[i][j+1]
            [ii,jj] = self.find_mesh_point(pt, C_off)
            pt1_old = pt1
            pt1 = C_off[ii][jj-1]

            #check for and handle intersection of shock wave with upstream (same-family) characteristic
            delPts = []
            if self.check_for_shock_char_intersect(pt1_old, pt1, pt_s_ups, waveAng=beta4): 
                pt1 = self.get_upstream_point_for_shock(C_off[ii], pt_s_ups, beta4) #finds new pt1
                if pt1 is None: 
                    pt1 = pt1_old
                    break #will be none if the next place for shock to hit is the wall 
                [iii, jjj] = self.find_mesh_point(pt1_old, C_on)
                [delPts.append(pt) for pt in C_on[iii][jjj+1:]]

            pt_a = pt3p
            beta4_old, def4_old = beta4, def4
            [pt4_dwn, pt4_ups, def4, beta4, pt3p, shockObj] = shock.interior_shock_point(pt_s_ups, pt_s_dwn, beta4, def4, pt1, pt_a, self.pcTOL, self.delta, self.gasProps, shockDir)
            
            #check for and handle intersection of shock wave with downstream (same-family) characteristic
            while self.check_for_shock_char_intersect(pt3p, pt_a, pt_s_ups, shockPt2=pt4_ups):
                if pt_a.isWall: 
                    [pt4_dwn, pt4_ups, def4, beta4, pt3p, pt_a, shockObj, C_on, C_off] = self.handle_from_wall_shock_char_intersect(pt_s_ups, pt_s_dwn, beta4_old, def4_old, pt_a, pt1, pt3p, C_on, C_off, shockDir) 
                    ii += 1 #update ii due to new characteristic being inserted between s and 4
                else:  
                    [pt4_dwn, pt4_ups, def4, beta4, pt3p, pt_a, shockObj, C_on, C_off] = self.handle_interior_shock_char_intersect(pt_s_ups, pt_s_dwn, beta4_old, def4_old, pt1, pt_a, C_on, C_off, shockDir)

            #print(f"\tshock wave angle {math.degrees(beta4)} deg")
            pt4_dwn = Mesh_Point(pt4_dwn.x, pt4_dwn.y, pt4_dwn.u, pt4_dwn.v, self.working_region+1, isShock=True)
            pt4_ups = Mesh_Point(pt4_ups.x, pt4_ups.y, pt4_ups.u, pt4_ups.v, self.working_region, isShock=True)
            pt3p = Mesh_Point(pt3p.x, pt3p.y, pt3p.u, pt3p.v, self.working_region+1)
            self.shock_object_list[-1].append(shockObj)
            
            self.shock_upst_segments_obj.append(pt4_ups)
            self.shockPts_backside.append(pt4_dwn)
            self.shockPts_frontside.append(pt4_ups)

            self.triangle_obj.append([pt4_ups, pt1, pt_s_ups])
            self.triangle_obj.append([pt3p, pt4_ups, pt_a])

            C_off[ii][jj] = pt4_ups #replace old point with shock point
            C_off[ii].append(pt4_dwn)
            C_off[ii].append(pt3p) #add new point  

            C_on[i][j+1] = pt4_ups #replace old point with shock point
            C_on[i+1].append(pt4_dwn)
            C_on[-1].append(pt3p)     
            if len(delPts) > 0: self.delete_mesh_points(delPts, C_on=C_on, C_off=C_off)

            pt_s_ups, pt_s_dwn = pt4_ups, pt4_dwn
            #check if next point is on the wall 
            count += 1
            if C_off[ii][jj-1].isWall:
                break
        
        #TO WALL SHOCK POINT####################################################
        [i,j] = self.find_mesh_point(pt_s_ups, C_on)
        pt = C_on[i][j+1]
        [ii,jj] = self.find_mesh_point(pt, C_off)
        pt1_old = pt1
        [iii,jjj] = self.find_mesh_point(pt1_old, C_on) #old point 1 location
        delPts = []
        try: 
            #will fail if a handle_shock_char_intersect occurred on previous interior point
            pt1 = C_on[iii][jjj+1]
            if self.check_for_shock_char_intersect(pt1_old, pt1, pt_s_ups, waveAng=beta4): 
                delPts.append(pt1)
                pt1 = C_off[ii-1][0]
        except: pass #use previous point-1 which should be a wall point
        
        #generate to-wall shock point
        pt_a = pt3p        
        pt_s_ups, pt_s_dwn = pt4_ups, pt4_dwn
        beta4_old, def4_old = beta4, def4
        [pt4_dwn, pt4_ups, def4, beta4, reflec, pt3p, shockObj] = shock.to_wall_shock_point(pt_s_ups, pt_s_dwn, beta4, def4, pt1, pt_a, y_x_f, dydx_f, self.pcTOL, self.delta, self.gasProps, shockDir)

        #check for same-family shock characteristic intersection on downstream side
        if self.check_for_shock_char_intersect(pt3p, pt_a, pt_s_ups, shockPt2=pt4_ups):
                [pt4_dwn, pt4_ups, def4, beta4, reflec, pt3p, pt_a, shockObj, C_on, C_off] = self.handle_wall_shock_char_intersect(pt_s_ups, pt_s_dwn, beta4_old, def4_old, pt1, pt_a, C_on, C_off, y_x_f, dydx_f, shockDir)
                
        pt4_dwn = Mesh_Point(pt4_dwn.x, pt4_dwn.y, pt4_dwn.u, pt4_dwn.v, self.working_region+1, isShock=True, isWall=True)
        pt4_ups = Mesh_Point(pt4_ups.x, pt4_ups.y, pt4_ups.u, pt4_ups.v, self.working_region, isShock=True, isWall=True)
        pt3p = Mesh_Point(pt3p.x, pt3p.y, pt3p.u, pt3p.v, self.working_region+1)
        self.shock_object_list[-1].append(shockObj)

        self.shock_upst_segments_obj.append(pt4_ups)
        self.shockPts_backside.append(pt4_dwn)
        self.shockPts_frontside.append(pt4_ups)

        self.triangle_obj.append([pt4_ups, None, pt_s_ups])
        self.triangle_obj.append([pt3p, pt4_ups, pt_a])

        C_off[ii][jj] = pt4_ups #replace old point with shock point
        C_off[ii].append(pt4_dwn)
        C_off[ii].append(pt3p) #add new point  

        C_on[i][j+1] = pt4_ups #replace old point with shock point
        C_on[-2].append(pt4_dwn)
        C_on[-1].append(pt3p)
        if len(delPts) > 0: self.delete_mesh_points(delPts, C_on=C_on, C_off=C_off)

        self.delete_mesh_points([pt for pt in C_on[i] if pt.isShock == False], C_on=C_on, C_off=C_off) #bit of a dirty fix
        self.impending_shock_reflec = reflec

        if shockDir == "neg": self.C_neg, self.C_pos = C_on, C_off
        elif shockDir == "pos": self.C_pos, self.C_neg = C_on, C_off

        self.shock_object_list.append([]) #add new container for next shock
        self.working_region += 1 #onto new region:

    def get_upstream_point_for_shock(self, charPts, pt_s, beta):
        """
        finds pt1 input for shock operators given existing shock point, wave angle, and list of points on a characteristic line
        returns a single mesh point
        """
        retPt = None
        for i,pt2 in enumerate(charPts):
            if i == 0: continue 
            pt1 = charPts[i-1]
            m = (pt2.y-pt1.y)/(pt2.x-pt1.x)
            a = np.array([[1, -m],[1, -math.tan(beta)]])
            b = np.array([pt1.y - m*pt1.x, pt_s.y - math.tan(beta)*pt_s.x])
            y,x = np.linalg.solve(a,b)
            r2s = np.array([pt2.y-y, pt2.x-x])
            r1s = np.array([pt1.y-y, pt1.x-x])
            if np.dot(r2s, r1s) <= 0: 
                #if point is inbetween return downstream-most point
                if np.dot(np.array([pt1.u, pt1.v]), np.array([pt2.x-pt1.x, pt2.y-pt1.y])) >= 0:
                    retPt = pt1
                else: 
                    retPt = pt2
                break 
 
        return retPt

    def check_for_shock_char_intersect(self, charPt1, charPt2, shockPt1, shockPt2=None, waveAng=None):
        """
        checks for a intersection of a shock wave and a same-family characteristic line. Returns true if an intersection occurs. Otherwise returns false
        """
        if shockPt2 is not None:
            def ccw(A,B,C):
                return (C[1]-A[1])*(B[0]-A[0]) > (B[1]-A[1])*(C[0]-A[0])
            # Return true if line segments AB and CD intersect
            def intersect(A,B,C,D):
                #return true is segments A-B and C-D intersect 
                return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D) 

            A,B = [charPt1.x, charPt1.y],[charPt2.x, charPt2.y]
            C,D = [shockPt2.x, shockPt2.y],[shockPt1.x, shockPt1.y] 
            return intersect(A,B,C,D)
        elif waveAng is not None: 

            m = (charPt1.y-charPt2.y)/(charPt1.x-charPt2.x)
            a = np.array([[1, -m],[1, -math.tan(waveAng)]])
            b = np.array([charPt1.y - m*charPt1.x, shockPt1.y - math.tan(waveAng)*shockPt1.x])
            y,x = np.linalg.solve(a,b)
            r2s = np.array([charPt2.y-y, charPt2.x-x])
            r1s = np.array([charPt1.y-y, charPt1.x-x])
            if np.dot(r2s, r1s) <= 0: 
                #if point is inbetween return downstream-most point
                return True
            
            if np.linalg.norm(r2s) < 1e-5 or np.linalg.norm(r1s) < 1e-5:
                #if point is really close to an endpoint, will interpret as intersect
                return True
            
            return False

    def check_for_coalescence(self, charList1, charList2):
        """
        checks if two characteristic lines intersect. If there is an intersect
        returns i1, i2: indices in charList1 and charList2 after which the 
        intersection occurs. Returns False otherwide. 
        """
        def do_segments_intersect(p1, p2, p3, p4):
            def ccw(A, B, C):
                return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x)

            return ccw(p1, p3, p4) != ccw(p2, p3, p4) and ccw(p1, p2, p3) != ccw(p1, p2, p4)
                
        for i, p1 in enumerate(charList1[:-1]):
            p2 = charList1[i + 1]

            for j, p3 in enumerate(charList2[:-1]):
                p4 = charList2[j + 1]

                if do_segments_intersect(p1, p2, p3, p4):
                    return i, j

        return None,None

    def handle_from_wall_shock_char_intersect(self, pt_s_ups, pt_s_dwn, beta_s, def_s, pt_a, pt1, pt3p, C_on, C_off, shockDir):
        """
        generates a downstream wall point to deal with an intersection of a shock segment and same family characteristic when point a is a wall point
        """
        if shockDir == "neg":
            offDir = "pos"
            y_x, dydx = self.geom.y_cowl, self.geom.dydx_cowl
        elif shockDir == "pos":
            offDir = "neg"
            y_x, dydx = self.geom.y_centerbody, self.geom.dydx_centerbody

        #get position intersection of a-3' line and shock wave
        m = (pt3p.y-pt_a.y)/(pt3p.x-pt_a.x)
        a = np.array([[1, -m],[1, -math.tan(beta_s)]])
        b = np.array([pt_a.y - m*pt_a.x, pt_s_dwn.y - math.tan(beta_s)*pt_s_dwn.x])
        y,x = np.linalg.solve(a,b)
        u,v = linear_interpolate(x, pt_a.u, pt3p.u, pt_a.x, pt3p.x), linear_interpolate(x, pt_a.v, pt3p.v, pt_a.x, pt3p.x)
        intersecPt = Mesh_Point(x,y,u,v,self.working_region,)

        [x_w, y_w, u_w, v_w] = moc_op.direct_wall(intersecPt, y_x, dydx, self.gasProps, self.delta, self.pcTOL, self.funcs, offDir)
        newPt = Mesh_Point(x_w,y_w,u_w,v_w,self.working_region,isWall=True)

        pt_a_upd = newPt
        [pt4_dwn, pt4_ups, def4, beta4, pt3p, shockObj] = shock.interior_shock_point(pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a_upd, self.pcTOL, self.delta, self.gasProps, shockDir)
        
        i,_ = self.find_mesh_point(pt_a, C_on)
        ii,_ = self.find_mesh_point(pt_a, C_off)

        C_on[i].append(intersecPt)
        C_on.append([newPt])
        C_off.insert(ii+1, [intersecPt])
        C_off[ii+1].append(newPt)
        self.triangle_obj.append([intersecPt, pt_a, newPt])

        return [pt4_dwn, pt4_ups, def4, beta4, pt3p, pt_a_upd,shockObj, C_on, C_off]
    
    def handle_interior_shock_char_intersect(self, pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a, C_on, C_off, shockDir):
        """
        generates a downstream characteristic to handle intersection of a shock segment and same-family characteristic. To be called after check_for_shock_char_intersect returns
        true
        """
        if shockDir == "neg":
            offDir = "pos"
            y_x, dydx = self.geom.y_cowl, self.geom.dydx_cowl
        elif shockDir == "pos":
            offDir = "neg"
            y_x, dydx = self.geom.y_centerbody, self.geom.dydx_centerbody
            
        i,j = self.find_mesh_point(pt_a, C_on)
        pointList = C_on[i][2:j+1]
        pt = C_on[i][1]
        ii,_ = self.find_mesh_point(pt, C_off)
        [x3, y3, u3, v3] = moc_op.direct_wall(pt, y_x, dydx, self.gasProps, self.delta, self.pcTOL, self.funcs, offDir) #create first wall point
        init_point = Mesh_Point(x3, y3, u3, v3, self.working_region+1, isWall=True)
        C_on.append([init_point]), C_off[ii].append(init_point), self.triangle_obj.append([init_point, pt, None])

        if shockDir=="neg":
            C_on, C_off = self.compute_next_neg_char(init_point, pointList, onChars=C_on, offChars=C_off, continueChar=True, terminate_wall=False, check_for_intersect=False)
        elif shockDir=="pos":
            C_on, C_off = self.compute_next_pos_char(init_point, pointList, onChars=C_on, offChars=C_off, continueChar=True, terminate_wall=False, check_for_intersect=False)
        
        #check if new mach line coalesces with a same family mach line before 
        #reaching the shock point
        line1, line2 = C_on[-1], C_on[-2]
        i1,i2 = self.check_for_coalescence(line1, line2)

        while None not in [i1,i2]:
            #if an intersection was detected enter this loop until it's fixed
            #delete invalid points from last characteristic
            self.delete_mesh_points(line1[i1+1:], C_on=C_on, C_off=C_off)
            #direct wall point
            pt = C_on[-1][1]
            I,_ = self.find_mesh_point(pt, C_off)
            [x3, y3, u3, v3] = moc_op.direct_wall(pt, y_x, dydx, self.gasProps, self.delta, self.pcTOL, self.funcs, offDir) #create first wall point
            init_point = Mesh_Point(x3, y3, u3, v3, self.working_region+1, isWall=True)
            C_on.append([init_point]), C_off[I].append(init_point), self.triangle_obj.append([init_point, pt, None])
            
            #compute next char
            pointList = C_on[-2][2:] + line2[i2+1:] #create point list from 
            if shockDir=="neg":
                C_on, C_off = self.compute_next_neg_char(init_point, pointList, onChars=C_on, offChars=C_off, continueChar=True, terminate_wall=False, check_for_intersect=False)
            elif shockDir=="pos":
                C_on, C_off = self.compute_next_pos_char(init_point, pointList, onChars=C_on, offChars=C_off, continueChar=True, terminate_wall=False, check_for_intersect=False)
            
            #check for coalescence
            line1, line2 = C_on[-1], C_on[-2]
            i1, i2 = self.check_for_coalescence(line1, line2)

        pt_a_upd = C_on[-1][-1]
        [pt4_dwn, pt4_ups, def_4, beta4, pt3p, shockObj] = shock.interior_shock_point(pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a_upd, self.pcTOL, self.delta, self.gasProps, shockDir)

        return [pt4_dwn, pt4_ups, def_4, beta4, pt3p, pt_a_upd, shockObj, C_on, C_off]

    def handle_wall_shock_char_intersect(self, pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a, C_on, C_off, y_x, dydx, shockDir):
        """
        generates a downstream characteristic to handle intersection of a shock segment and same-family characteristic at the wall. To be called after check_for_shock_char_intersect returns
        true
        """
        if shockDir == "neg":
            offDir = "pos"
            y_i = self.geom.y_cowl
            dydx_i = self.geom.dydx_cowl

        elif shockDir == "pos":
            offDir = "neg"
            y_i = self.geom.y_centerbody
            dydx_i = self.geom.dydx_centerbody

        i,j = self.find_mesh_point(pt_a, C_on)
        pointList = C_on[i][2:j+1]
        pt = C_on[i][1] #first interior point
        ii,_ = self.find_mesh_point(pt, C_off)
        [x3, y3, u3, v3] = moc_op.direct_wall(pt, y_i, dydx_i, self.gasProps, self.delta, self.pcTOL, self.funcs, offDir) #create first wall point
        init_point = Mesh_Point(x3, y3, u3, v3, self.working_region+1, isWall=True)
        C_on.append([init_point]), C_off[ii].append(init_point), self.triangle_obj.append([init_point, pt, None])
        
        if shockDir == "neg":
            C_on, C_off = self.compute_next_neg_char(init_point, pointList, onChars=C_on, offChars=C_off, continueChar=True, terminate_wall=False, check_for_intersect=False)
        elif shockDir == "pos":
            C_on, C_off = self.compute_next_pos_char(init_point, pointList, onChars=C_on, offChars=C_off, continueChar=True, terminate_wall=False, check_for_intersect=False)

        pt_a_upd = C_on[-1][-1]
        [pt4_dwn, pt4_ups, def4, beta4, delta_thet_w, pt3p, shockObj] = shock.to_wall_shock_point(pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a_upd, y_x, dydx, self.pcTOL, self.delta, self.gasProps, shockDir)
        
        return [pt4_dwn, pt4_ups, def4, beta4, delta_thet_w, pt3p, pt_a_upd, shockObj, C_on, C_off]

    def compile_mesh_points(self):
        """
        dumps all mesh points into one bucket, assigns them all indices, and makes a new triangle list of point indices
        need to run this before plotting the mesh
        Useful for debugging
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
        for tri in self.triangle_obj: #!TEMPORAY IMPROVE LATER
            self.triangle.append([pt.i if pt is not None else None for pt in tri])

        if self.init_method in ["STRAIGHT IDL", "MACH LINE"]:
            self.p0_ratio_by_region = [self.idl_p0_p0f]

        if self.shockMesh: 
            
            self.shock_segs = []
            [self.shock_segs.append(pt.i) for pt in self.shock_upst_segments_obj if pt.i is not None]
            
            #compute total pressure of each shock region
            if len(self.shock_object_list) > 0: #if at least one shock 
                for shock in self.shock_object_list:
                    if len(shock) == 0: continue
                    tot_press_loss = [s.p02_p01 for s in shock]
                    avg_tot_press_loss = sum(tot_press_loss)/len(shock)
                    p0_p0f = avg_tot_press_loss
                    for ratio in self.p0_ratio_by_region: 
                        p0_p0f *= ratio
                    self.p0_ratio_by_region.append(p0_p0f)

    def get_point_properties_and_minmax(self):
        """
        """
        self.minmax_p_p0f = [0,0]
        self.minmax_T_T0 = [0,0]
        self.minmax_rho_rho0f = [0,0]
        self.minmax_V = [0,0]
        for pt in self.meshPts:
            pt.get_point_properties(self)

            if pt.p_p0f < self.minmax_p_p0f[0]: self.minmax_p_p0f[0] = pt.p_p0f
            if pt.p_p0f > self.minmax_p_p0f[-1]: self.minmax_p_p0f[-1] = pt.p_p0f

            if pt.T_T0 < self.minmax_T_T0[0]: self.minmax_T_T0[0] = pt.T_T0
            if pt.T_T0 > self.minmax_T_T0[-1]: self.minmax_T_T0[-1] = pt.T_T0

            if pt.rho_rho0f < self.minmax_rho_rho0f[0]: 
                self.minmax_rho_rho0f[0] = pt.rho_rho0f
            if pt.rho_rho0f > self.minmax_rho_rho0f[-1]: 
                self.minmax_rho_rho0f[-1] = pt.rho_rho0f

            V = math.sqrt(pt.u**2 + pt.v**2)
            if V < self.minmax_V[0]: self.minmax_V[0] = V
            if V > self.minmax_V[-1]: self.minmax_V[-1] = V



    def find_mesh_point(self, pt, C_posneg):
        """
        gets the index of point in C_posneg (self.C_pos or self.C_neg)
        returns [i,j]
            i = line index 
            j = point index
        """
        for i,char in enumerate(C_posneg):
            for j,p in enumerate(char): 
                if p == pt: return [i,j]

        return [None, None]

    def get_max_indices(self, C_posneg):
        return len(C_posneg), len(C_posneg[-1])

    def find_point_index(self, pt, C_posneg, init=None):
        """
        Search through a 2D list of points and return the indices of a given point.
        Works for non-rectangular lists with i_init,j_init being any location. 

        Args:
            pt (object): The point object to be found.
            C_posneg (list): The 2D list of points to search through.
            init (list): [i_init, j_init]
                i_init (int): The starting index for the first dimension of the search range.
                j_init (int): The starting index for the second dimension of the search range.

        Returns:
            tuple: The indices of the given point in the list. Returns None if the point is not found.
        """

        if init==None: 
            i_init, j_init = len(C_posneg)-1, len(C_posneg[-1])-1
        else: 
            i_init, j_init = init[0], init[1]

        search_range = [(i_init, j_init)]
    
        while search_range:
            i, j = search_range.pop(0)
    
            if C_posneg[i][j] == pt:
                return i, j
    
            # Expand search range
            for new_i in range(i-1, i+2):
                for new_j in range(j-1, j+2):
                    if (0 <= new_i < len(C_posneg)
                        and 0 <= new_j < len(C_posneg[0])
                        and C_posneg[new_i][new_j] is not None
                        and (new_i, new_j) not in search_range):
                        search_range.append((new_i, new_j))
    
        # Point not found
        return None        

    def trim_mesh_after_intersect(self, pt2, pt1, charDir):
        """
        deletes a portion of the mesh which exists downstream of a crossed characteristic
        pt2,pt1 are parent points to the hypothetic pt3 which exhibited the intersection 
            pt2 = upper point 
            pt1 = lower point 
        charDir is the family of the crossed characteristic ("pos" or "neg")
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

    def delete_mesh_points(self, delPts, C_on=None, C_off=None):
        """
        given a list of mesh points, this function will delete them from the mesh, line segments included
        """
        def delete_from_charList(C1, C2):
            for i,char in enumerate(C1):
                newchar = [pt for pt in char if pt not in delPts]
                C1[i] = newchar        

            for i,char in enumerate(C2):
                newchar = [pt for pt in char if pt not in delPts]
                C2[i] = newchar

        #delete from pos & neg char list 
        if C_on is None and C_off is None: 
            delete_from_charList(self.C_pos, self.C_neg)
        else: 
            delete_from_charList(C_on, C_off)

        #delete from triangle list
        removeTriInd = []
        for i,tri in enumerate(self.triangle_obj): 
            for pt in tri:
                if pt in delPts and pt is not None: 
                    removeTriInd.append(i)
                    break
        self.triangle_obj = [tri for i,tri in enumerate(self.triangle_obj) if i not in removeTriInd]

        #get rid of empty lists in both characteristics lists 
        emptyLists = [i for i,char in enumerate(self.C_pos) if len(char) == 0]
        self.C_pos = [char for i,char in enumerate(self.C_pos) if i not in emptyLists]
        emptyLists = [i for i,char in enumerate(self.C_neg) if len(char) == 0]
        self.C_neg = [char for i,char in enumerate(self.C_neg) if i not in emptyLists]

    def compile_wall_points(self, y_upper, y_lower): 
        """
        stores upper and lower wall mesh points
        """
        self.wallPtsUpper = []
        self.wallPtsLower = []

        for pt in self.meshPts:
            if pt.isWall:
                if y_upper(pt.x) is not None:
                    if abs(y_upper(pt.x) - pt.y) < 1e-8:
                        self.wallPtsUpper.append(pt)
                if y_lower(pt.x) is not None: 
                    if abs(y_lower(pt.x) - pt.y) < 1e-8:
                        self.wallPtsLower.append(pt)  

    def compute_local_mass_flow(self):
        
        def calculate_mass_flow(dl):
            """
            calculates the ratio of mass flow to freestream total density across a data line
            dl: list of mesh points dl[0] and dl[-1] should be wall points for this calculation to be meaningful 
            TODO: need static density at each point
            """
            #first, delete sequences of points which have same coordinate (i.e. pairs of shock points)
            dl_corrected = [dl[0]]
            for i,pt2 in enumerate(dl[1:]):
                pt1 = dl[i-1] 
                if np.linalg.norm(np.array([pt2.y - pt1.y, -1*(pt2.x - pt1.x)]), ord=2) != 0:
                    dl_corrected.append(pt2) 
            
            #iterate through corrected data line, and summate mass flow rate
            mdot = 0
            for i,pt2 in enumerate(dl_corrected):
                if i == 0: continue 
                pt1 = dl_corrected[i-1]

                V = np.array([0.5*(pt2.u + pt1.u), 0.5*(pt2.v + pt1.v)]) #average velocity 
                rho_avg = 0.5*(pt1.rho_rho0f + pt2.rho_rho0f)
                nHat = np.array([pt2.y - pt1.y, -1*(pt2.x - pt1.x)])
                nHat = np.divide(nHat, np.linalg.norm(nHat, ord=2))

                if self.delta == 1:
                    A = math.pi*(pt1.y + pt2.y)*math.sqrt((pt1.x - pt2.x)**2 + (pt1.y - pt2.y)**2)

                elif self.delta == 0: 
                    A = math.sqrt((pt1.x - pt2.x)**2 + (pt1.y - pt2.y)**2)

                mdot += abs(np.dot(np.multiply(rho_avg, V), np.multiply(nHat, A)))

            return mdot
        
        if self.idl[0].isWall == False or self.idl[-1].isWall == False: 
            raise ValueError("Initial Data Line Must Span Wall-to-Wall!")
        
        idl_mass_flow = calculate_mass_flow(self.idl)

        self.mesh_mass_flow = [[[],[]],[[],[]]] #[[[x's],[m's]],[[x's],[m's]]]

        for i,char in enumerate(self.C_pos):
            if len(char) <= 2: continue
            if char[0].isWall and char[-1].isWall:
                mflow = calculate_mass_flow(char)
                self.mesh_mass_flow[0][1].append(mflow/idl_mass_flow)
                self.mesh_mass_flow[0][0].append(i)

        for i,char in enumerate(self.C_neg):
            if len(char) < 2: continue
            if char[0].isWall and char[-1].isWall:
                mflow = calculate_mass_flow(char)
                self.mesh_mass_flow[1][1].append(mflow/idl_mass_flow)     
                self.mesh_mass_flow[1][0].append(i)       
            

class Mesh_Point: 
    """
    generates and manipulates individual mesh point objects
    """
    def __init__(self,x,y,u,v, region ,ind=None,isWall=False, isIdl=False, isShock=False):
        """
        x,y,u,v: position and velocity components
        ind: numerical index for point
        isWall: set to True if point exists on the wall boundary 
        isIdl: set to True if point belongs to the initial data line
        """
        self.x,self.y,self.u,self.v = x,y,u,v 
        self.i = ind
        self.reg = region
        self.isWall = isWall #is the point on the boundary? 
        self.isIdl = isIdl #is the point on the initial data line? 
        self.isShock = isShock #is the point on a shock? 

    def get_point_properties(self, mesh): 
        """
        Gets flow properties at an individual mesh point. Calculates temperature, pressure, density, mach number, etc. 
        """
        #unpacking
        gam, a0, T0 = mesh.gasProps.gam, mesh.gasProps.a0, mesh.gasProps.T0
        #p0 = mesh.tot_press_by_region[self.reg]
        V = math.sqrt(self.u**2 + self.v**2)
        a = math.sqrt(a0**2 - 0.5*(gam-1)*V**2)
        self.mach = V/a #mach number 
        self.T = T0/(1+0.5*(gam-1)*(V/a)**2) #static temperature
        self.T_T0 = (1 + 0.5*(gam-1)*self.mach**2)**-1
        #self.p = p0/((1 + 0.5*(gam-1)*(V/a)**2)**(gam/(gam-1))) #static pressure 
        self.p_p0 = ((1 + 0.5*(gam-1)*self.mach**2)**(gam/(gam-1)))**-1
        #self.rho = self.p/(mesh.gasProps.R*self.T) #ideal gas law
        self.rho_rho0 = ((1 + 0.5*(gam-1)*self.mach**2)**(1/(gam-1)))**-1

        self.p_p0f = self.p_p0*mesh.p0_ratio_by_region[self.reg]
        self.rho_rho0f = self.rho_rho0*mesh.p0_ratio_by_region[self.reg] #p and rho are directly proportional via ideal gas law 

def linear_interpolate(x, z1, z3, x1, x3):
    """
    simple 1-d linear interpolation between x1 and x3 with interpolated value between z1 and z3
    """
    return ((z3-z1)/(x3-x1))*(x-x1) + z1 #linear interpolation function