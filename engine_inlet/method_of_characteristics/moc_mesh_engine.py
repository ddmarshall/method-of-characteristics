import method_of_characteristics.unit_processes as moc_op
import method_of_characteristics.oblique_shock_point_irrot as shock
import math
import numpy as np
"""
Module responsible for generating method-of-characteristics mesh and mesh points
TODO: for shock mesh, make function to perform inverse wall operations in compression corners and add to shockPoint list
"""

class Mesh:

    def __init__(self, idl, Geom, gasProps, delta, pcTOL, kill_func, explicit_shocks=False):
        self.C_pos, self.C_neg = [],[] #containers for characteristics lines
        self.triangle_obj = [] #container for point line segments
        self.shock_upst_segments_obj = [] #container for shock line segments
        self.shockPts_frontside = [] 
        self.shockPts_backside = []
        self.funcs = moc_op.operator_funcs() #operator functions
        self.gasProps = gasProps
        self.delta = delta
        self.pcTOL = pcTOL
        self.geom = Geom
        self.f_kill = [kill_func, False]
        self.alternateChar = False
        self.numPointsGen = 0
        self.hasIntersected = False
        self.shockMesh = explicit_shocks

        self.idl = []
        for i,x in enumerate(idl.x):
            self.idl.append(Mesh_Point(x, idl.y[i], idl.u[i], idl.v[i], None, isIdl=True))
            #TODO add check for endpoints on the wall and set booleans to true
        if hasattr(idl, 'cowlPoint'):
            self.idl[0].isWall = True 
        if hasattr(idl, 'cbPoint'):
            self.idl[-1].isWall = True

        if explicit_shocks: 
            #adding first shock point
            #shockPt_init = Mesh_Point(self.idl[0].x,self.idl[0].y, self.idl[0].u, self.idl[0].v, isWall=True, isIdl=True, isShock=True)
            shockPt_init = self.idl[0]
            shockPt_init.isShock = True
            #!HARD CODED (FORCING TOP POINT IN IDL TO BE SHOCK ORIGIN)
            self.shockPts_frontside.append(shockPt_init)
            self.shock_upst_segments_obj.append(shockPt_init)
            self.generate_mesh_with_shocks() #generate mesh with shocks

        else: 
            self.generate_mesh()
        
        self.compile_mesh_points()
        [pt.get_point_properties(self.gasProps) for pt in self.meshPts]
        self.compile_wall_points(self.geom.y_cowl, self.geom.y_centerbody)

    def generate_mesh(self):
        """
        generates shock-less characteristic mesh until the kill function is triggered 
        Handles shock waves implicitly i.e. trims mesh after same-family intersection and does NOT include oblique shock calculations
        """
        charDir = "neg" #!hard coded for now... starting direction
        self.generate_initial_mesh_from_idl(charDir) #generate mesh from initial data line
        self.compile_mesh_points()
        while self.f_kill[1] == False: 
            #!following sections were written for initiating with neg lines. Untested for starting with positive
            if charDir == "neg" and self.alternateChar == False:
                #TODO use generate_mesh_from_line here instead 
                
                #generate initial above-wall point 
                pt1 = self.C_neg[-1][1] #2nd point in most recent characteristic
                [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
                pt3 = Mesh_Point(x3,y3,u3,v3, isWall=True)
                self.numPointsGen += 1 
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

    def generate_mesh_with_shocks(self):
        """
        generates characteristic mesh using oblique shock points
        """
        charDir = "neg"
        self.generate_initial_mesh_from_idl(charDir)
        
        if charDir == "neg": charList = self.C_neg
        elif charDir == "pos": charList = self.C_pos
        if self.check_for_passed_shock_point(charList) == True:
            self.compute_wall_to_wall_shock(charDir)

        
        while self.f_kill[1] == False:

            break 

    def generate_initial_mesh_from_idl(self, charDir):
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

    def generate_mesh_from_line(self, dl, charDir):
        """
        Generates the characteristic mesh from an existic characteristic (either pos or neg) terminates when all opposite characteristics originating 
        from dl have been generated
        dl = list of mesh points from positive/negative characteristic spanning wall to wall 
        """
             
        if charDir == "pos":
            #top wall solution
            pt1 = dl[1]
            [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
            self.numPointsGen += 1
            if self.f_kill[0](self) == True:
                self.f_kill[1] = True
                return 
            pt3 = Mesh_Point(x3,y3,u3,v3,isWall=True)
            self.triangle_obj.append([pt3, None, pt1])
            self.C_neg.append([pt3])
            i,_ = self.find_mesh_point(pt1, self.C_pos)
            self.C_pos[i].append(pt3)
        
        elif charDir == "neg":
            #bottom wall solution
            pt2 = dl[1]
            [x3,y3,u3,v3] = moc_op.direct_wall(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs, "neg")
            self.numPointsGen += 1
            if self.f_kill[0](self) == True:
                self.f_kill[1] = True
                return 
            pt3 = Mesh_Point(x3,y3,u3,v3,isWall=True)
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

    def compute_next_neg_char(self, init_point, prev_n_char, continueChar=False, terminate_wall=True, check_for_intersect=True):
        """
        Generates the next leading negative characteristic by advancing the mesh along the previous one
        init_point could be wall or idl point
        !NOTE TO SELF: CHANGES HERE NEED TO BE REFLECTED IN TWIN FUNCTION
        """
           
        if continueChar is False:
            self.C_neg.append([init_point])
        for pt in prev_n_char:
            pt2 = init_point #above
            pt1 = pt #below 
            [x3, y3, u3, v3] = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.pcTOL, self.funcs)
            self.numPointsGen += 1
            pt3 = Mesh_Point(x3, y3, u3, v3)
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
                i,j = self.find_mesh_point(pt2, self.C_neg)
                self.C_neg[i].append(pt3)
            else: 
                self.C_neg[-1].append(pt3)

            self.triangle_obj.append([pt3, pt2, pt1])

            #find point 1 and append the new point 
            i,_ = self.find_mesh_point(pt1, self.C_pos)
            self.C_pos[i].append(pt3)

            init_point = pt3 #update initial point

        #terminating wall point
        if terminate_wall: 
            if continueChar:
                pt2 = self.C_pos[-1][-1]
                [x3,y3,u3,v3] = moc_op.direct_wall(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs, "neg")
                self.numPointsGen += 1
                pt3 = Mesh_Point(x3, y3, u3, v3, isWall=True)
                i,_ = self.find_mesh_point(pt2, self.C_neg)
                self.C_neg[i].append(pt3)

            else:
                pt2 = self.C_neg[-1][-1]
                [x3,y3,u3,v3] = moc_op.direct_wall(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs, "neg")
                self.numPointsGen += 1
                pt3 = Mesh_Point(x3, y3, u3, v3, isWall=True)
                self.C_neg[-1].append(pt3)

            self.C_pos.append([pt3])
            self.triangle_obj.append([pt3, pt2, None])

    def compute_next_pos_char(self, init_point, prev_p_char, continueChar=False, terminate_wall=True, check_for_intersect=True):
        """
        Generates the next leading positive characteristic by advancing the mesh along the previous one
        init_point could be wall or idl point
        """
        if continueChar is False:
            self.C_pos.append([init_point])
        for pt in prev_p_char:
            pt2 = pt #above
            pt1 = init_point #below 
            [x3, y3, u3, v3] = moc_op.interior_point(pt1, pt2, self.gasProps, self.delta, self.pcTOL, self.funcs)
            self.numPointsGen += 1
            pt3 = Mesh_Point(x3, y3, u3, v3)
            if check_for_intersect: 
                if self.check_for_int_intersect(pt3, pt2, pt1, "pos"):
                    self.hasIntersected = True 
                    self.trim_mesh_after_intersect(pt2, pt1, "pos")
                    init_point = pt1
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

        #terminating wall point
        if terminate_wall:
            if continueChar:
                pt1 = self.C_neg[-1][-1]
                [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
                self.numPointsGen += 1
                pt3 = Mesh_Point(x3, y3, u3, v3, isWall=True)
                i,_ = self.find_mesh_point(pt1, self.C_pos)
                self.C_pos[i].append(pt3)

            else: 
                pt1 = self.C_pos[-1][-1]
                [x3,y3,u3,v3] = moc_op.direct_wall(pt1, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos")
                self.numPointsGen += 1
                pt3 = Mesh_Point(x3, y3, u3, v3, isWall=True)
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
                pt3 = Mesh_Point(x3,y3,u3,v3, isWall=True)
                self.numPointsGen += 1 
                self.triangle_obj.append([pt3, None, pt1]) 
                i,_ = self.find_mesh_point(pt1, self.C_pos)
                self.C_pos[i].append(pt3) #add to positive
                init_point = pt3
                prev_n_char = self.C_neg[-1][2:]
                self.compute_next_neg_char(init_point, prev_n_char)

            elif charDir == "pos":
                
                pt2 = self.C_pos[-1][1] #2nd point in most recent characteristic
                [x3,y3,u3,v3] = moc_op.direct_wall(pt2, self.geom.y_centerbody, self.geom.dydx_centerbody, self.gasProps, self.delta, self.pcTOL, self.funcs, "neg")
                pt3 = Mesh_Point(x3,y3,u3,v3, isWall=True)
                self.numPointsGen += 1 
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

    def compute_wall_to_wall_shock(self, shockDir):
        """
        computes a shock wave from one wall to another, deleting and appending mesh points at it moves along
        !Following function only works for negative shock direction 
        TODO need a way to handle shock intersecting same family characteristic (upstream of shock)
        """
        #generate from-wall shock point
        pt_w_ups = self.shockPts_frontside[-1]
        [i,j] = self.find_mesh_point(pt_w_ups, self.C_neg)
        pt = self.C_neg[i][j+1]
        [ii,jj] = self.find_mesh_point(pt, self.C_pos)
        pt1 = self.C_pos[ii][jj-1]

        y_x, dydx = self.geom.y_cowl, self.geom.dydx_cowl
        
        [pt4_dwn, pt4_ups, def_4, beta4, ptw_dwn, pt3p] = shock.wall_shock_point(pt_w_ups, y_x, dydx, pt1, self.pcTOL, self.delta, self.gasProps, shockDir)        
        pt4_dwn = Mesh_Point(pt4_dwn.x, pt4_dwn.y, pt4_dwn.u, pt4_dwn.v, isShock=True)
        pt4_ups = Mesh_Point(pt4_ups.x, pt4_ups.y, pt4_ups.u, pt4_ups.v, isShock=True)
        ptw_dwn = Mesh_Point(ptw_dwn.x, ptw_dwn.y, ptw_dwn.u, ptw_dwn.v, isShock=True)
        pt3p = Mesh_Point(pt3p.x, pt3p.y, pt3p.u, pt3p.v, isWall=True)
        
        self.shock_upst_segments_obj.append(pt4_ups) #add shock pt4_ups to segment list 
        [self.shockPts_backside.append(pt) for pt in [ptw_dwn, pt4_dwn]]
        self.shockPts_frontside.append(pt4_ups)

        self.triangle_obj.append([pt4_ups, pt1, pt_w_ups]) #add mesh segments to triangle list 
        self.triangle_obj.append([pt3p, pt4_ups, None]) 

        self.C_pos[ii][jj] = pt4_ups #replace old point with shock point
        self.C_pos[ii].append(pt3p) #add new point  
        self.C_neg[i][j+1] = pt4_ups #replace old point with shock point
        self.C_neg.append([pt3p])
        
        #generate interior shock points 
        upcoming_wall = False
        pt_s_ups, pt_s_dwn = pt4_ups, pt4_dwn 
        while upcoming_wall == False: 
            [i,j] = self.find_mesh_point(pt_s_ups, self.C_neg)
            pt = self.C_neg[i][j+1]
            [ii,jj] = self.find_mesh_point(pt, self.C_pos)
            pt1_old = pt1
            pt1 = self.C_pos[ii][jj-1]

            #check for and handle intersection of shock wave with upstream (same-family) characteristic
            delPts = []
            if self.check_for_shock_char_intersect(pt1_old, pt1, pt_s_ups, waveAng=beta4):
                pt1, delPts = self.get_upstream_point_for_shock(self.C_pos[ii], pt_s_ups, beta4) #finds new pt1
                if pt1 is None: 
                    pt1 = pt1_old
                    break #will be none if the next place for shock to hit is the wall 
                [iii, jjj] = self.find_mesh_point(pt1_old, self.C_neg)
                [delPts.append(pt) for pt in self.C_neg[iii][jjj+1:]]

            pt_a = pt3p
            [pt4_dwn, pt4_ups, def_4, beta4, pt3p] = shock.interior_shock_point(pt_s_ups, pt_s_dwn, beta4, def_4, pt1, pt_a, self.pcTOL, self.delta, self.gasProps, shockDir)
            
            #check for and handle intersection of shock wave with downstream (same-family) characteristic
            if self.check_for_shock_char_intersect(pt3p, pt_a, pt_s_ups, shockPt2=pt4_ups):
                print("\tInterior shock segment and same-family mach line intersection detected! Applying fix...")
                [pt4_dwn, pt4_ups, def_4, beta4, pt3p, pt_a] = self.handle_shock_char_intersect(pt_s_ups, pt_s_dwn, beta4, def_4, pt1, pt_a, shockDir)
            
            print(f"\tshock wave angle {math.degrees(beta4)} deg")
            pt4_dwn = Mesh_Point(pt4_dwn.x, pt4_dwn.y, pt4_dwn.u, pt4_dwn.v, isShock=True)
            pt4_ups = Mesh_Point(pt4_ups.x, pt4_ups.y, pt4_ups.u, pt4_ups.v, isShock=True)
            pt3p = Mesh_Point(pt3p.x, pt3p.y, pt3p.u, pt3p.v)
            
            self.shock_upst_segments_obj.append(pt4_ups)
            self.shockPts_backside.append(pt4_dwn)
            self.shockPts_frontside.append(pt4_ups)

            self.triangle_obj.append([pt4_ups, pt1, pt_s_ups])
            self.triangle_obj.append([pt3p, pt4_ups, pt_a])

            self.C_pos[ii][jj] = pt4_ups #replace old point with shock point
            self.C_pos[ii].append(pt3p) #add new point  
            self.C_neg[i][j+1] = pt4_ups #replace old point with shock point
            self.C_neg[-1].append(pt3p)     
            if len(delPts) > 0: self.delete_mesh_points(delPts)

            pt_s_ups, pt_s_dwn = pt4_ups, pt4_dwn

            #check if next point is on the wall 
            if self.C_pos[ii][jj-1].isWall:
                upcoming_wall = True   
             
        #setup for to-wall shock point
        [i,j] = self.find_mesh_point(pt_s_ups, self.C_neg)
        pt = self.C_neg[i][j+1]
        [ii,jj] = self.find_mesh_point(pt, self.C_pos)
        pt1_old = pt1
        pt1 = self.C_pos[ii][jj-1]

        #check for and handle intersection of shock wave with upstream (same-family) characteristic
        delPts = []
        if self.check_for_shock_char_intersect(pt1_old, pt1, pt_s_ups, waveAng=beta4):
            charPts = []
            for char in self.C_neg:
                if char[-1].isWall: 
                    charPts.append(char[-1])
            pt1, delPts = self.get_upstream_point_for_shock(charPts, pt_s_ups, beta4) #find a new pt1
            [iii, jjj] = self.find_mesh_point(pt1_old, self.C_neg)
            [delPts.append(pt) for pt in self.C_neg[iii][jjj+1:]]


        #generate to-wall shock point
        pt_a = pt3p
        y_x, dydx = self.geom.y_centerbody, self.geom.dydx_centerbody
        pt_s_ups, pt_s_dwn = pt4_ups, pt4_dwn
        [pt4_dwn, pt4_ups, def4_upd, beta4, delta_thet_w, pt3p] = shock.to_wall_shock_point(pt_s_ups, pt_s_dwn, beta4, def_4, pt1, pt_a, y_x, dydx, self.pcTOL, self.delta, self.gasProps, shockDir)
        #TODO check if new wall point causes shock to intersect characteristic 
        if self.check_for_shock_char_intersect(pt3p, pt_a, pt_s_ups, shockPt2=pt4_ups):
                print("\tWall shock segment and same-family mach line intersection detected! Applying fix...")

        pt4_dwn = Mesh_Point(pt4_dwn.x, pt4_dwn.y, pt4_dwn.u, pt4_dwn.v, isShock=True, isWall=True)
        pt4_ups = Mesh_Point(pt4_ups.x, pt4_ups.y, pt4_ups.u, pt4_ups.v, isShock=True, isWall=True)
        pt3p = Mesh_Point(pt3p.x, pt3p.y, pt3p.u, pt3p.v)

        self.shock_upst_segments_obj.append(pt4_ups)
        self.shockPts_backside.append(pt4_dwn)
        self.shockPts_frontside.append(pt4_ups)

        self.triangle_obj.append([pt4_ups, None, pt_s_ups])
        self.triangle_obj.append([pt3p, pt4_ups, pt_a])

        self.C_pos[ii][jj] = pt4_ups #replace old point with shock point
        self.C_pos[ii].append(pt3p) #add new point  
        self.C_neg[i][j+1] = pt4_ups #replace old point with shock point
        self.C_neg[-1].append(pt3p)
        if len(delPts) > 0: self.delete_mesh_points(delPts) #delete points whihch are no longer needed

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
 
        delPts = [pt for ind,pt in enumerate(charPts) if ind>i]
        return retPt, delPts

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
            
            return False

    def handle_shock_char_intersect(self, pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a, shockDir):
        """
        generates a downstream characteristic to handle intersection of a shock segment and same-family characteristic. To be called after check_for_shock_char_intersect returns
        true
        !Currently written for negative shock direction
        """
        if shockDir == "neg":
            i,j = self.find_mesh_point(pt_a, self.C_neg)
            pointList = self.C_neg[i][2:j+1]
            pt = self.C_neg[i][1] #first interior point
            [x3, y3, u3, v3] = moc_op.direct_wall(pt, self.geom.y_cowl, self.geom.dydx_cowl, self.gasProps, self.delta, self.pcTOL, self.funcs, "pos") #create first wall point
            init_point = Mesh_Point(x3, y3, u3, v3, isWall=True)
            self.C_neg.append([init_point]), self.triangle_obj.append([init_point, pt, None])
            self.compute_next_neg_char(init_point, pointList, terminate_wall=False, check_for_intersect=False)
            
            pt_a_upd = self.C_neg[-1][-1]
            [pt4_dwn, pt4_ups, def_4, beta4, pt3p] = shock.interior_shock_point(pt_s_ups, pt_s_dwn, beta_s, def_s, pt1, pt_a_upd, self.pcTOL, self.delta, self.gasProps, shockDir)
            return [pt4_dwn, pt4_ups, def_4, beta4, pt3p, pt_a_upd]
             
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

        if self.shockMesh: 
            self.shock_segs = []
            for pt in self.shock_upst_segments_obj:
                self.shock_segs.append(pt.i)
            
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

    def delete_mesh_points(self, delPts):
        """
        given a list of mesh points, this function will delete them from the mesh, line segments included
        """
        #delete from pos & neg char list 
        for i,char in enumerate(self.C_pos):
            newchar = [pt for pt in char if pt not in delPts]
            self.C_pos[i] = newchar        
        
        for i,char in enumerate(self.C_neg):
            newchar = [pt for pt in char if pt not in delPts]
            self.C_neg[i] = newchar   

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

    def calculate_mass_flow(self, dl, delta):
        """
        calculates the total mass flow rate across a data line
        dl: list of mesh points dl[0] and dl[-1] should be wall points for this calculation to be meaningful 
        TODO: need static density at each point
        """
        mdot = 0 
        for i,pt2 in enumerate(dl):
            if i == 0: continue 
            pt1 = dl[i-1]

            V = np.array([0.5*(pt2.u + pt1.u), 0.5*(pt2.v + pt1.v)]) #average velocity 
            rho_avg = 0.5*(pt1.rho + pt2.rho)
            nHat = np.array([pt2.y - pt1.y, -1*(pt2.x - pt1.x)])
            nHat = np.divide(nHat, np.linalg.norm(nHat, ord=2))

            if delta == 1:
                A = math.pi*(pt1.y + pt2.y)*math.sqrt((pt1.x - pt2.x)**2 + (pt1.y - pt2.y)**2)

            elif delta == 0: 
                A = math.sqrt((pt1.x - pt2.x)**2 + (pt1.y - pt2.y)**2)
            
            mdot += np.dot(np.multiply(rho_avg, V), np.multiply(nHat, A))
        
        return mdot

    def compile_wall_points(self, y_upper, y_lower): 
        """
        stores upper and lower wall mesh points
        """
        self.wallPtsUpper = []
        self.wallPtsLower = []

        for pt in self.meshPts:
            if pt.isWall: 
                if abs(y_upper(pt.x) - pt.y) < 1e-8:
                    self.wallPtsUpper.append(pt)
                elif abs(y_lower(pt.x) - pt.y) < 1e-8:
                    self.wallPtsLower.append(pt)  


class Mesh_Point: 
    """
    generates and manipulates individual mesh point objects
    """
    def __init__(self,x,y,u,v,ind=None,isWall=False, isIdl=False, isShock=False):
        """
        x,y,u,v: position and velocity components
        ind: numerical index for point
        isWall: set to True if point exists on the wall boundary 
        isIdl: set to True if point belongs to the initial data line
        """
        self.x,self.y,self.u,self.v = x,y,u,v 
        self.i = ind
        self.isWall = isWall #is the point on the boundary? 
        self.isIdl = isIdl #is the point on the initial data line? 
        self.isShock = isShock #is the point on a shock? 

    def get_point_properties(self, gasProps): 
        """
        Gets flow properties at an individual mesh point. Calculates temperature, pressure, density, mach number, etc. 
        """
        #unpacking
        gam, a0, T0, p0 = gasProps.gam, gasProps.a0, gasProps.T0, gasProps.p0
        V = math.sqrt(self.u**2 + self.v**2)
        a = math.sqrt(a0**2 + 0.5*(gam-1)*V**2)
        self.mach = V/a #mach number 
        self.T = T0/(1+0.5*(gam-1)*(V/a)**2) #static temperature
        self.p = p0/((1 + 0.5*(gam-1)*(V/a)**2)**(gam/(gam-1))) #static pressure 

        self.rho = self.p/(gasProps.R*self.T) #ideal gas law