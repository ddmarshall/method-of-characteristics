import unittest 
import unit_processes as unit_processes 
"""
unit testing script to verify output of various method of characteristics operators is within expected bounds
"""
class Test_Operators(unittest.TestCase):
    
    #Testing Interior Operator 
    def test_interior(self): 

        #create input objects
        class arbitrary(): 
            def __init__(self): return

        pt1, pt2, gasProps = arbitrary(), arbitrary(), arbitrary()

        #adding appropriate object attributes
        pt1.x, pt1.y, pt1.u, pt1.v = 0.131460, 0.040118, 2473.4, 812.8
        pt2.x, pt2.y, pt2.u, pt2.v = 0.135683, 0.037123, 2502.8, 737.6
        gasProps.R, gasProps.gamma, gasProps.T_0 = 320, 1.2, 3000

        res = unit_processes.interior(pt1, pt2, gasProps, 0.05) #running interior operator on single test point
        expected_res = [2510.1, 780.2, 0.14118, 0.04056] 
        self.assertAlmostEqual(res[0], expected_res[0], places=1) #u4
        self.assertAlmostEqual(res[1], expected_res[1], places=1) #v4
        self.assertAlmostEqual(res[2], expected_res[2], places=4) #x4
        self.assertAlmostEqual(res[3], expected_res[3], places=4) #y4 

    #Testing Direct Wall Operator
    def test_dirWall(self):

        class Arbitrary():
            def __init__(self): return

        wall = lambda x : (22.1852 + 0.71568*(x*1000) - 0.0010787*(x*1000)**2)/1000 #wall y function 
        wall_dydx = lambda x : 0.071568*1000 - 2*1000*0.0010787*x #wall slope function 
        wall_xBounds = (0, 0.25) #wall minimum and maximum x (m)

        pt2, gasProps = Arbitrary(), Arbitrary()
        pt2.x, pt2.y, pt2.u, pt2.v = 0.135683, 0.037123, 2502.8, 737.6
        gasProps.R, gasProps.gamma, gasProps.T_0 = 320, 1.2, 3000

        res = unit_processes.dir_wall_above(pt2, wall, wall_xBounds, wall_dydx, gasProps, 0.05)
        expected_res = [0.063485, 0.063273, 1977.4, 1144.4]

        self.assertAlmostEqual(res[0], expected_res[0], places=4)
        self.assertAlmostEqual(res[1], expected_res[1], places=4)
        self.assertAlmostEqual(res[2], expected_res[2], places=1)
        self.assertAlmostEqual(res[3], expected_res[3], places=1)

    #Testing Inverse Wall Operator 
    def test_invWall(self):

        class Arbitrary(): 
            def __init__(self): return

        pt1, pt3, pt4, gasProps = Arbitrary(), Arbitrary(), Arbitrary(), Arbitrary()
        pt1.x, pt1.y, pt1.u, pt1.v = 0.005496, 0.02602, 1578.3, 705.7
        pt3.x, pt3.y, pt3.u, pt3.v = 0.005085, 0.02608, 1577.5, 702.3
        pt4.x, pt4.y = 0.005283, 0.02617
        gasProps.R, gasProps.gamma, gasProps.T_0 = 320, 1.2, 3000

        #testing above inverse wall operator: 
        res = unit_processes.inv_wall_above(pt1, pt3, pt4, gasProps, 0.05)
        expected_res = [1583.4, 738.3]
        self.assertAlmostEqual(res[0], expected_res[0], places=1)
        self.assertAlmostEqual(res[1], expected_res[1], places=1) 
    
    #Testing Symmetry Point Operator 
    def test_sym(self):

        #create input objects
        class arbitrary(): 
            def __init__(self): return

        pt1, gasProps = arbitrary(), arbitrary()

        #adding appropriate object attributes
        pt1.x, pt1.y, pt1.u, pt1.v = 0.079625, 0.001290, 2306.1, 35.7
        gasProps.R, gasProps.gamma, gasProps.T_0 = 320, 1.2, 3000

        res = unit_processes.axis(pt1, gasProps, 0.05) #running interior operator on single test point
        expected_res = [0.083308, 2332.4] 
        self.assertAlmostEqual(res[0], expected_res[0], places=4) #x4
        self.assertAlmostEqual(res[1], expected_res[1], places=1) #u4
    
if __name__ == "__main__":
    unittest.main() 