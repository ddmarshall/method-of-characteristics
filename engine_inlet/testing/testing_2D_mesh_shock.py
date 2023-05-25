import sys 
import os
import unittest
sys.path.append(os.getcwd())
import method_of_characteristics.oblique_shock as obs

class Test_2D_Mesh_Shock(unittest.TestCase):

    def test_2D_uniform_shock(self):
        """
        Tests in-mesh shock wave in uniform flow through a duct with a single 
        compression corner. Compares results to classic oblique shock results 
        """
