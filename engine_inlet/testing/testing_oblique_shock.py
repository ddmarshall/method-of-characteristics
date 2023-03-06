import sys 
import os
import unittest
import math
sys.path.append(os.getcwd()) #add path to taylor maccoll module
import method_of_characteristics.oblique_shock as obs

class Test_Oblique_Shock(unittest.TestCase):



    def test_solve_weak_shock(self):

        M1 = 3 
        thet = math.radians(10)
        gam = 1.4
        shock = obs.Oblique_Shock(M1, gam, thet=thet) #create object 
        beta_w_exp, beta_s_exp = math.radians(27.3826906), math.radians(86.4082502)
        p2_p1_exp, rho2_rho1_exp, T2_T1_exp, p02_p01_exp = 2.05447215,1.65458799,1.24168201,0.96308338

        self.assertAlmostEqual(shock.beta_w, beta_w_exp,        places=4)
        self.assertAlmostEqual(shock.beta_s, beta_s_exp,        places=4)
        self.assertAlmostEqual(shock.p2_p1, p2_p1_exp,          places=4)
        self.assertAlmostEqual(shock.rho2_rho1, rho2_rho1_exp,  places=4)
        self.assertAlmostEqual(shock.T2_T1, T2_T1_exp,          places=4)
        self.assertAlmostEqual(shock.p02_p01, p02_p01_exp,      places=4)



    """
    def test_get_flow_deflection(self):

        beta = math.radians(37.2101360)
        M = 3
        gam = 1.4
        thet = obs.get_flow_deflection(beta, M, gam)
        thet_exp = math.radians(8)

        self.assertAlmostEqual(thet, thet_exp, places=1)
    """



if __name__ == "__main__":
    unittest.main()
