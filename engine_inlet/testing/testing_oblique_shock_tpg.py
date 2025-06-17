import sys 
import os
import unittest
import math
import numpy as np
import matplotlib.pyplot as plt
#print(os.getcwd())
#sys.path.append("./method_of_characteristics")
#sys.path.append('C:/Users/emeve/Documents/GitHub/method-of-characteristics/engine_inlet/method_of_characteristics')
import engine_inlet.method_of_characteristics.oblique_shock as obs
import engine_inlet.method_of_characteristics.oblique_shock_tpg as obs_tpg
import engine_inlet.method_of_characteristics.gas_state_TPG as Gas_state
from gas_state import gas_state 

gas = Gas_state.ThermallyPerfectGas


class Test_Oblique_Shock(unittest.TestCase):

    def test_solve_weak_shock_using_known_beta(self):

        M1 = 1.7 
        theta = math.radians(10)
        gam = 1.4
        beta = math.radians(30)
        
        Rho = 1.225
        T = 40
        Tv = 3500 
        R = 287
        P = Rho * 287 * T
        delta = 1
        h0 = gas._h0_given_T_M(R, Tv, T, M1)

        T0 = gas._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = gas._a(self, T0)
                else: 
                    self.a0 = gas._a(self, T0)
        GS = gasProps(1.4, 287, T0, p0, delta, Tv)
        shock = obs.Oblique_Shock(M1, GS, beta=beta) #create object
        tpg_shock = obs_tpg.Oblique_Shock( M1, GS, beta=beta)

        #
        #beta_w_exp, beta_s_exp = math.radians(27.3826906), math.radians(86.4082502)
        #p2_p1_exp, rho2_rho1_exp, T2_T1_exp, p02_p01_exp = 2.05447215,1.65458799,1.24168201,0.96308338
        #
        self.assertAlmostEqual(shock.beta, tpg_shock.beta,        places=4)
        #self.assertAlmostEqual(shock.beta_s, beta_s_exp,        places=3)
        self.assertAlmostEqual(shock.p2_p1, tpg_shock.p2_p1,          places=4)
        self.assertAlmostEqual(shock.rho2_rho1, tpg_shock.rho2_rho1,  places=4)
        self.assertAlmostEqual(shock.T2_T1, tpg_shock.T2_T1,          places=4)
        #self.assertAlmostEqual(shock.p02_p01, p02_p01_exp,      places=4)
        self.assertAlmostEqual(shock.deflec, tpg_shock.deflec,        places=4)



    def test_solve_weak_shock_using_known_deflec(self):

            M1 = 1.7 
            deflec = math.radians(10)
            gam = 1.4
            beta = math.radians(30)
            R = 287
            Rho = 1.225
            T = 40
            Tv = 3500 
            P = Rho * 287 * T
            delta = 1

            h0 = gas._h0_given_T_M(R, Tv, T, M1)
            T0 = gas._T0_given_h0(R, Tv, h0, T)
            p0 = 101325 
            class gasProps:
                def __init__(self, gam, R, T0, p0, delta, Tv ): 
                    self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                    self.p0_p0f = 1
                    if delta == 1: 
                        #Enable Vibrational Mode
                        self.Tv = Tv
                        self.p0f = p0
                        self.a0 = gas._a(self, T0)
                    else: 
                        self.a0 = gas._a(self, T0)
            GS = gasProps(1.4, 287, T0, p0, delta, Tv)

            shock_deflec = obs.Oblique_Shock(M1, GS, deflec=deflec) #create object
            tpg_shock_deflec = obs_tpg.Oblique_Shock( M1, GS, deflec=deflec)

            #
            #beta_w_exp, beta_s_exp = math.radians(27.3826906), math.radians(86.4082502)
            #p2_p1_exp, rho2_rho1_exp, T2_T1_exp, p02_p01_exp = 2.05447215,1.65458799,1.24168201,0.96308338
            #
            self.assertAlmostEqual(shock_deflec.beta, tpg_shock_deflec.beta,        places=4)
            self.assertAlmostEqual(shock_deflec.p2_p1, tpg_shock_deflec.p2_p1,          places=4)
            self.assertAlmostEqual(shock_deflec.rho2_rho1, tpg_shock_deflec.rho2_rho1,  places=4)
            self.assertAlmostEqual(shock_deflec.T2_T1, tpg_shock_deflec.T2_T1,          places=4)
            self.assertAlmostEqual(shock_deflec.deflec, float(tpg_shock_deflec.deflec),        places=4)
            #self.assertAlmostEqual(shock.p02_p01, p02_p01_exp,      places=4)
    def test_high_temperature_internal_test(self):
        #Test accuracy at higher values by comparing the results between both solution schemes
        #using beta of 10 deg, solve shock solver for deflection angle
        #That deflection angle will be used to solve for beta in the other solver: vlaue should come back to 10 deg
            M1 = 4.5
            beta = math.radians(30)
            gam = 1.4
            Rho = 1.225
            T = 400
            Tv = 3500 
            R = 287
            P = Rho * 287 * T
            delta = 1
            
            h0 = gas._h0_given_T_M(R, Tv, T, M1)
            T0 = gas._T0_given_h0(R, Tv, h0, T)
            p0 = 101325 
            class gasProps:
                def __init__(self, gam, R, T0, p0, delta, Tv ): 
                    self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                    self.p0_p0f = 1
                    if delta == 1: 
                        #Enable Vibrational Mode
                        self.Tv = Tv
                        self.p0f = p0
                        self.a0 = gas._a(self, T0)
                    else: 
                        self.a0 = gas._a(self, T0)
            GS = gasProps(1.4, 287, T0, p0, delta, Tv)
            tpg_shock_beta = obs_tpg.Oblique_Shock( M1, GS, beta=beta)

            deflec = tpg_shock_beta.deflec 
            tpg_shock_deflec = obs_tpg.Oblique_Shock( M1, GS, deflec=deflec)

            self.assertAlmostEqual(tpg_shock_beta.beta, tpg_shock_deflec.beta, places=4)

    def test_external_NASA_data_deflec_M2_4_T390(self): 
        Tv = 3500
        R = 287 
        gam = 1.4
        delta = 1

        M = np.array([2.4, 2.6, 2.8])
        deflec = np.array([14.7165, 15.4762, 16.8338, 18.2636, 19.6948, 20.8362, 
                           21.9599, 23.211, 24.2635, 25.2072, 26.1337, 26.879, 
                           27.5535, 28.2122, 28.5284])
        deflec = np.deg2rad(deflec)
        Beta_Graph = np.array([37.9818, 38.7612, 40.3589, 41.9914, 43.9485, 45.5143,
                               47.1526, 49.0771, 50.8971, 52.5388, 54.3612, 55.9345, 
                               57.7616, 60.0942, 62.325])
        T = float(390 * 5 / 9) #Rankine conversion to K
        Rho = 1.225
        P = Rho * R * T

        h0 = gas._h0_given_T_M(R, Tv, T, M[0])
        T0 = gas._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = gas._a(self, T0)
                else: 
                    self.a0 = gas._a(self, T0)
        GS = gasProps(1.4, 287, T0, p0, delta, Tv)

        Beta = np.zeros(len(deflec))

        for i in range(len(deflec)):
            #shock = obs.Oblique_Shock(M[1], gam, deflec=deflec) #create object
            tpg_shock = obs_tpg.Oblique_Shock( M[0], GS, deflec=deflec[i])
            
            Beta[i] = tpg_shock.beta
            #print(i)
            #print(np.abs((np.rad2deg(Beta[i]) - Beta_Graph[i]) / Beta_Graph[i]) * 100)
    

    def test_external_NASA_data_beta_M2_4_T390(self): 
        Tv = 3500
        R = 287 
        gam = 1.4
        delta = 1

        M = np.array([2.4, 2.6, 2.8])
        deflec = np.array([14.7165, 15.4762, 16.8338, 18.2636, 19.6948, 20.8362, 
                           21.9599, 23.211, 24.2635, 25.2072, 26.1337, 26.879, 
                           27.5535, 28.2122, 28.5284])
        deflec = np.deg2rad(deflec)
        Beta_Graph = np.array([37.9818, 38.7612, 40.3589, 41.9914, 43.9485, 45.5143,
                               47.1526, 49.0771, 50.8971, 52.5388, 54.3612, 55.9345, 
                               57.7616, 60.0942, 62.325])
        T = float(390 * 5 / 9) #Rankine conversion to K
        Rho = 1.225
        P = Rho * R * T

        h0 = gas._h0_given_T_M(R, Tv, T, M[0])
        T0 = gas._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = gas._a(self, T0)
                else: 
                    self.a0 = gas._a(self, T0)
        GS = gasProps(1.4, 287, T0, p0, delta, Tv)

        Beta = np.zeros(len(deflec))

        for i in range(len(deflec)):
            #shock = obs.Oblique_Shock(M[1], gam, deflec=deflec) #create object
            if i > 16:
                with self.assertRaises(ValueError):
                    tpg_shock = obs_tpg.Oblique_Shock( M[0], GS, deflec=deflec[i])
            else: 
                tpg_shock = obs_tpg.Oblique_Shock( M[0], GS, deflec=deflec[i])
                Beta[i] = np.rad2deg(tpg_shock.beta)
                #print(i)
                #print(np.abs((Beta[i]) - Beta_Graph[i]) / Beta_Graph[i] * 100)
            
    def test_external_NASA_data_deflec_M7_T390(self): 
        Tv = 3500
        R = 287 
        gam = 1.4
        delta = 1

        M = np.array([2.4, 2.6, 2.8, 7])
        deflec = np.array([5.5619, 6.5195, 7.8935, 9.8917, 11.4745, 13.7986,
                           16.251, 19.0196, 20.7576, 23.047, 25.2191, 27.147, 
                           28.1975, 29.6098, 31.0585, 32.5796, 33.8926, 35.0886, 
                           36.2212, 37.3358, 38.5597, 40.0101, 41.2799, 41.951,
                           42.4422, 43.169, 43.9599, 44.7975, 45.5116, 45.8744])
        deflec = np.deg2rad(deflec)
        Beta_Graph = np.array([12.3409, 13.0814, 14.3361, 16.2465, 17.876, 20.1951,
                               23.0168, 26.085, 28.1807, 30.789, 33.454, 35.6359, 
                               37.005, 38.7641, 40.5946, 42.532, 44.1847, 45.9118, 
                               47.5319, 49.1343, 51.0232, 53.2505, 55.337, 56.4247, 
                               57.5338, 58.9992, 60.6979, 62.7745, 65.3946, 67.9672])
        T = float(390 * 5 / 9) #Rankine conversion to K
        Rho = 1.225
        P = Rho * R * T

        h0 = gas._h0_given_T_M(R, Tv, T, M[3])
        T0 = gas._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = gas._a(self, T0)
                else: 
                    self.a0 = gas._a(self, T0)
        GS = gasProps(1.4, 287, T0, p0, delta, Tv)

        Beta = np.zeros(len(deflec))

        for i in range(len(deflec)):
            #shock = obs.Oblique_Shock(M[1], gam, deflec=deflec) #create object

            if i > 25:
                with self.assertRaises(ValueError):
                    tpg_shock = obs_tpg.Oblique_Shock( M[3], GS, deflec=deflec[i])
            else: 
                tpg_shock = obs_tpg.Oblique_Shock( M[3], GS, deflec=deflec[i])
                Beta[i] = tpg_shock.beta
                print(i)
                print(np.abs((np.rad2deg(Beta[i]) - Beta_Graph[i]) / Beta_Graph[i]) * 100)
   
   
    def test_external_NASA_data_deflec_M2_4_T500(self): 
        
        Tv = 3500
        R = 287 
        gam = 1.4
        delta = 1

        M = np.array([2.4, 2.6, 2.8, 7])
        deflec = np.array([0.6577, 2.9853, 5.4454, 8.1053, 9.9028, 12.1622, 14.171, 15.993, 
                           17.5026, 20.1824, 22.0767, 23.7956, 25.3762, 26.8563, 27.8352, 28.4006])
        deflec = np.deg2rad(deflec)
        Beta_Graph = np.array([25.068, 26.846, 28.96, 31.251, 32.987, 35.176, 37.39, 39.277,
                               40.987, 44.709, 47.399, 50.114, 52.955, 55.896, 58.885, 61.899])
        T = float(500 * 5 / 9) #Rankine conversion to K
        Rho = 1.225
        P = Rho * R * T

        h0 = gas._h0_given_T_M(R, Tv, T, M[0])
        T0 = gas._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = gas._a(self, T0)
                else: 
                    self.a0 = gas._a(self, T0)
        GS = gasProps(1.4, 287, T0, p0, delta, Tv)

        Beta = np.zeros(len(deflec))

        for i in range(len(deflec)):
            #shock = obs.Oblique_Shock(M[1], gam, deflec=deflec) #create object

            if i > 16:
                with self.assertRaises(ValueError):
                    tpg_shock = obs_tpg.Oblique_Shock( M[0], GS, deflec=deflec[i])
            else: 
                tpg_shock = obs_tpg.Oblique_Shock( M[0], GS, deflec=deflec[i])
                Beta[i] = tpg_shock.beta
                #print(i)
                #print(np.abs((np.rad2deg(Beta[i]) - Beta_Graph[i]) / Beta_Graph[i]) * 100)
   

    def test_external_NASA_data_deflec_M7_T500(self): 
        
        Tv = 3500
        R = 287 
        gam = 1.4
        delta = 1

        M = np.array([2.4, 2.6, 2.8, 7])
        deflec = np.array([0.0551, 1.8815, 5.1428, 7.4657, 9.101, 11.51, 13.1825, 14.6673, 16.6766, 18.6227, 20.4571, 21.966, 
                           23.463, 26.2821, 28.3647, 30.8838, 32.9908, 34.4868, 35.8204, 37.6646, 
                           39.4336, 41.0398, 42.7204, 44.2125, 45.2302, 45.8102, 46.0259])
        deflec = np.deg2rad(deflec)
        Beta_Graph = np.array([8.238, 9.443, 12.012, 13.95, 15.536, 17.926, 19.586, 21.222, 23.285, 25.498,  
                               27.435, 29.321, 31.057, 34.352, 36.917, 40.061, 42.827, 44.814, 46.675, 49.39, 
                                52.105, 54.795, 57.711, 60.778, 63.441, 65.853, 67.46, 68.59])
        T = float(500 * 5 / 9) #Rankine conversion to K
        Rho = 1.225
        P = Rho * R * T

        h0 = gas._h0_given_T_M(R, Tv, T, M[3])
        T0 = gas._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = gas._a(self, T0)
                else: 
                    self.a0 = gas._a(self, T0)
        GS = gasProps(1.4, 287, T0, p0, delta, Tv)

        Beta = np.zeros(len(deflec))

        for i in range(len(deflec)):
            #shock = obs.Oblique_Shock(M[1], gam, deflec=deflec) #create object

            if i > 22:
                with self.assertRaises(ValueError):
                    tpg_shock = obs_tpg.Oblique_Shock( M[3], GS, deflec=deflec[i])
            else: 
                tpg_shock = obs_tpg.Oblique_Shock( M[3], GS, deflec=deflec[i])
                Beta[i] = tpg_shock.beta
                #print(i)
                #print(np.abs((np.rad2deg(Beta[i]) - Beta_Graph[i]) / Beta_Graph[i]) * 100)
   
    def test_external_NASA_data_deflec_M2_4_T630(self): 
        
        Tv = 3500
        R = 287 
        gam = 1.4
        delta = 1

        M = np.array([2.4, 2.6, 2.8, 7])
        deflec = np.array([0.3346, 3.3855, 7.2051, 10.0589, 12.4759, 15.6263, 18.7216, 20.661, 
                           22.9264, 25.1066, 26.6063, 27.9148, 28.8568 ])
        deflec = np.deg2rad(deflec)
        Beta_Graph = np.array([24.6, 27.052, 30.181, 32.721, 35.364, 38.759, 42.29, 
                               45.009, 48.062, 51.966, 55.558, 59.202, 62.429])
        T = float(630 * 5 / 9) #Rankine conversion to K
        Rho = 1.225
        P = Rho * R * T

        h0 = gas._h0_given_T_M(R, Tv, T, M[0])
        T0 = gas._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = gas._a(self, T0)
                else: 
                    self.a0 = gas._a(self, T0)
        GS = gasProps(1.4, 287, T0, p0, delta, Tv)

        Beta = np.zeros(len(deflec))

        for i in range(len(deflec)):
            #shock = obs.Oblique_Shock(M[1], gam, deflec=deflec) #create object

            if i > 11:
                with self.assertRaises(ValueError):
                    tpg_shock = obs_tpg.Oblique_Shock( M[0], GS, deflec=deflec[i])
            else: 
                tpg_shock = obs_tpg.Oblique_Shock( M[0], GS, deflec=deflec[i])
                Beta[i] = tpg_shock.beta
                #print(i)
                #print(np.abs((np.rad2deg(Beta[i]) - Beta_Graph[i]) / Beta_Graph[i]) * 100)
   
    def test_external_NASA_data_deflec_M7_T630(self): 
        
        Tv = 3500
        R = 287 
        gam = 1.4
        delta = 1

        M = np.array([2.4, 2.6, 2.8, 7])
        deflec = np.array([0.527, 3.7009, 6.0316, 8.2316, 10.6361, 13.4603, 16.6513, 18.6558, 20.7784, 
                           23.6229, 26.5641, 29.085, 34.7496, 37.6027, 39.4542, 40.7129, 42.9703, 44.6563, 45.6755, 
                           46.4676, 46.7436])
        deflec = np.deg2rad(deflec)
        Beta_Graph = np.array([8.554, 10.625, 12.467, 14.503, 16.826, 19.765, 23.014, 25.22, 27.693, 30.981, 
                               34.445, 37.467, 44.74, 48.617, 51.194, 53.152, 57.129, 60.532, 62.836, 
                               65.093, 66.47])
        T = float(630 * 5 / 9) #Rankine conversion to K
        Rho = 1.225
        P = Rho * R * T

        h0 = gas._h0_given_T_M(R, Tv, T, M[3])
        T0 = gas._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = gas._a(self, T0)
                else: 
                    self.a0 = gas._a(self, T0)
        GS = gasProps(1.4, 287, T0, p0, delta, Tv)

        Beta = np.zeros(len(deflec))

        for i in range(len(deflec)):
            #shock = obs.Oblique_Shock(M[1], gam, deflec=deflec) #create object

            if i > 16:
                with self.assertRaises(ValueError):
                    tpg_shock = obs_tpg.Oblique_Shock( M[3], GS, deflec=deflec[i])
            else: 
                tpg_shock = obs_tpg.Oblique_Shock( M[3], GS, deflec=deflec[i])
                Beta[i] = tpg_shock.beta
                print(i)
                print(np.abs((np.rad2deg(Beta[i]) - Beta_Graph[i]) / Beta_Graph[i]) * 100)
   
    def test_create_graph_T390(self): 
        Tv = 3500
        R = 287 
        gam = 1.4
        delta = 1

        M = np.array([2, 2.4, 3, 4, 5 , 7, 10])
        theta_graph = np.linspace(0.01,50,200)
        theta = np.deg2rad(theta_graph)
        T = float(390 * 5 / 9) #Rankine conversion to K
        Rho = 1.225
        P = Rho * R * T
        beta = np.zeros(len(theta))
        color = ['#3A913F', '#A4D65E', '#F2C75C', '#F8E08E', '#54585A', '#CAC7A7', '#B7CDC2', '#789F90']

        for i in range(len(M)):

            for j in range(len(theta)): 
                try:
                    h0 = gas._h0_given_T_M(R, Tv, T, M[i])
                    T0 = gas._T0_given_h0(R, Tv, h0, T)
                    p0 = 101325 
                    class gasProps:
                        def __init__(self, gam, R, T0, p0, delta, Tv ): 
                            self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                            self.p0_p0f = 1
                            if delta == 1: 
                                #Enable Vibrational Mode
                                self.Tv = Tv
                                self.p0f = p0
                                self.a0 = gas._a(self, T0)
                            else: 
                                self.a0 = gas._a(self, T0)
                    GS = gasProps(1.4, 287, T0, p0, delta, Tv) 
                    tpg_shock = obs_tpg.Oblique_Shock( M[i], GS, deflec=theta[j])
                    beta[j] = np.rad2deg(tpg_shock.beta)
                except ValueError:
                    break
            

            #exit interior loop, plot the line 
            mylabel = "M = " + str(M[i])
            plt.plot(theta_graph[0:j], beta[0:j], label=mylabel, color=color[i], linewidth=2.5)
        deflec_graph_M7 = np.array([5.5619, 6.5195, 7.8935, 9.8917, 11.4745, 13.7986,
                           16.251, 19.0196, 20.7576, 23.047, 25.2191, 27.147, 
                           28.1975, 29.6098, 31.0585, 32.5796, 33.8926, 35.0886, 
                           36.2212, 37.3358, 38.5597, 40.0101, 41.2799, 41.951,
                           42.4422, 43.169, 43.9599, 44.7975, 45.5116, 45.8744])

        Beta_Graph_M7 = np.array([12.3409, 13.0814, 14.3361, 16.2465, 17.876, 20.1951,
                               23.0168, 26.085, 28.1807, 30.789, 33.454, 35.6359, 
                               37.005, 38.7641, 40.5946, 42.532, 44.1847, 45.9118, 
                               47.5319, 49.1343, 51.0232, 53.2505, 55.337, 56.4247, 
                               57.5338, 58.9992, 60.6979, 62.7745, 65.3946, 67.9672])

        deflec_graph_M2_4 = np.array([14.7165, 15.4762, 16.8338, 18.2636, 19.6948, 20.8362, 
                           21.9599, 23.211, 24.2635, 25.2072, 26.1337, 26.879, 
                           27.5535, 28.2122, 28.5284])

        Beta_Graph_M2_4 = np.array([37.9818, 38.7612, 40.3589, 41.9914, 43.9485, 45.5143,
                               47.1526, 49.0771, 50.8971, 52.5388, 54.3612, 55.9345, 
                               57.7616, 60.0942, 62.325])        
        plt.plot(deflec_graph_M2_4, Beta_Graph_M2_4, 'x', markerfacecolor='none', markeredgecolor='black', markersize=10, label="NACA Data",)
        plt.plot(deflec_graph_M7,   Beta_Graph_M7,  'x', markerfacecolor='none', markeredgecolor='black', markersize=10)
        plt.rcParams['font.family'] = 'Times New Roman'
        plt.rcParams['font.size'] = 34
        plt.rcParams["mathtext.fontset"] = "custom"
        plt.rcParams["mathtext.rm"] = 'Times New Roman'
        plt.rcParams["mathtext.it"] = "Times New Roman:italic"
        plt.rcParams["mathtext.bf"] = "Times New Roman:bold"
        plt.rcParams["mathtext.default"] = "rm"
        plt.xlim(0, 48)
        plt.ylim(0,70)
        plt.xlabel('Deflection Angle, \u03B8 [\u00B0]', fontsize=34, fontname = "Times New Roman")
        plt.ylabel('Wave Angle, \u03B2 [\u00B0]', fontsize=34, fontname = "Times New Roman")
        plt.tick_params(axis='both', which='major', labelsize=34, labelfontfamily = "Times New Roman")
        plt.rcParams.update({'font.family': 'Times New Roman'})
        plt.grid(True)
        plt.legend()
        #show the plot    
        plt.show()





""""
# Sample data
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)
y3 = np.tan(x)

# Plot the lines
plt.plot(x, y1, label='Sine', color='blue', linestyle='-')
plt.plot(x, y2, label='Cosine', color='green', linestyle='--')
plt.plot(x, y3, label='Tangent', color='red', linestyle=':')

# Add titles and labels
plt.title('Multiple Line Plot')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')

# Show legend
plt.legend()

# Display the plot
plt.show()

"""
if __name__ == "__main__":


    unittest.main()