"""
This file runs unit tests on the Taylor_Maccoll_TPG.py module and other supporting functions 
"""
import sys 
import os
import unittest
import math
import numpy as np
import scipy
from scipy.optimize import fsolve
from scipy.optimize import least_squares
sys.path.append(os.getcwd()) #add path to taylor maccoll module
import engine_inlet.taylor_maccoll_cone.taylor_maccoll as taylor_maccoll_cpg
import engine_inlet.taylor_maccoll_cone.taylor_maccoll_tpg as taylor_maccoll_tpg
import engine_inlet.method_of_characteristics.gas_state_TPG as Gas_state 
from taylor_maccoll_tpg_Lampe_variant import TaylorMaccoll_Cone_TPG_Lampe as taylor_maccoll_Lampe
from  gas_state import gas_state

GS = Gas_state.ThermallyPerfectGas

class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = self.GS._a(self, T0)
                else: 
                    self.a0 = self.GS._a(self, T0)

class Test_Taylor_Maccoll(unittest.TestCase):
    def test_h_rootfinder(self): 
        R = 287
        Tv = 3500 
        T0 = 1000
        T_target = 500

        h0 = GS._h(R, Tv, T0)

        #Find the speed based on total enthalpy and target temperature
        V = np.sqrt(2 * (GS._h(R, Tv, T0) - GS._h(R, Tv, T_target))) #m/s 

        h = h0 - V**2 / 2

        def f(T):
            return h - GS._h(R, Tv, T)
        
        initial_guess = float(T0 / 2)
        bounds = (0, T0)
        T = least_squares(f, initial_guess, bounds=bounds)
        #print(T) 

    
    def test_Taylor_Maccoll_ODE_Solver_against_CPG(self): 
        #Build the gas state value holder 
        gam = 1.4       
        Rho = 1.225
        T = 133 + 1/3
        Tv = 3500 
        P = Rho * 287 * T
        delta = 1
        GS_init = gas_state(P, Rho, T, Tv, gam, delta)
        #instantiate object: 
        cone_ang = 0.197
        M_inf = 2.5
        R = 287

        h0 = GS._h0_given_T_M(R, Tv, T, M_inf)
        T0 = GS._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = GS._a(self, T0)
                else: 
                    self.a0 = GS._a(self, T0)
        gas = gasProps(1.4, 287, T0, p0, delta, Tv)

       
 
        cone_flow_cpg = taylor_maccoll_cpg.TaylorMaccoll_Cone(cone_ang, M_inf, gas)
        cone_flow_tpg = taylor_maccoll_tpg.TaylorMaccoll_Cone(cone_ang, M_inf, gas)
        rel_err = abs((cone_flow_cpg.shock_ang - cone_flow_tpg.shock_ang) / cone_flow_cpg.shock_ang) * 100
        print(rel_err) 
        self.assertAlmostEqual(cone_flow_cpg.shock_ang, cone_flow_tpg.shock_ang,  places=3)



    def test_Taylor_Maccoll_ODE_Solver_against_highTemp(self): 
        #Build the gas state value holder 
        gam = 1.4       
        Rho = 1.225
        T = 3000
        R = 287
        Tv = 3500 
        P = Rho * 287 * T
        delta = 1
        GS_init = gas_state(P, Rho, T, Tv, gam, delta, R=R)
        #instantiate object: 
        cone_ang = math.radians(15)
        M_inf = 4.5 #4.5
        

        h0 = GS._h0_given_T_M(R, Tv, T, M_inf)

        T0 = GS._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = GS._a(self, T0)
                else: 
                    self.a0 = GS._a(self, T0)
        gas = gasProps(1.4, 287, T0, p0, delta, Tv)

        cone_flow_cpg = taylor_maccoll_cpg.TaylorMaccoll_Cone(cone_ang, M_inf, gas)
        cone_flow_tpg = taylor_maccoll_tpg.TaylorMaccoll_Cone(cone_ang, M_inf, gas)
        rel_err = abs((cone_flow_cpg.shock_ang - cone_flow_tpg.shock_ang) / cone_flow_cpg.shock_ang) * 100
        print(rel_err) 
        #self.assertAlmostEqual(cone_flow_cpg.shock_ang, cone_flow_tpg.shock_ang,  places=3)


    # def test_Taylor_Maccoll_against_Data_from_Lampe(self): 
    #     #Page 69
    #     #Build the gas state value holder 
    #     gam = 1.4       
    #     Rho = 1.225
    #     R = 287
    #     Tv = 3500 
    #     delta = 1
    #     M_inf = 8
    #     #instantiate object: 
    #     cone_ang = np.deg2rad([50.88, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 28, 29, 30, 30.5])
    #     Data_set_5_0 = [41.81031, 41.81202, 41.87005, 42.10652, 43.62511, 45.02309, 46.82965, 49.01042, 51.61332, 54.74925, 58.63933, 59.59896, 61.77161, 64.75091, 67.59622]
    #     Data_set_1_5 = [41.81178, 41.86707, 42.09204, 42.62391, 43.55700, 44.90082, 46.65177, 48.76213, 51.25908, 54.24169, 57.86897, 58.74501, 60.56109, 63.01373, 64.53310]
    #     Data_set_0_5 = [41.81161, 41.86481, 42.08150, 42.59778, 43.50216, 44.81982, 46.52185, 48.61556, 51.04943, 53.93270, 57.37846, 58.22053, 59.95754, 62.09982, 63.37443]
    #     #Get the initial conditions given Mach and Theta_v0
    #     Theta_V0 = 1.5
    #     T0 = 0.5 * (Tv / Theta_V0)

    #     #h0 = GS._h(GS, R, Tv, T0, M_inf)

    #     T = GS._T_given_h0_M(GS, R, Tv, T0, M_inf)




    #     P = Rho * 287 * T
    #     GS_init = gas_state(P, Rho, T, Tv, gam, delta)
    #     class gasProps:
    #         def __init__(self, gam, R, T0, T): 
    #             self.gam, self.R, self.T0, self.T = gam, R, T0, T

    #     gas = gasProps(1.4, 287, T0, T)
    #     Shock_ang_calc_tpg = np.zeros(len(cone_ang))
    #     Shock_ang_calc_lampe = np.zeros(len(cone_ang))
    #     perc_error_tpg = np.zeros(len(cone_ang))
    #     perc_error_lampe = np.zeros(len(cone_ang))
    #     for i in range(6, len(cone_ang)):
    #         cone_flow_tpg = taylor_maccoll_tpg.TaylorMaccoll_Cone_TPG(cone_ang[i], M_inf, gas, GS_init)
    #         cone_flow_lampe = taylor_maccoll_Lampe(cone_ang[i], M_inf, gas, GS_init)
    #         Shock_ang_calc_tpg[i] = math.degrees(cone_flow_tpg.shock_ang)
    #         Shock_ang_calc_lampe[i] = math.degrees(cone_flow_lampe.shock_ang)
    #         print(i)
    #         if Theta_V0 == 5: 
    #             perc_error_tpg[i] = abs(Shock_ang_calc_tpg[i] - Data_set_5_0[i]) / Data_set_5_0[i] * 100 
    #             perc_error_lampe[i] = abs(Shock_ang_calc_lampe[i] - Data_set_5_0[i]) / Data_set_5_0[i] * 100 
    #         elif Theta_V0 == 1.5: 
    #             perc_error_lampe[i] = abs(Shock_ang_calc_lampe[i] - Data_set_1_5[i]) / Data_set_1_5[i] * 100 
    #             perc_error_tpg[i] = abs(Shock_ang_calc_tpg[i] - Data_set_1_5[i]) / Data_set_1_5[i] * 100
    #         elif Theta_V0 == 0.5:
    #             perc_error_lampe[i] = abs(Shock_ang_calc_lampe[i] - Data_set_0_5[i]) / Data_set_0_5[i] * 100 
    #             perc_error_tpg[i] = abs(Shock_ang_calc_tpg[i] - Data_set_0_5[i]) / Data_set_0_5[i] * 100
            

    #     print(Theta_V0)
    #     print('Time to check the data against TPG that i made')
    #     print(perc_error_tpg)
    #     print('Time to check the data against replicated Lampe')
    #     print(perc_error_lampe)

   

    def test_Taylor_Maccoll_against_Data_from_CFD(self): 
        #Page 69
        #Build the gas state value holder 
        gam = 1.4       
        Rho = 1.225
        R = 287
        Tv = 3101.6 
        delta = 1
        M_inf = [8.0, 8.0, 8.0, 8.0, 10, 10, 10, 10, 10, 10, 12, 12, 12, 12, 12, 12, 
                 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 18, 18, 18, 20, 20]
        #instantiate object: 
        cone_ang = np.deg2rad([50.88, 45.6, 39.6152, 33.9152, 50.0973,
                               45.0073, 39.1105, 33.4661, 22.5840, 12.0913, 
                               49.6799, 44.6730, 38.8176, 33.1895, 22.3076, 11.7726, 
                               45.6730, 44.4660, 38.6325, 33.0121, 22.1243, 11.5592, 
                               38.5085, 32.8918, 21.9961, 11.4065, 
                               32.8066, 21.9031, 11.2932, 
                               21.8336, 11.2059])
        Data_set_5_0 = [58.2043, 51.5167, 44.4732, 38.0707, 56.5326, 50.2616,
                        43.4111, 37.0907, 25.3029, 14.3641, 
                        55.6491, 49.5592, 42.7980, 36.5127, 24.7220, 13.6750, 
                        50.5531, 49.1257, 42.4120, 36.1431, 24.3389, 13.2182, 
                        42.1539, 35.8930, 24.0718, 12.8939,
                        35.7162, 23.8785, 12.6543, 
                        23.7343, 12.4707]
        #Get the initial conditions given Mach and Theta_v0
        T = 300 #change to 10 degrees if you want to see the CPG difference
        

        #h0 = GS._h(GS, R, Tv, T0, M_inf)

        




        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = GS._a(self, T0)
                else: 
                    self.a0 = GS._a(self, T0)

        

        
        Shock_ang_calc_tpg = np.zeros(len(cone_ang))
        Shock_ang_calc_cpg = np.zeros(len(cone_ang))
        perc_error_tpg = np.zeros(len(cone_ang))
        perc_error_cpg = np.zeros(len(cone_ang))
        for i in range(0, len(cone_ang)):
            h0 = GS._h0_given_T_M(R, Tv, T, M_inf[i])
            T0 = GS._T0_given_h0(R, Tv, h0, T)
            gas = gasProps(1.4, 287, T0, p0, delta, Tv)
            cone_flow_tpg = taylor_maccoll_tpg.TaylorMaccoll_Cone(cone_ang[i], M_inf[i], gas)
            #cone_flow_cpg = taylor_maccoll_cpg.TaylorMaccoll_Cone(cone_ang[i], M_inf[i], gas)
            Shock_ang_calc_tpg[i] = math.degrees(cone_flow_tpg.shock_ang)
            #Shock_ang_calc_cpg[i] = math.degrees(cone_flow_cpg.shock_ang)
            print(i)
            perc_error_tpg[i] = abs(Shock_ang_calc_tpg[i] - Data_set_5_0[i]) / Data_set_5_0[i] * 100 
           # perc_error_cpg[i] = abs(Shock_ang_calc_cpg[i] - Data_set_5_0[i]) / Data_set_5_0[i] * 100 
            
            

        #print(Theta_V0)
        print('Time to check the data against TPG that i made')
        print(perc_error_tpg)
        #print('Time to check the data against CPG that shay did')
        #print(perc_error_cpg)
    
    def test_Taylor_Maccoll_against_temp_from_CFD(self): 
        gam = 1.4       
        Rho = 1.225
        R = 287
        Tv = 3101.6 
        T = 300

        delta = 1
        P = Rho * 287 * T
        M_inf = 8
        cone_ang = np.radians(50.881)

        h0 = GS._h0_given_T_M(R, Tv, T, M_inf)
        T0 = GS._T0_given_h0(R, Tv, h0, T)
        p0 = 101325 
        class gasProps:
            def __init__(self, gam, R, T0, p0, delta, Tv ): 
                self.gam, self.R, self.T0, self.p0f, self.delta, self.Tv = gam, R, T0, p0, delta, Tv
                self.p0_p0f = 1
                if delta == 1: 
                    #Enable Vibrational Mode
                    self.Tv = Tv
                    self.p0f = p0
                    self.a0 = GS._a(self, T0)
                else: 
                    self.a0 = GS._a(self, T0)
        gas = gasProps(1.4, 287, T0, p0, delta, Tv)
        

        cone_flow_tpg = taylor_maccoll_tpg.TaylorMaccoll_Cone(cone_ang, M_inf, gas)
        Shock_ang_calc_tpg = math.degrees(cone_flow_tpg.shock_ang)

        u, v = cone_flow_tpg.f_veloc_uv(math.radians(50.8810))
        T_cs = cone_flow_tpg.shock_values.T2_T1 * T
        rho_cs = cone_flow_tpg.shock_values.rho2_rho1 * Rho
        T_isen  = GS._T_given_T0_V(R, Tv, cone_flow_tpg.T0, u, v)
        rho_isen = GS._isen_rho_given_T(Tv, T_isen, T_cs, rho_cs )
        P_isen = rho_isen * R * T_isen
        print((rho_isen / Rho - 6.5955) / 6.5955* 100 ) 
        print((T_isen / T - 8.80809) / 8.80809  * 100 )
        print((P_isen / P - 58.0946) / 58.0946  * 100 )


if __name__ == "__main__":
    unittest.main()