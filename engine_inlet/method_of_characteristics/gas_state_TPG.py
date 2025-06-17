

import math 
import scipy.optimize
import numpy as np 



class GeneralGasBase:
    """
    Base class for implementing a general gas model.

    This class provides the most general implementations of the thermodynamic relations
    needed for a gas model. By overriding the functions specific to the gas model of
    interest a gas model can easily be created. This model will not be efficient.
    """

    pass


_R_UNIV = np.float64(8.31446261815324e3)

_GAS_PROP_DATA = {'air': ('Air',               # string name
                          np.float64(28.966),  # molecular weight [kg/kmol]
                          np.float64(1.4),     # CP ratio of specific heats [-]
                          ),
                  'N2':  ('Nitrogen Gas',       # string name
                          np.float64(28.0134),  # molecular weight [kg/kmol]
                          np.float64(1.4),      # CP ratio of specific heats [-]
                          ),
                  'He':  ('Helium',              # string name
                          np.float64(4.002602),  # molecular weight [kg/kmol]
                          np.float64(5/3),       # CP ratio of specific heats [-]
                          ),
                  'Ne':  ('Neon',              # string name
                          np.float64(20.179),  # molecular weight [kg/kmol]
                          np.float64(5/3),     # CP ratio of specific heats [-]
                          ),
                  'Ar':  ('Argon',             # string name
                          np.float64(39.948),  # molecular weight [kg/kmol]
                          np.float64(5/3),     # CP ratio of specific heats [-]
                          ),
                  'H2':  ('Hydrogen Gas',     # string name
                          np.float64(2.016),  # molecular weight [kg/kmol]
                          np.float64(1.405),  # CP ratio of specific heats [-]
                          ),
                  'O2':  ('Oxygen Gas',         # string name
                          np.float64(31.9988),  # molecular weight [kg/kmol]
                          np.float64(1.4),      # CP ratio of specific heats [-]
                          ),
                  'CO2': ('Carbon Dioxide',   # string name
                          np.float64(44.01),  # molecular weight [kg/kmol]
                          np.float64(1.289),  # CP ratio of specific heats [-]
                          ),
                  }

class CaloricallyPerfectGas: 
      def _Mach(gam, u, v): 
        return np.sqrt((2/(gam-1))/(-1 + 1/(math.sqrt(u**2 + v**2)**2)))#anderson eq 13.81 rearranged\
      
      @staticmethod
      def _a(gasProps,  T):
        
        return np.sqrt(gasProps.gam * gasProps.R * T)
      
      @staticmethod 
      def _a_given_speed(gasProps, u, v): 
           gam = gasProps.gam
           a0 = gasProps.a0
           return np.sqrt(a0**2 - 0.5*(gam-1)*(u**2 + v**2))
      
      @staticmethod
      def _p0_p(gam, M):
        return (1 + M**2*(gam-1)/2)**(gam/(gam-1)) #isentropic total:static pressure ratio
      
      @staticmethod
      def _p2_p1(gam, M): 
        return 1 + (2*gam/(gam+1))*(M**2 - 1) #normal shock static pressure ratio
      
      @staticmethod
      def _T0_T(gam, M): 
        return 1 + (gam-1)*M**2/2 #isentropic total:static temperature ratio
      
      @staticmethod
      def _T2_T1(gam, M): 
        return (1+2*gam*(M**2-1)/(gam+1))*(2 + (gam-1)*M**2)/((gam+1)*M**2)
      
      @staticmethod
      def _rho0_rho(gam, M): 
        return (1 + (gam-1)*M**2/2)**(1/(gam-1)) #isentropic total:static density ratio
      
      
      @staticmethod
      def _get_state_after_TMC(TMC, ang): 
        """"
        M: Mach at the ray
        T: Temperature along the ray
        V: magnitude of velocity along the ray

        """
        u, v = TMC.f_veloc_uv(ang)
        V = np.sqrt(u**2 + v**2)

        u_shock, v_shock = TMC.f_veloc_uv(TMC.shock_ang)
        V_shock = np.sqrt(u_shock**2 + v_shock**2)
        
        gam = TMC.gam
        T0 = TMC.T0 
        R = TMC.R

        a0 = np.sqrt(gam * R * T0)
        a = math.sqrt(a0**2 - 0.5*(gam-1)*V**2)
        M = V / a 
        M_inf = TMC.M_inf

        a_shock = math.sqrt(a0**2 - 0.5*(gam-1)*V_shock**2)
        M_shock = V_shock / a_shock

        
        T_T0 = 1 / CaloricallyPerfectGas._T0_T(gam, M)
        T = T0 * T_T0

        p01_p1 = CaloricallyPerfectGas._p0_p(gam, M_inf)
        p02_p2 = CaloricallyPerfectGas._p0_p(gam, M_shock)
        
        M1_n = M_inf*math.sin(TMC.shock_ang)#get freestream normal mach component 
        p2_p1 = CaloricallyPerfectGas._p2_p1(gam,M1_n) #static pressure ratio across shock
        p02_p01 = p02_p2*p2_p1/p01_p1
        p_p0 = 1 / CaloricallyPerfectGas._p0_p(gam, M) #pressure ratio at point wrt the local stagnation

        rho_rho0 = 1 / CaloricallyPerfectGas._rho0_rho(gam, M)

        return M, T, V, T_T0, p02_p01, p_p0, rho_rho0
      
      @staticmethod
      def _get_state_after_2D_Shock(OBS): 
        """"
        M: Mach at the ray
        T: Temperature along the ray
        V: magnitude of velocity along the ray

        """
        M = OBS.M2
        M_inf = OBS.M1
        V = M * OBS.a2
        T_T0 = OBS.T2_T0
        T = T_T0 * OBS.T0

        p2_p1 = OBS.p2_p1 #static pressure ratio across shock
        p02_p01 = OBS.p02_p01
        p_p0 = 1 / OBS.p02_p2 #pressure ratio at point wrt the local stagnation

        rho_rho0 = 1 / CaloricallyPerfectGas._rho0_rho(OBS.gam, M)

        return M, T, V, T_T0, p02_p01, p_p0, rho_rho0
      
      @staticmethod
      def _get_mesh_point_properties(gasProps, p0_p0f, u, v):
        gam = gasProps.gam
        T0 = gasProps.T0
        a0 = gasProps.a0 

        #speed
        V = np.sqrt(u**2 + v**2)
        a = np.sqrt(a0**2 - 0.5*(gam-1)*V**2)
        M = V/a

        #Temperature
        T = T0/(1+0.5*(gam-1)* M **2) #static temperature
        T_T0 = T / T0

        #Total Pressure
        p_p0 = ((1 + 0.5*(gam-1)* M **2)**(gam/(gam-1)))**-1

        #Density
        rho_rho0 = ((1 + 0.5*(gam-1)* M **2)**(1/(gam-1)))**-1

        #Wrt freestream values 
        p_p0f = p_p0 * p0_p0f
        rho_rho0f = rho_rho0 * p0_p0f #p0 and rho0 are directly proportional via ideal gas law 


     
        return V, M, T, T_T0, p_p0, rho_rho0, p_p0f, rho_rho0f



class ThermallyPerfectGas:
    @staticmethod
    def _energy(R, Tv, T):
            return R * ( 5 / 2 * T + (Tv / (np.exp(Tv/T) - 1)) )
    
    @staticmethod
    def _a(gasProps, T):
        Tv = gasProps.Tv 
        R = gasProps.R

        return np.sqrt(ThermallyPerfectGas._gamma(R, Tv, T) * R * T )
    @staticmethod
    def _a_given_speed(gasProps, u, v): 
         Tv = gasProps.Tv
         R = gasProps.R
         T0 = gasProps.T0

         T = ThermallyPerfectGas._T_given_T0_V(R, Tv, T0, u, v)
         return np.sqrt(ThermallyPerfectGas._gamma(R, Tv, T) * R * T)

    @staticmethod
    def _Cv(R, Tv, T): 
            return R * ( 5 / 2 + (Tv / (2 * T * math.sinh(Tv / (2 * T))))**2 )
        
    @staticmethod
    def _Cp(R, Tv, T):
            return R * (1 + 5 / 2 + (Tv / (2 * T * math.sinh(Tv / (2 * T))))**2 )
    
    @staticmethod
    def _gamma(R, Tv, T):
        #    return   self._Cp(R,Tv,T) / self._Cv(R, Tv, T)
            return  (7 / 2 + (Tv / (2 * T * math.sinh(Tv / (2 * T))))**2 ) / ( 5 / 2 + (Tv / (2 * T * math.sinh(Tv / (2 * T))))**2 )

    def _alpha_v(T): 
            return 1 / T 
    
    def _k_T(P): 
            return 1 / P 
    
    
    @staticmethod
    def _h(R, Tv, T):
          return R * ( 7 / 2 * T + (Tv / (np.exp(Tv/T) - 1)) )

    @staticmethod
    def _h0_given_T_M(R, Tv, T, M ): 
          return ThermallyPerfectGas._h(R, Tv, T) + M**2 / 2 * ( ThermallyPerfectGas._gamma(R, Tv, T) * R * T )
    
    @staticmethod
    def _T0_given_h0(R, Tv, h0, T_guess):

        initial_guess = float(T_guess / 2)
        def f(T):
            return h0 - ThermallyPerfectGas._h(R, Tv, T)

        T0 = scipy.optimize.newton(f, initial_guess)
        return T0
    
    @staticmethod
    def _T_given_T0_M(R, Tv, T0, M): 
        #Get stagnation enthalpy
        h0 = ThermallyPerfectGas._h(R, Tv, T0)

        def f(T):
                return h0 - ThermallyPerfectGas._h(R,Tv,T) - ThermallyPerfectGas._gamma(R, Tv, T) * R * T * M**2 / 2
        
        initial_guess = float(T0 / 2)
        bounds = (0, T0)
        root = scipy.optimize.least_squares(f, initial_guess, bounds=bounds)
        return root.x[0]
    
    @staticmethod
    def _T_given_T0_V(R, Tv, T0, u, v):
        h0 = ThermallyPerfectGas._h(R, Tv, T0)
        h = h0 - (u**2 + v**2) / 2 
        def f(T):
                return h - ThermallyPerfectGas._h(R, Tv, T)    
        initial_guess = float(T0 / 2)
        bounds = (0, T0)
        root = scipy.optimize.least_squares(f, initial_guess, bounds=bounds)
        T = root.x[0]
        return T

    def _isen_rho_given_T(Tv, T2, T1, rho1):
        Theta_1 = Tv / (2 * T1) #Temperature ratio for T1
        Theta_2 = Tv / (2 * T2) # Temperature ratio for T2
        def psi(Theta): 
              return Theta / np.tanh(Theta) - np.log(np.sinh(Theta))
        rho_rat = (3 + 2) / 2 * np.log(T2 / T1) + psi(Theta_2) - psi(Theta_1)
        rho2 =  rho1 * np.e**rho_rat
        return rho2     

    @staticmethod
    def _get_state_after_TMC(TMC, ang): 
          #Function for unpacking Taylor Maccoll Flow 
          #Get values necessary for the IDL 
          u, v = TMC.f_veloc_uv(ang)
          R = TMC.R
          T0 = TMC.T0 
          Tv = TMC.Tv 

          #Need the values cross the shock 
          T_shock = TMC.shock_values.T2_T1 * TMC.shock_values.T1
          rho_shock = TMC.shock_values.rho2_rho1 * TMC.shock_values.rho1
          p02_p01 = TMC.shock_values.p02_p01
          

          T = ThermallyPerfectGas._T_given_T0_V(R, Tv, T0, u, v)
          rho = ThermallyPerfectGas._isen_rho_given_T(Tv, T, T_shock, rho_shock)
          rho0 = ThermallyPerfectGas._isen_rho_given_T(Tv, T0, T_shock, rho_shock)
          P = rho * R * T
          P0 = rho0 * R * T0

          T_T0 = T / T0
          rho_rho0 = rho / rho0
          P_P0 = rho_rho0 * T_T0 #ideal gas law relationship

          V = np.sqrt(u**2 + v**2)
          a = np.sqrt(ThermallyPerfectGas._gamma(R, Tv, T) * R * T)
          M = V / a

          return M, T, V, T_T0, p02_p01, P_P0, rho_rho0
    
    @staticmethod
    def _get_mesh_point_properties(gasProps, p0_p0f, u, v):
         #unpack gasProps
         R = gasProps.R
         T0 = gasProps.T0 
         Tv = gasProps.Tv 
         p0 = p0_p0f * gasProps.p0f
         rho0 = p0 / (R * T0)


         V = np.sqrt(u**2 + v**2)


         T = ThermallyPerfectGas._T_given_T0_V(R, Tv, T0, u, v)
         T_T0 = T / T0
         gam = ThermallyPerfectGas._gamma(R, Tv, T) 

         a = np.sqrt(gam * R * T)
         M = V / a

         rho = ThermallyPerfectGas._isen_rho_given_T(Tv, T, T0, rho0)

         rho_rho0 = rho / rho0 
         p_p0 = rho * R * T / p0 #while adiabatic is assumed, the stagnation ratios are equivalent
         
         p_p0f = rho * R * T / gasProps.p0f
         rho_rho0f = rho_rho0 * p0_p0f
         return V, M, T, T_T0, p_p0, rho_rho0, p_p0f, rho_rho0f
        
    def _get_state_after_2D_Shock(OBS): 
        """"
        M: Mach at the ray
        T: Temperature along the ray
        V: magnitude of velocity along the ray

        """
        M = OBS.M2
        M_inf = OBS.M1
        V = M * OBS.a2
        T_T0 = OBS.T2_T0
        T = T_T0 * OBS.T0

        p2_p1 = OBS.p2_p1 #static pressure ratio across shock
        p02_p01 = OBS.p02_p01
        p_p0 = p2_p1 * OBS.p1 / (OBS.p02_p01 * OBS.p01) #pressure ratio at point wrt the local stagnation

        rho_rho0 = OBS.rho2_rho1 * OBS.rho1 / (OBS.rho02_rho01 * OBS.rho01)

        return M, T, V, T_T0, p02_p01, p_p0, rho_rho0







            
    pass