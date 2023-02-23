import taylor_maccoll as tmc 
import math

cone_ang = math.radians(12.5)
M_inf = 2
        
class gasProps:
    def __init__(self, gam, R, T0): 
        self.gam, self.R, self.T0 = gam, R, T0

gas = gasProps(1.4, 287.05, 288.15)

coneSol = tmc.TaylorMaccoll_Cone(cone_ang, M_inf, gas)
