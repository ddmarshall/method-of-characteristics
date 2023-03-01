import math
import scipy.optimize
import unit_processes as up
"""
This file tests out some methods for incorporating shock points in the characteristic mesh
"""
def get_oblique_shock_angle(M, thet, gam): 
    """
    calculates the shock wave angle required for a given upstream mach number, flow deflection from upstream angle
    Inputs
        M: upstream mach number 
        thet: (rad) flow deflection (positive number) from upstream direction
        gam: ratio of specific heats
    Return: 
        beta: (rad) shock wave angle 
    """
    def thetBetaM(beta, thet, M, gam):
        return (2/math.tan(beta))*(M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta)) + 2) - math.tan(thet)

    alpha = math.asin(1/M)
    beta_weak = scipy.optimize.root_scalar(thetBetaM, args=(thet,M,gam), method='bisect', bracket=[alpha, 0.5*(alpha + math.pi/2)])
    beta_strong = scipy.optimize.root_scalar(thetBetaM, args=(thet,M,gam), method='bisect', bracket=[0.5*(alpha + math.pi/2), math.pi/2])

    return beta_weak.root, beta_strong.root

if __name__ == "__main__":

    thet = math.radians(2.223009757)
    M = 1.306751141
  
    beta_w, beta_s = get_oblique_shock_angle(M, thet, 1.4)
    pass