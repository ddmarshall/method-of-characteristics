import math
import scipy.optimize
#computes the shock wave angle from upstream mach number of flow deflection angle

def get_oblique_shock_angle(M, thet, gam): 
    """
    calculates the shock wave angle required for a given upstream mach number, flow deflection from upstream angle
    Inputs
        M: upstream mach number 
        thet: (rad) flow deflection (positive number) from upstream direction
        gam: ratio of specific heats
    Return: 
        thet: (rad) shock wave angle 
    """
    def thetBetaM(beta, thet, M, gam):
        return (2/math.tan(beta))*(M**2*math.sin(beta)**2 - 1)/(M**2*(gam + math.cos(2*beta)) + 2) - math.tan(thet)

    alpha = math.asin(1/M)
    thet = scipy.optimize.fsolve(thetBetaM, args=(thet,M,gam), x0=0.5*(alpha + math.pi/2))
    return thet

if __name__ == "__main__":
    M = 2
    thet = math.radians(16)
    gam = 1.4
    thet = get_oblique_shock_angle(M, thet, gam)
    print(f"shock wave angle: {math.degrees(thet)} deg")