"""
Example inlet from Anderson BH fig 18 pg 118 M_inf = 2.5
"""
class inletGeom: 
    def __init__(self): 
        self.cone_ang_deg = 12.5
        self.centerbody_bounds = (0, 3.9)
        self.cowl_bounds = (2,4.1)
        def y_centerbody(x):
            if x >= 0 and x < 2.8: 
                return 0.22169466264*x
            elif x >= 2.8 and x <= 3.8:
                A,B,C,D = 0.12523609835933303, -0.18091239851849594, -0.012534962196626642, -0.10011920449716483
                return A*(x-2.8)**5 + B*(x-2.8)**4 + C*(x-2.8)**3 + D*(x-2.8)**2 + 0.22169466264*x
            else: 
                return None
        self.y_centerbody = y_centerbody
        def y_cowl(x):
            if x >= 2 and x <= 4.1:
                A,B,C,D = 0.014656593603382383, -0.155835602414445, 0.48384724402657875, 0.534568305777872
                return A*x**3 + B*x**2 + C*x + D
            else:
                return None 
        self.y_cowl = y_cowl