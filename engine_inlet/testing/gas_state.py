
class gas_state: 

    def __init__(self, P, Rho, T, Tv, gam, delta, R=None):

        self.P = P
        self.rho = Rho
        self.T = T
        self.Tv = Tv
        self.gam = gam
        self.delta = delta
        if R is not None: self.R = R
        else: self.R = None
        pass

    pass