
#------------------------------------------------------------
#------------------------------------------------------------
class PhotoAtomicXS:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, energy, microXS):
        self.energy = energy # keV
        self.microXS = microXS # barn
        self.mac = 0.0 # mass attenuation coefficient, cm2 / g
