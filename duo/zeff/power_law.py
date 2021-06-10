import numpy as np
import duo.core.element as duoelement
import duo.core.material as material

#------------------------------------------------------------
#------------------------------------------------------------
class PowerLaw:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        self.m = 1.0

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateZeff(self, mat):
        sum = 0.0
        temp = 0.0

        for Z, ec in mat.elementList.items():
            element = mat.elementTable.GetElementByZ(Z)
            temp += ec.atomicFraction * Z

        for Z, ec in mat.elementList.items():
            element = mat.elementTable.GetElementByZ(Z)
            a = ec.atomicFraction * Z / temp
            sum += a * np.power(Z, self.m)

        return np.power(sum, 1.0 / self.m)

#------------------------------------------------------------
#------------------------------------------------------------
class Mayneord(PowerLaw):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        super().__init__()

        self.m = 2.94
