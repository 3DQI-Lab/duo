import numpy as np
import duo.core.element as duoelement
import duo.core.material as material
import duo.core.element_table as element_table
import duo.zeff.torikoshi_coeff_calculator as torikoshi_coeff_calculator

#------------------------------------------------------------
#------------------------------------------------------------
class Torikoshi:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, dataPath):
        self.dataPath = dataPath
        self.elementTable = None

        self.Initialize()

        self.tcc = torikoshi_coeff_calculator.TorikoshiCoeffCalculator(self.elementTable)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        # load xs data
        self.elementTable = element_table.ElementTable(self.dataPath)

    #------------------------------------------------------------
    # based on endfb xs data
    #------------------------------------------------------------
    def CalculateZeffAtE(self, mat, energy):
        my_xs_tt = mat.CalculateElectronXSAtE(energy)

        #------------------------------------------------------------
        # plain linear interpolation
        #------------------------------------------------------------
        # # iterate all elements (Z = 1 ~ 100)
        # for Z, element in self.elementTable.elementList.items():
            # y1 = element.CalculateElectronXSAtE(energy)

            # if my_xs_tt < y1:
                # x0 = Z - 1
                # x1 = Z

                # element_prev = mat.elementTable.GetElementByZ(x0)
                # y0 = element_prev.CalculateElectronXSAtE(energy)

                # return (x1 - x0) / (y1 - y0) * (my_xs_tt - y0) + x0

        # # if the desired interval is not found
        # return 0.0

        coeffList = np.array([self.tcc.a5,
                              self.tcc.a4,
                              self.tcc.a3,
                              self.tcc.a2,
                              self.tcc.a1,
                              self.tcc.a0 - my_xs_tt])

        results = np.roots(coeffList)

        root = -1.0
        for result in results:
            if np.isreal(result) and result > 0.0:
                root = np.real(result)

        if root < 0.0:
            raise de.DuoException("--> Torikoshi: root not found.")

        return root

