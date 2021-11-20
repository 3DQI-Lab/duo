import numpy as np
import duo.core.element as duoelement
import duo.core.material as material
import duo.core.element_table as element_table
import duo.zeff.torikoshi_coeff_calculator as torikoshi_coeff_calculator
import numpy.polynomial.polynomial as poly

#------------------------------------------------------------
#------------------------------------------------------------
class Torikoshi:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        self.com = com
        self.elementTable = com.elementTable
        self.tcc = torikoshi_coeff_calculator.TorikoshiCoeffCalculator(self.elementTable)

    #------------------------------------------------------------
    # based on endfb xs data
    #------------------------------------------------------------
    def CalculateZeffAtE(self, mat, energy):
        my_xs_tt = mat.CalculateElectronXSAtE(energy)

        self.tcc.ParameterizeAtE(energy)

        results = self.tcc.FindRoots(my_xs_tt)

        root = -1.0
        allRoots = []
        for result in results:
            if np.isreal(result) and result > 0.0 and result <= 100.0:
                root = np.real(result)
                allRoots.append(root)

        if len(allRoots) > 1:
            root = np.amin(allRoots)

        if root < 0.0:
            raise de.DuoException("--> Torikoshi: root not found.")

        return root

