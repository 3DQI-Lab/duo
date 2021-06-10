import numpy as np
import duo.core.element as duoelement
import duo.core.material as material
import duo.core.element_table as element_table
import scipy.interpolate as interpolate
import duo.zeff.taylor_coeff_calculator as taylor_coeff_calculator
import duo.zeff.common as common
import duo.core.duo_exception as de

#------------------------------------------------------------
#------------------------------------------------------------
class Taylor:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method = "original"):
        self.elementTable = com.elementTable
        self.method = method

        self.tcc = None
        if self.method == "original":
            self.tcc = taylor_coeff_calculator.TaylorCoeffCalculator(self.elementTable)
        elif self.method == "bourque":
            self.tcc = taylor_coeff_calculator.AltTaylorCoeffCalculator(self.elementTable)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateZeffAtE(self, mat, energy):
        my_xs_tt = mat.CalculateTotalMicroXSAtEPerAtom(energy)

        # # linear interpolation
        # # iterate all elements (Z = 1 ~ 100)
        # for Z, element in self.elementTable.elementList.items():
            # y1 = element.CalculateTotalMicroXSAtE(energy)

            # if my_xs_tt < y1:
                # x0 = Z - 1
                # x1 = Z

                # element_prev = mat.elementTable.GetElementByZ(x0)
                # y0 = element_prev.CalculateTotalMicroXSAtE(energy)

                # return (x1 - x0) / (y1 - y0) * (my_xs_tt - y0) + x0

        # # if the desired interval is not found
        # return 0.0

        self.tcc.ParameterizeAtE(energy)

        results = self.tcc.FindRoots(my_xs_tt)

        root = -1.0
        allRoots = []
        for result in results:
            if np.isreal(result) and result > 0.0:
                root = np.real(result)
                allRoots.append(root)

        if len(allRoots) > 1:
            # print("--> WARNING: Bourque::CalculateZeffAtE(): more than one roots found!!!", allRoots)
            root = np.amin(allRoots)

        if root < 0.0:
            raise de.DuoException("--> Taylor: root not found.")

        return root


