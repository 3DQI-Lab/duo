import numpy as np
import duo.core.element as duoelement
import duo.core.material as material
import duo.core.element_table as element_table
import duo.core.material as material
import scipy
import scipy.optimize
import scipy.interpolate as interpolate

#------------------------------------------------------------
#------------------------------------------------------------
class TaylorCoeffCalculator:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, elementTable):
        self.elementTable = elementTable

        self.bs = None # b-spline object

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ParameterizeAtE(self, energy):
        # iterate all elements (Z = 1 ~ 100)
        ZList = []
        xs_tt_list = []
        for Z, element in self.elementTable.elementList.items():
            ZList.append(Z)

            xs_tt = element.CalculateTotalMicroXSAtE(energy)
            xs_tt_list.append(xs_tt)

        # find the B-spline representation of 1-D curve
        (t, c, k) = interpolate.splrep(ZList, xs_tt_list, k = 3)
        # vector of knots
        # the B-spline coefficients
        # the degree of the spline

        self.bs = interpolate.BSpline(t, c, k, extrapolate = True)

    #------------------------------------------------------------
    # base on parameterized xs data, for the same energy
    # that has been passed to ParameterizeAtE(self, energy)
    #------------------------------------------------------------
    def CalculateXS(self, Z):
        return self.bs(Z)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def FindRoots(self, my_xs_tt):
        newT = self.bs.t
        newC = self.bs.c - my_xs_tt
        newK = self.bs.k
        results = interpolate.sproot((newT, newC, newK))

        return results




#------------------------------------------------------------
#------------------------------------------------------------
def Func6(Z, a0, a1, a2, a3, a4, a5):
    result = a5 * np.power(Z, 5) + a4 * np.power(Z, 4) + a3 * np.power(Z, 3) +\
             a2 * np.power(Z, 2) + a1 * np.power(Z, 1) + a0

    # result = a5 * np.power(Z, 6) + a4 * np.power(Z, 5) + a3 * np.power(Z, 4) +\
             # a2 * np.power(Z, 3) + a1 * np.power(Z, 2) + a0 * Z
    return result

#------------------------------------------------------------
#------------------------------------------------------------
class AltTaylorCoeffCalculator(TaylorCoeffCalculator):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, elementTable):
        self.elementTable = elementTable

        self.a0 = 0.0
        self.a1 = 0.0
        self.a2 = 0.0
        self.a3 = 0.0
        self.a4 = 0.0
        self.a5 = 0.0

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ParameterizeAtE(self, energy):
        ZList = np.arange(1, 52 + 1)

        xs_list = []
        for Z in ZList:
            element = self.elementTable.elementList[Z]
            xs = element.CalculateTotalMicroXSAtE(energy)
            xs_list.append(xs)

        result = scipy.optimize.curve_fit(Func6,
                 ZList,
                 xs_list)
        self.a0 = result[0][0]
        self.a1 = result[0][1]
        self.a2 = result[0][2]
        self.a3 = result[0][3]
        self.a4 = result[0][4]
        self.a5 = result[0][5]

    #------------------------------------------------------------
    # base on parameterized xs data, for the same energy
    # that has been passed to ParameterizeAtE(self, energy)
    #------------------------------------------------------------
    def CalculateXS(self, Z):
        return Func6(Z,
                    self.a0,
                    self.a1,
                    self.a2,
                    self.a3,
                    self.a4,
                    self.a5)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def FindRoots(self, my_xs_tt):
        coeffList = np.array([self.a5,
                              self.a4,
                              self.a3,
                              self.a2,
                              self.a1,
                              self.a0 - my_xs_tt])

        results = np.roots(coeffList)

        return results



