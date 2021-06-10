import numpy as np
import duo.core.element as duoelement
import duo.core.material as material
import duo.core.element_table as element_table
import duo.core.material as material
import numpy.polynomial.polynomial as poly
import numpy.polynomial.chebyshev as chebyshev
import scipy.interpolate as interpolate

#------------------------------------------------------------
#------------------------------------------------------------
class BourqueCoeffCalculator:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, elementTable):
        self.elementTable = elementTable

        self.aList = []

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ParameterizeAtE(self, energy):
        ZList = np.arange(1, 52 + 1)

        xs_list = []
        for Z in ZList:
            element = self.elementTable.elementList[Z]
            xs = element.CalculateElectronXSAtE(energy)
            xs_list.append(xs)

        self.aList = poly.polyfit(ZList, xs_list, 9)


    #------------------------------------------------------------
    # base on parameterized xs data, for the same energy
    # that has been passed to ParameterizeAtE(self, energy)
    #------------------------------------------------------------
    def CalculateElectronXS(self, Z):
        return poly.polyval(Z, self.aList)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def FindRoots(self, my_xs_tt):
        coeffList = np.copy(self.aList)
        coeffList[0] -= my_xs_tt
        results = poly.polyroots(coeffList)

        return results

#------------------------------------------------------------
#------------------------------------------------------------
class ChebyshevBourqueCoeffCalculator(BourqueCoeffCalculator):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, elementTable):
        super().__init__(elementTable)


    #------------------------------------------------------------
    #------------------------------------------------------------
    def ParameterizeAtE(self, energy):
        ZList = np.arange(1, 52 + 1)

        xs_list = []
        for Z in ZList:
            element = self.elementTable.elementList[Z]
            xs = element.CalculateElectronXSAtE(energy)
            xs_list.append(xs)

        self.aList = chebyshev.chebfit(ZList, xs_list, 9)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateElectronXS(self, Z):
        return chebyshev.chebval(Z, self.aList)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def FindRoots(self, my_xs_tt):
        coeffList = np.copy(self.aList)
        coeffList[0] -= my_xs_tt
        results = chebyshev.chebroots(coeffList)

        return results

#------------------------------------------------------------
#------------------------------------------------------------
class BSplineBourqueCoeffCalculator(BourqueCoeffCalculator):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, elementTable):
        super().__init__(elementTable)
        self.bs = None # b-spline object

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ParameterizeAtE(self, energy):
        ZList = np.arange(1, 52 + 1)

        xs_list = []
        for Z in ZList:
            element = self.elementTable.elementList[Z]
            xs = element.CalculateElectronXSAtE(energy)
            xs_list.append(xs)

        # find the B-spline representation of 1-D curve
        (t, c, k) = interpolate.splrep(ZList, xs_list, k = 3)
        # vector of knots
        # the B-spline coefficients
        # the degree of the spline

        self.bs = interpolate.BSpline(t, c, k, extrapolate = True)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateElectronXS(self, Z):
        return self.bs(Z)

    #------------------------------------------------------------
    # Known issue: for purity Z=1 (H, H2, H3 ...), interpolate.sproot()
    # is somehow unable to find the root where Z=1
    #------------------------------------------------------------
    def FindRoots(self, my_xs_tt):
        newT = self.bs.t
        newC = self.bs.c - my_xs_tt
        newK = self.bs.k
        results = interpolate.sproot((newT, newC, newK))



        return results

    #------------------------------------------------------------
    # debugging purpose
    #------------------------------------------------------------
    def EvaluateSpline(self, my_xs_tt, Z):
        newT = self.bs.t
        newC = self.bs.c - my_xs_tt
        newK = self.bs.k

        result = interpolate.splev(Z, (newT, newC, newK))
        return result







