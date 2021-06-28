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
    """Class that manages the parameters of Bourque's formalism,
    using polynomial for curve fitting.

    """


    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, elementTable):
        self.elementTable = elementTable

        self.aList = []

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ParameterizeAtE(self, energy):
        """Calculate parameters at the given photon energy

        :param energy: Photon energy in keV.
        :type energy: float.

        """
        ZList = np.arange(1, 52 + 1)

        xs_list = []
        for Z in ZList:
            element = self.elementTable.elementList[Z]
            xs = element.CalculateElectronXSAtE(energy)
            xs_list.append(xs)

        self.aList = poly.polyfit(ZList, xs_list, 9)


    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateElectronXS(self, Z):
        """Given Z, calculate electron microscopic cross-section.
        This method must be used after a call to :meth:`~ParameterizeAtE`,
        which calculates parameters based on the given energy.

        :param Z: Atomic number.
        :type Z: int.

        """
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
    """Subclass of :class:`BourqueCoeffCalculator`, but using Chebyshev polynomial instead for curve fitting.

    """

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
    """Subclass of :class:`BourqueCoeffCalculator`, but using B-spline instead for curve fitting.

    """

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
    #------------------------------------------------------------
    def FindRoots(self, my_xs_tt):
        """Known issue: for single-element material Z=1 (H, H2, H3 ...), interpolate.sproot()
        is somehow unable to find the root where Z=1
        """
        newT = self.bs.t
        newC = self.bs.c - my_xs_tt
        newK = self.bs.k
        results = interpolate.sproot((newT, newC, newK))



        return results

    #------------------------------------------------------------
    #------------------------------------------------------------
    def EvaluateSpline(self, my_xs_tt, Z):
        """Debugging purpose

        """
        newT = self.bs.t
        newC = self.bs.c - my_xs_tt
        newK = self.bs.k

        result = interpolate.splev(Z, (newT, newC, newK))
        return result







