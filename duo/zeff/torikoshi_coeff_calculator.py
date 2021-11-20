import numpy as np
import duo.core.element as duoelement
import duo.core.material as material
import duo.core.element_table as element_table
import duo.core.material as material
import numpy.polynomial.polynomial as poly

#------------------------------------------------------------
#------------------------------------------------------------
class TorikoshiCoeffCalculator:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, elementTable):
        self.elementTable = elementTable

        self.fList = []
        self.gList = []
        self.aList = []

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ParameterizeAtE(self, energy):
        """Calculate parameters at the given photon energy

        :param energy: Photon energy in keV.
        :type energy: float.

        """

        ZList = np.arange(1, 20 + 1)

        # the new numpy.polynomial.polynomial module
        # uses the following order:
        # p(x) = p[0] +
        #        p[1] * x +
        #        p[deg - 1] * x^(deg - 1) +
        #        p[deg] * x^deg
        # the polynomial has a total of (deg + 1) terms
        # p[0] is the constant term
        # p[deg] is the coef of the highest degree term

        # derive parameters for F(Z, E)
        fDegree = 4 #
        xs_list = []
        for Z in ZList:
            element = self.elementTable.elementList[Z]
            xs = element.CalculatePEMicroXSAtE(energy) / np.power(Z, 5.0)
            xs_list.append(xs)
        self.fList = poly.polyfit(ZList, xs_list, fDegree)

        # derive parameters for G(Z, E)
        gDegree = 3
        xs_list = []
        for Z in ZList:
            element = self.elementTable.elementList[Z]
            xs = element.CalculateCSMicroXSAtE(energy) + element.CalculateRLMicroXSAtE(energy)
            xs /= Z
            xs_list.append(xs)
        self.gList = poly.polyfit(ZList, xs_list, gDegree)

        highestDeg = fDegree + 4
        self.aList = [0 for idx in range(highestDeg + 1)]

        # combine coefficients
        for idx in range(gDegree + 1):
            self.aList[idx] += self.gList[idx]

        for idx in range(fDegree + 1):
            self.aList[idx + 4] += self.fList[idx]


    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateElectronXS(self, Z):
        """Given Z, calculate electron microscopic cross-section.
        This method must be used after a call to :meth:`~ParameterizeAtE`,
        which calculates parameters based on the given energy.

        :param Z: Atomic number.
        :type Z: int.

        """
        return np.power(Z, 4.0) * poly.polyval(Z, self.fList) + poly.polyval(Z, self.gList)


    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateF(self, Z):
        """Given Z, calculate F(energy, Z) function defined in Torikoshi 2003.
        This method must be used after a call to :meth:`~ParameterizeAtE`,
        which calculates parameters based on the given energy.

        :param Z: Atomic number.
        :type Z: int.

        """

        return poly.polyval(Z, self.fList)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateG(self, Z):
        """Given Z, calculate G(energy, Z) function defined in Torikoshi 2003.
        This method must be used after a call to :meth:`~ParameterizeAtE`,
        which calculates parameters based on the given energy.

        :param Z: Atomic number.
        :type Z: int.

        """

        return poly.polyval(Z, self.gList)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def FindRoots(self, my_xs_tt):
        coeffList = np.copy(self.aList)
        coeffList[0] -= my_xs_tt
        results = poly.polyroots(coeffList)

        return results
