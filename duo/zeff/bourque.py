import numpy as np
import duo.core.element as duoelement
import duo.core.material as material
import duo.core.element_table as element_table
import duo.zeff.bourque_coeff_calculator as bourque_coeff_calculator
import numpy.polynomial.polynomial as poly
import scipy.interpolate as interpolate
import duo.zeff.common as common
import duo.core.duo_exception as de

#------------------------------------------------------------
#------------------------------------------------------------
def Func6(Z, d0, d1, d2, d3, d4, d5):
    result = d5 * np.power(Z, 5) + d4 * np.power(Z, 4) + d3 * np.power(Z, 3) +\
             d2 * np.power(Z, 2) + d1 * np.power(Z, 1) + d0
    return result

#------------------------------------------------------------
#------------------------------------------------------------
class Bourque:
    """ Critical: By default Bourque parameterizes electron microscopic cross-section (exs) using polynomial approximation.
        At low-energy range, when exs changes abruptly due to absorption, curve fitting would have large error.

        :ivar method: Type of Bourque coefficient calculator to apply.

            - "original": Using polynomial.

            - "chebyshev": Using Chebyshev polynomial.

            - "bspline": Using cubic B-Spline.

        :ivar bcc: Bourque coefficient calculator object.
        :ivar dList: Parameters of DER(Z) polynomial curve fitting.
        :ivar bs: Parameters of Z(DER) B-spline curve fitting.
    """

    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method = "bspline"):
        """Initialization method.

        """
        self.elementTable = com.elementTable
        self.method = method

        self.bcc = None
        if self.method == "original":
            self.bcc = bourque_coeff_calculator.BourqueCoeffCalculator(self.elementTable)
            # print("--> Bourque (original)")
        elif self.method == "chebyshev":
            self.bcc = bourque_coeff_calculator.ChebyshevBourqueCoeffCalculator(self.elementTable)
            # print("--> Bourque (Chebyshev)")
        elif self.method == "bspline":
            self.bcc = bourque_coeff_calculator.BSplineBourqueCoeffCalculator(self.elementTable)
            # print("--> Bourque (B-spline)")

        self.dList = []
        # self.fList = []

        self.bs = None # b-spline object

        self.com = com

        self.derMax = None
        self.derMin = None

        self.Ehigh = 0.0
        self.Elow = 0.0

        # construct water material
        self.water = material.Material("Water", self.com.elementTable)
        self.water.AddElement(1 , weightFraction =  0.111894)
        self.water.AddElement(8 , weightFraction =  0.888106)
        self.water.Commit()
        self.water.density = 1.0

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateZeffAtE(self, mat, energy):
        """Given material and energy, calculate Zeff.
        This method calls :meth:`duo.zeff.bourque_coeff_calculator.BourqueCoeffCalculator.ParameterizeAtE`.
        The curve fitting depends on :attr:`.method` which defaults to "bspline".
        ENDFB cross-section data is used.
        """
        my_xs_tt = mat.CalculateElectronXSAtE(energy)

        # mat.ShowInfo(energy)

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

        #------------------------------------------------------------
        # polynomial parameterization
        #------------------------------------------------------------
        self.bcc.ParameterizeAtE(energy)

        results = self.bcc.FindRoots(my_xs_tt)

        root = -1.0
        allRoots = []
        for result in results:
            if np.isreal(result) and result > 0.0 and result <= 100.0:
                root = np.real(result)
                allRoots.append(root)

        if len(allRoots) > 1:
            # print("--> WARNING: Bourque::CalculateZeffAtE(): more than one roots found!!!", allRoots)
            root = np.amin(allRoots)

        if root < 0.0:
            raise de.DuoException("--> Bourque: root not found.")

        return root

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateDualEnergyRatio(self, Ehigh, Elow, Z):
        """Given Ehigh, Elow, Z, calculate DER using parametric equation of exs in Z

        """
        self.bcc.ParameterizeAtE(Elow)
        exs_low = self.bcc.CalculateElectronXS(Z)

        self.bcc.ParameterizeAtE(Ehigh)
        exs_high = self.bcc.CalculateElectronXS(Z)

        der = exs_low / exs_high

        der *= self.water.CalculateMacAtE(Ehigh) / self.water.CalculateMacAtE(Elow)

        return der

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ParameterizeDualEnergyRatioAndZ(self, Ehigh, Elow):
        """
        Given Ehigh and Elow, derive

        * The parametric equation of DER as a function of Z.

        * The parametric equation of Z as a function of DER.

        Bourque states that Z and gamma are bijective in [1, 38]
        We noticed they are bijective in [1, 36] using the ENDFB library.

        Following Bourque's method, we establish a relation (:attr:`.bs`) between
        Z and DER, and use DER of an unknown material to predict its Z.
        The result would still be the same if we directly use exs_low / exs_high
        without considering water.
        """
        self.Ehigh = Ehigh
        self.Elow = Elow

        derList = []
        ZList = np.arange(1, 36 + 1)

        for Z in ZList:
            der = self.CalculateDualEnergyRatio(self.Ehigh, self.Elow, Z)
            derList.append(der)

        self.dList = poly.polyfit(ZList, derList, 9)

        self.derMax = np.amax(derList)
        self.derMin = np.amin(derList)

        print("    der max = ", self.derMax)
        print("    der min = ", self.derMin)

        # polynomial fitting does not look well enough for Z(DER)
        # self.fList = poly.polyfit(derList, ZList, 9)

        # use spline instead
        # find the B-spline representation of 1-D curve
        (t, c, k) = interpolate.splrep(derList, ZList, k = 3)
        # vector of knots
        # the B-spline coefficients
        # the degree of the spline

        self.bs = interpolate.BSpline(t, c, k, extrapolate = True)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateDualEnergyRatio2(self, Z):
        """Given Z, calculate DER using parametric equation of DER in Z
        """
        return poly.polyval(Z, self.dList)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateZeff2(self, der):
        # return poly.polyval(der, self.fList)
        return self.bs(der)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateGivenCTNumber(self, imageHighKVP, imageLowKVP):
        muHighKVP = self.com.ConvertCTNumberToMu(imageHighKVP, self.com.muWater_Ehigh, self.com.muAir_Ehigh)
        muLowKVP  = self.com.ConvertCTNumberToMu(imageLowKVP , self.com.muWater_Elow , self.com.muAir_Elow )

        return self.Calculate(muHighKVP, muLowKVP)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Calculate(self, muHighKVP, muLowKVP):
        derList = muLowKVP / muHighKVP
        derList *= self.water.CalculateMacAtE(self.Ehigh) / self.water.CalculateMacAtE(self.Elow)

        print("--> before clamp")
        print("    der max = ", np.amax(derList))
        print("    der min = ", np.amin(derList))

        # clamp
        if self.derMax == None or self.derMin == None:
            raise de.DuoException("--> Bourque: derMax or derMin is None.")

        imageMaxValue = self.derMax
        imageMinValue = self.derMin
        derList = np.where(derList > imageMaxValue,  imageMaxValue, derList)
        derList = np.where(derList < imageMinValue,  imageMinValue, derList)

        # self.com.PlotSingleImage(derList, "der.pdf", markWrongData = True)

        print("--> after clamp")
        print("    der max = ", np.amax(derList))
        print("    der min = ", np.amin(derList))

        imageZeff = self.bs(derList)

        # adjust out-of-range data
        upperBoundZ = self.bs(self.derMax)
        lowerBoundZ = self.bs(self.derMin)

        imageZeff = np.where(derList > self.derMax, upperBoundZ, imageZeff)
        imageZeff = np.where(derList < self.derMin, lowerBoundZ, imageZeff)

        return imageZeff

    #------------------------------------------------------------
    #------------------------------------------------------------
    def EvaluateSpline(self, mat, energy, Z):
        """Debugging purpose
        """
        self.bcc.ParameterizeAtE(energy)

        my_xs_tt = mat.CalculateElectronXSAtE(energy)

        result = self.bcc.EvaluateSpline(my_xs_tt, Z)

        msg = "Cubic B-Spline evaluation = {:.16}, xs_tt = {:.16}".format(result, my_xs_tt)
        print(msg)

