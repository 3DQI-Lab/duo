import duo.core.photoatomic_xs as xs
import duo.core.constant as constant
import duo.core.search_interp as si

#------------------------------------------------------------
#------------------------------------------------------------
class Element:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, Z = 0):
        self.Z = Z
        self.symbol = ""
        self.A = 0.0
        self.AWR = 0.0 # atomic weight ratio relative to neutron

        # a dictionary
        # each element is: lib_id (a string) : xsList (a list of PhotoAtomicXS)
        self.xsTable = {}

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetAFromAWR(self):
        # derive A from AWR
        self.A = self.AWR * constant.Endfb.neutronMass

    # #------------------------------------------------------------
    # #------------------------------------------------------------
    # def GetMacFromTotalMicroXS(self):
        # # derive mass attenuation coefficient from microscopic total cross-section
        # for xs in self.xsList:
            # xs.mac = constant.Endfb.Avogadro / self.A * xs.totalMicroXS * 1e-24 # cm2 / g

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Show(self):
        print("-->", self.Z, self.symbol, self.A)
        for key, value in self.xsTable.items():
            result = "    xs name: {0:30s} number of data: {1:d}".format(key, len(value))
            print(result)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateTotalMicroXSAtE(self, energy):
        """It turns out photoatomic total xs should use lin-lin interpolation!!!
        We compared lin-lin with log-log. Using lin-lin, an almost perfect match
        between calculated and reference data was observed.
        ref: https://www-nds.iaea.org/exfor/servlet/efhelp/interp.html

        According to Endfb photoatomic data notes:
        WARNING  - As a result the total cross section
                   MUST NOT be interpolated to define the
                   total between tabulated energies. The
                   ONLY consist way to define the total
                   between tabulated energies is to
                   interpolate all of the partials and
                   add them up.
        So it is incorrect to directly interpolate self.xsTable["total_ref"]
        The total microscopic cross-section must be calculated on the fly!!!

        """

        result = 0.0

        for key, value in self.xsTable.items():
            if key != "total_ref":

                # linear search
                # temp = si.InterpXSLinearSearch(value, energy)

                # binary search
                temp = si.InterpXSBinarySearch(value, energy)

                # print("-->", key, temp)
                result += temp

        return result

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculatePEMicroXSAtE(self, energy):
        result = 0.0

        for key, value in self.xsTable.items():
            if ("[incoherent]" not in key) and\
               ("[coherent]" not in key) and\
               ("pair production" not in key) and\
               ("total_ref" not in key):

                # binary search
                temp = si.InterpXSBinarySearch(value, energy)

                # print("-->", key, temp)
                result += temp

        return result

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateCSMicroXSAtE(self, energy):
        result = 0.0

        for key, value in self.xsTable.items():
            if "[incoherent]" in key:

                # binary search
                temp = si.InterpXSBinarySearch(value, energy)

                # print("-->", key, temp)
                result += temp

        return result

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateRLMicroXSAtE(self, energy):
        result = 0.0

        for key, value in self.xsTable.items():
            if "[coherent]" in key:

                # binary search
                temp = si.InterpXSBinarySearch(value, energy)

                # print("-->", key, temp)
                result += temp

        return result

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateElectronXSAtE(self, energy):
        result = self.CalculateTotalMicroXSAtE(energy)
        result /= self.Z
        return result

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetCSList(self):
        for key, value in self.xsTable.items():
            if "[incoherent]" in key:
                return value

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetRLList(self):
        for key, value in self.xsTable.items():
            if "[coherent]" in key:
                return value

