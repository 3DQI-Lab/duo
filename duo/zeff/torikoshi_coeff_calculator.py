import numpy as np
import duo.core.element as duoelement
import duo.core.material as material
import duo.core.element_table as element_table
import duo.core.material as material
import scipy
import scipy.optimize

#------------------------------------------------------------
#------------------------------------------------------------
def FFunc(T, a, b, c, d, f, g):
    Z, energy = T
    result = a * Z * Z +\
             b * Z * energy +\
             c * energy * energy +\
             d * Z +\
             f * energy +\
             g
    return result

#------------------------------------------------------------
#------------------------------------------------------------
def GFunc(T, a, b, c, d, f, g):
    Z, energy = T
    result = a * Z * Z +\
             b * Z * energy +\
             c * energy * energy +\
             d * Z +\
             f * energy +\
             g
    return result

#------------------------------------------------------------
#------------------------------------------------------------
class TorikoshiCoeffCalculator:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, elementTable):
        self.elementTable = elementTable

        self.a_F = 0.0
        self.b_F = 0.0
        self.c_F = 0.0
        self.d_F = 0.0
        self.f_F = 0.0
        self.g_F = 0.0

        self.a_G = 0.0
        self.b_G = 0.0
        self.c_G = 0.0
        self.d_G = 0.0
        self.f_G = 0.0
        self.g_G = 0.0

        self.Parameterize()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Parameterize(self):
        startZ = 1
        stopZ = 20
        ZList = np.arange(startZ, stopZ + 1)

        startEnergy = 50.0
        stopEnergy = 100.0
        numPoints = 81
        energyList = np.linspace(startEnergy, stopEnergy, numPoints)

        # derive parameters for F(Z, E)
        ZFullList = []
        energyFullList = []
        xsFullList = []

        for Z in ZList:
            element = self.elementTable.elementList[Z]
            for energy in energyList:
                ZFullList.append(Z)
                energyFullList.append(energy)

                xs = element.CalculatePEMicroXSAtE(energy) / np.power(Z, 5.0)
                xsFullList.append(xs)


        result = scipy.optimize.curve_fit(FFunc,
                    (ZFullList, energyFullList),
                    xsFullList)
        self.a_F = result[0][0]
        self.b_F = result[0][1]
        self.c_F = result[0][2]
        self.d_F = result[0][3]
        self.f_F = result[0][4]
        self.g_F = result[0][5]

        print("    a_F = ", self.a_F)
        print("    b_F = ", self.b_F)
        print("    c_F = ", self.c_F)
        print("    d_F = ", self.d_F)
        print("    f_F = ", self.f_F)
        print("    g_F = ", self.g_F)



        # derive parameters for G(Z, E)
        ZFullList = []
        energyFullList = []
        xsFullList = []

        for Z in ZList:
            element = self.elementTable.elementList[Z]
            for energy in energyList:
                ZFullList.append(Z)
                energyFullList.append(energy)

                xs = element.CalculateCSMicroXSAtE(energy) + element.CalculateRLMicroXSAtE(energy)
                xs /= Z
                xsFullList.append(xs)


        result = scipy.optimize.curve_fit(GFunc,
                    (ZFullList, energyFullList),
                    xsFullList)
        self.a_G = result[0][0]
        self.b_G = result[0][1]
        self.c_G = result[0][2]
        self.d_G = result[0][3]
        self.f_G = result[0][4]
        self.g_G = result[0][5]

        print("    a_G = ", self.a_G)
        print("    b_G = ", self.b_G)
        print("    c_G = ", self.c_G)
        print("    d_G = ", self.d_G)
        print("    f_G = ", self.f_G)
        print("    g_G = ", self.g_G)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateElectronXS(self, Z, energy):
        return np.power(Z, 4.0) * FFunc((Z, energy),
                                        self.a_F,
                                        self.b_F,
                                        self.c_F,
                                        self.d_F,
                                        self.f_F,
                                        self.g_F) +\
                                GFunc((Z, energy),
                                      self.a_G,
                                      self.b_G,
                                      self.c_G,
                                      self.d_G,
                                      self.f_G,
                                      self.g_G)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateF(self, Z, energy):
        return FFunc((Z, energy),
                    self.a_F,
                    self.b_F,
                    self.c_F,
                    self.d_F,
                    self.f_F,
                    self.g_F)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateG(self, Z, energy):
        return GFunc((Z, energy),
                    self.a_G,
                    self.b_G,
                    self.c_G,
                    self.d_G,
                    self.f_G,
                    self.g_G)


