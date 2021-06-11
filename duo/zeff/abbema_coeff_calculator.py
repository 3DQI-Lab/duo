import os
import sys
import numpy as np
import duo.zeff.nist as nist
import duo.core.element_table as element_table
import duo.core.material as material
import duo.core.mixture as mixture
import duo.core.dicom_decoder as dicom_decoder
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import duo.core.timer as timer
import scipy
import scipy.optimize
import pydicom
from matplotlib import cm
import duo.zeff.enhance as enhance
import ctypes
import duo.zeff.common as common





#------------------------------------------------------------
#------------------------------------------------------------
class AbbemaCoeffCalculator:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        self.b = 0.0
        self.f = 0.0
        self.j = 0.0
        self.c = 0.0
        self.g = 0.0
        self.k = 0.0
        self.a = 0.0
        self.d = 0.0
        self.h = 0.0

        self.elementTable = com.elementTable

        self.Initialize()

        self.Show()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        pass

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Show(self):
        print("--> Parameter")
        print("    a = {0:.20f}".format(self.a))
        print("    b = {0:.20f}".format(self.b))
        print("    c = {0:.20f}".format(self.c))

        print("    d = {0:.20f}".format(self.d))
        print("    f = {0:.20f}".format(self.f))
        print("    g = {0:.20f}".format(self.g))

        print("    h = {0:.20f}".format(self.h))
        print("    j = {0:.20f}".format(self.j))
        print("    k = {0:.20f}".format(self.k))

    #------------------------------------------------------------
    # Given E and Z, calculate PE micro xs according to Abbema parameterization
    #------------------------------------------------------------
    def CalculatePEMicroXS(self, E, Z):
        return self.a * np.power(E, -self.b) * np.power(Z, self.c)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateCSMicroXS(self, E, Z):
        return self.h * np.exp(-self.j * E) * np.power(Z, self.k)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateRLMicroXS(self, E, Z):
        return self.d * np.power(E, -self.f) * np.power(Z, self.g)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateTotalMicroXS(self, E, Z):
        result = self.CalculatePEMicroXS(E, Z) + self.CalculateCSMicroXS(E, Z) + self.CalculateRLMicroXS(E, Z)
        return result









#------------------------------------------------------------
# Reproduce Abbema parameters using curve fitting
#
# The parameters b,
# f and j were calculated using cross section data of oxygen
# in the energy range of 50-100 keV.
#
# The parameters c, g
# and k were calculated at a mean effective energy of
# 60.61 keV (mean of 51.93 and 69.28 keV) and for atomic
# numbers ranging from 6 to 20.
#
# Finally, the parameters a,
# d and h were determined at the cross section of oxygen
# at 60.61 keV.
#------------------------------------------------------------
class ReproducedAbbemaCoeffCalculator(AbbemaCoeffCalculator):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        super().__init__(com)

        # reference result
        # a = 30.89327040265194312951
        # b = 3.30064464966206916330
        # c = 4.46154843728991234997
        # d = 4.15364737197781863642
        # f = 1.85157275088228523430
        # g = 2.53799461281200366969
        # h = 0.67492767468508230166
        # j = 0.00198782715752225062
        # k = 0.93731306997511498746

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        super().Initialize()

        self.Derive1()
        self.Derive2()
        self.Derive3()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Derive1(self):
        Z = 8
        element = self.elementTable.elementList[Z]

        # PE: get b
        xs_pe_list = []
        startEnergy = 50.0
        stopEnergy = 100.0
        numPoints = 100
        energy_pe_list = np.linspace(startEnergy, stopEnergy, numPoints)
        for energy in energy_pe_list:
            xsPE = element.CalculatePEMicroXSAtE(energy)
            xs_pe_list.append(xsPE)

        result = scipy.optimize.curve_fit(lambda x, b, t : t * np.power(x, -b),
                    energy_pe_list,
                    xs_pe_list)
        self.b = result[0][0]
        t = result[0][1]

        fit_xs_pe_list = []
        for energy in energy_pe_list:
            y = t * np.power(energy, -self.b)
            fit_xs_pe_list.append(y)


        # CS: get j
        energy_cs_list = []
        xs_cs_list = []
        for xs in element.GetCSList():
            if xs.energy >= startEnergy and xs.energy <= stopEnergy:
                energy_cs_list.append(xs.energy)
                xs_cs_list.append(xs.microXS)

        # initial guess from Cai's draft paper
        t_init = 0.688 * np.power(Z, 0.928)
        j_init = 0.00198

        result = scipy.optimize.curve_fit(lambda x, j, t : t * np.exp(-j * x),
                    energy_cs_list,
                    xs_cs_list,
                    p0 = (j_init, t_init)) # for exp, initial guess must be provided to prevent wrong results
        self.j = result[0][0]
        t = result[0][1]

        fit_xs_cs_list = []
        for energy in energy_cs_list:
            y = t * np.exp(- self.j * energy)
            fit_xs_cs_list.append(y)

        # RL: get f
        energy_rl_list = []
        xs_rl_list = []
        for xs in element.GetRLList():
            if xs.energy >= startEnergy and xs.energy <= stopEnergy:
                energy_rl_list.append(xs.energy)
                xs_rl_list.append(xs.microXS)

        result = scipy.optimize.curve_fit(lambda x, f, t : t * np.power(x, -f),
                    energy_rl_list,
                    xs_rl_list)
        self.f = result[0][0]
        t = result[0][1]

        fit_xs_rl_list = []
        for energy in energy_rl_list:
            y = t * np.power(energy, -self.f)
            fit_xs_rl_list.append(y)

        # plot
        fig = plt.figure(figsize = (10, 16))
        ax = fig.add_subplot(111)
        ax.semilogy(energy_pe_list, xs_pe_list, label = 'pe', linestyle = 'None', marker='x')
        ax.semilogy(energy_pe_list, fit_xs_pe_list, label = 'pe fit', linestyle = '-', marker='None')

        ax.semilogy(energy_cs_list, xs_cs_list, label = 'cs', linestyle = 'None', marker='^')
        ax.semilogy(energy_cs_list, fit_xs_cs_list, label = 'cs fit', linestyle = '-', marker='None')

        ax.semilogy(energy_rl_list, xs_rl_list, label = 'rl', linestyle = 'None', marker='o')
        ax.semilogy(energy_rl_list, fit_xs_rl_list, label = 'rl fit', linestyle = '-', marker='None')

        legend = ax.legend(loc = 'best', shadow = True, edgecolor = '#000000')
        ax.set_xlabel("energy [keV]")
        ax.set_ylabel("microscopic cross-section [barn]")
        plt.savefig("bfj_plot.pdf")

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Derive2(self):
        startZ = 6
        stopZ = 20

        energy = (51.93 + 69.28) / 2.0

        xs_pe_list = []
        xs_cs_list = []
        xs_rl_list = []

        ZList = np.arange(startZ, stopZ + 1)

        for Z in ZList:
            element = self.elementTable.elementList[Z]

            xsPE = element.CalculatePEMicroXSAtE(energy)
            xs_pe_list.append(xsPE)

            xsCS = element.CalculateCSMicroXSAtE(energy)
            xs_cs_list.append(xsCS)

            xsRL = element.CalculateRLMicroXSAtE(energy)
            xs_rl_list.append(xsRL)

        result = scipy.optimize.curve_fit(lambda x, c, t : t * np.power(x, c),
                    ZList,
                    xs_pe_list)
        self.c = result[0][0]
        t1 = result[0][1]

        result = scipy.optimize.curve_fit(lambda x, k, t : t * np.power(x, k),
                    ZList,
                    xs_cs_list)
        self.k = result[0][0]
        t2 = result[0][1]

        result = scipy.optimize.curve_fit(lambda x, g, t : t * np.power(x, g),
                    ZList,
                    xs_rl_list)
        self.g = result[0][0]
        t3 = result[0][1]


        fit_xs_pe_list = []
        fit_xs_cs_list = []
        fit_xs_rl_list = []
        for Z in ZList:
            y = t1 * np.power(Z, self.c)
            fit_xs_pe_list.append(y)

            y = t2 * np.power(Z, self.k)
            fit_xs_cs_list.append(y)

            y = t3 * np.power(Z, self.g)
            fit_xs_rl_list.append(y)

        # plot
        fig = plt.figure(figsize = (10, 16))
        ax = fig.add_subplot(111)
        ax.semilogy(ZList, xs_pe_list, label = 'pe', linestyle = 'None', marker='x')
        ax.semilogy(ZList, fit_xs_pe_list, label = 'pe fit', linestyle = '-', marker='None')

        ax.semilogy(ZList, xs_cs_list, label = 'cs', linestyle = 'None', marker='^')
        ax.semilogy(ZList, fit_xs_cs_list, label = 'cs fit', linestyle = '-', marker='None')

        ax.semilogy(ZList, xs_rl_list, label = 'rl', linestyle = 'None', marker='o')
        ax.semilogy(ZList, fit_xs_rl_list, label = 'rl fit', linestyle = '-', marker='None')

        legend = ax.legend(loc = 'best', shadow = True, edgecolor = '#000000')
        ax.set_xlabel("Z")
        ax.set_ylabel("microscopic cross-section [barn]")
        plt.savefig("cgk_plot.pdf")

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Derive3(self):
        Z = 8
        element = self.elementTable.elementList[Z]

        energy = (51.93 + 69.28) / 2.0

        self.a = element.CalculatePEMicroXSAtE(energy) / np.power(energy, -self.b) / np.power(Z, self.c)
        self.h = element.CalculateCSMicroXSAtE(energy) / np.exp(- energy * self.j) / np.power(Z, self.k)
        self.d = element.CalculateRLMicroXSAtE(energy) / np.power(energy, -self.f) / np.power(Z, self.g)

        print("    pe = {0:.16f}".format(element.CalculatePEMicroXSAtE(energy)))
        print("    pe fit = {0:.16f}".format(self.CalculatePEMicroXS(energy, Z)))

        print("    cs = {0:.16f}".format(element.CalculateCSMicroXSAtE(energy)))
        print("    cs fit = {0:.16f}".format(self.CalculateCSMicroXS(energy, Z)))

        print("    rl = {0:.16f}".format(element.CalculateRLMicroXSAtE(energy)))
        print("    rl fit = {0:.16f}".format(self.CalculateRLMicroXS(energy, Z)))







#------------------------------------------------------------
#------------------------------------------------------------
def FuncPE(T, a, b, c):
    Z, energy = T
    return a * np.power(energy, -b) * np.power(Z, c)

#------------------------------------------------------------
#------------------------------------------------------------
def FuncCS(T, h, j, k):
    Z, energy = T
    return h * np.exp(- energy * j) * np.power(Z, k)

#------------------------------------------------------------
#------------------------------------------------------------
def FuncRL(T, d, f, g):
    Z, energy = T
    return d * np.power(energy, -f) * np.power(Z, g)

#------------------------------------------------------------
# Improve Abbema coefficient calculation where curve fitting applies to
# a range of E and Z at the same time
#
# Critical note:
# From the xs plot where x is energy (linear), and y is xs (log),
# the new method seemingly causes more discrepancy for smaller xs values.
# But in actuality, this method does improve upon the original Abbema:
# the curve fitting occurs to linear space, not log space!!!
#------------------------------------------------------------
class ImprovedAbbemaCoeffCalculator(AbbemaCoeffCalculator):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        super().__init__(com)

        # reference result
        # when consider Z = 1 ~ 20
        # a = 18.18567267773311613155
        # b = 3.14004359694221291122
        # c = 4.45272467756578027576
        # d = 3.43276854093298400272
        # f = 1.79278994911922029409
        # g = 2.54404644915319000376
        # h = 0.62376590644053797607
        # j = 0.00147582777946576068
        # k = 0.95521150397177878588

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        super().Initialize()

        self.Derive1()
        self.Derive2()
        self.Derive3()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Derive1(self):
        startZ = 1
        stopZ = 20
        ZList = np.arange(startZ, stopZ + 1)

        startEnergy = 50.0
        stopEnergy = 100.0
        numPoints = 100
        energyList = np.linspace(startEnergy, stopEnergy, numPoints)

        ZFullList = []
        energyFullList = []
        xsFullList = []

        for Z in ZList:
            element = self.elementTable.elementList[Z]
            for energy in energyList:
                ZFullList.append(Z)
                energyFullList.append(energy)

                xs = element.CalculatePEMicroXSAtE(energy)
                xsFullList.append(xs)


        result = scipy.optimize.curve_fit(FuncPE,
                    (ZFullList, energyFullList),
                    xsFullList)
        self.a = result[0][0]
        self.b = result[0][1]
        self.c = result[0][2]

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Derive2(self):
        startZ = 1
        stopZ = 20
        ZList = np.arange(startZ, stopZ + 1)

        startEnergy = 50.0
        stopEnergy = 100.0
        numPoints = 100
        energyList = np.linspace(startEnergy, stopEnergy, numPoints)

        ZFullList = []
        energyFullList = []
        xsFullList = []

        for Z in ZList:
            element = self.elementTable.elementList[Z]
            for energy in energyList:
                ZFullList.append(Z)
                energyFullList.append(energy)

                xs = element.CalculateCSMicroXSAtE(energy)
                xsFullList.append(xs)

        # initial guess from Cai's draft paper
        h_init = 0.688
        j_init = 0.00198
        k_init = 0.928
        result = scipy.optimize.curve_fit(FuncCS,
                    (ZFullList, energyFullList),
                    xsFullList,
                    p0 = (h_init, j_init, k_init)) # for exp, initial guess must be provided to prevent wrong results
        self.h = result[0][0]
        self.j = result[0][1]
        self.k = result[0][2]

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Derive3(self):
        startZ = 1
        stopZ = 20
        ZList = np.arange(startZ, stopZ + 1)

        startEnergy = 50.0
        stopEnergy = 100.0
        numPoints = 100
        energyList = np.linspace(startEnergy, stopEnergy, numPoints)

        ZFullList = []
        energyFullList = []
        xsFullList = []

        for Z in ZList:
            element = self.elementTable.elementList[Z]
            for energy in energyList:
                ZFullList.append(Z)
                energyFullList.append(energy)

                xs = element.CalculateRLMicroXSAtE(energy)
                xsFullList.append(xs)


        result = scipy.optimize.curve_fit(FuncRL,
                    (ZFullList, energyFullList),
                    xsFullList)
        self.d = result[0][0]
        self.f = result[0][1]
        self.g = result[0][2]






#------------------------------------------------------------
# Parameters are from Abbema paper without any change
#------------------------------------------------------------
class OriginalAbbemaCoeffCalculator(AbbemaCoeffCalculator):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        super().__init__(com)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        super().Initialize()

        self.a = 22.3
        self.b = 3.302
        self.c = 4.62

        self.d = 3.71
        self.f = 1.856
        self.g = 2.60

        self.h = 0.672
        self.j = 0.00197
        self.k = 0.939






#------------------------------------------------------------
# Parameters are calculated by a researcher in Dr. Cai's group
#------------------------------------------------------------
class CaiModifiedAbbemaCoeffCalculator(AbbemaCoeffCalculator):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        super().__init__(com)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        super().Initialize()

        self.a = 20.7
        self.b = 3.303
        self.c = 4.66

        self.d = 3.54
        self.f = 1.857
        self.g = 2.63

        self.h = 0.688
        self.j = 0.00198
        self.k = 0.928
