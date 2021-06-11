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
import duo.zeff.abbema_coeff_calculator as acc
import duo.zeff.common as common

# root finding algorithm references
# https://blogs.sas.com/content/iml/2015/06/22/root-guess.html
# https://www.desmos.com/calculator

#------------------------------------------------------------
#------------------------------------------------------------
def Func(Z, coeff1, coeff2, coeff3, c, g, k):
    f = coeff1 * np.power(Z, c) +\
        coeff2 * np.power(Z, g) +\
        coeff3 * np.power(Z, k)

    return f

#------------------------------------------------------------
#------------------------------------------------------------
def FuncDerivative(Z, coeff1, coeff2, coeff3, c, g, k):
    f = coeff1 * c * np.power(Z, c - 1.0) +\
        coeff2 * g * np.power(Z, g - 1.0) +\
        coeff3 * k * np.power(Z, k - 1.0)

    return f

#------------------------------------------------------------
# Given high and low CT images, calculate ZeffAve image.
# To save time, pre-calculated Abbema coefficients are used
#------------------------------------------------------------
class Abbema:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        self.libPath = None
        self.zeffLib = None

        self.acc = acc.ImprovedAbbemaCoeffCalculator(com)

        self.com = com

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Calculate(self, imageHighKVP, imageLowKVP):
        imageZeff = np.zeros(imageHighKVP.shape)

        # initialGuess = 7.37
        initialGuess = 3.707852

        aHigh = self.acc.a * np.power(self.com.Ehigh, -self.acc.b)
        aLow  = self.acc.a * np.power(self.com.Elow,  -self.acc.b)

        dHigh = self.acc.d * np.power(self.com.Ehigh, -self.acc.f)
        dLow  = self.acc.d * np.power(self.com.Elow,  -self.acc.f)

        hHigh = self.acc.h * np.exp(-self.acc.j * self.com.Ehigh)
        hLow  = self.acc.h * np.exp(-self.acc.j * self.com.Elow)

        for i in range(imageHighKVP.shape[0]):
            for j in range(imageHighKVP.shape[1]):
                muHigh = self.com.ConvertCTNumberToMu(imageHighKVP[i][j], self.com.muWater_Ehigh, self.com.muAir_Ehigh)
                muLow  = self.com.ConvertCTNumberToMu(imageLowKVP[i][j] , self.com.muWater_Elow , self.com.muAir_Elow)

                coeff1 = muHigh * aLow - muLow * aHigh
                coeff2 = muHigh * dLow - muLow * dHigh
                coeff3 = muHigh * hLow - muLow * hHigh

                myArgs = (coeff1, coeff2, coeff3, self.acc.c, self.acc.g, self.acc.k)
                root = scipy.optimize.newton(Func, initialGuess, fprime = FuncDerivative, args = myArgs, disp = False, maxiter=50)
                imageZeff[i][j] = root

        return imageZeff

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateGPU(self, imageHighKVP, imageLowKVP):
        tmg = timer.TimerManager()
        tmg.StartOrResume("total")

        tmg.StartOrResume("newton gpu")
        numElement = imageHighKVP.shape[0] * imageHighKVP.shape[1]
        imageZeff1D = np.zeros(numElement, dtype = np.float64)
        imageHighKVP1D = imageHighKVP.astype(np.float64).flatten()
        imageLowKVP1D  = imageLowKVP.astype(np.float64).flatten()
        self.zeffLib.CalculateZeffWithNewton(imageZeff1D.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        imageHighKVP1D.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        imageLowKVP1D.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        numElement,
        self.com.Ehigh,
        self.com.Elow,
        self.acc.a,
        self.acc.b,
        self.acc.c,
        self.acc.d,
        self.acc.f,
        self.acc.g,
        self.acc.h,
        self.acc.j,
        self.acc.k,
        self.com.muWater_Ehigh,
        self.com.muWater_Elow,
        self.com.muAir_Ehigh,
        self.com.muAir_Elow)
        tmg.Stop("newton gpu")

        imageZeff = np.reshape(imageZeff1D, imageHighKVP.shape)

        return imageZeff

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateSinglePixel(self, pixelHighKVP, pixelLowKVP):
        result = self.Calculate(self, np.array([pixelHighKVP]), np.array([pixelLowKVP]))
        return result[0]

    #------------------------------------------------------------
    #------------------------------------------------------------
    def LoadDll(self, libPath):
        if sys.platform == "linux":
            self.zeffLib = ctypes.CDLL(libPath)
        elif sys.platform == "win32":
            self.zeffLib = ctypes.WinDLL(libPath)

        self.zeffLib.CalculateZeffWithNewton.argtypes = [ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.POINTER(ctypes.c_double),
                                                         ctypes.c_ulonglong,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double,
                                                         ctypes.c_double]
        self.zeffLib.CalculateZeffWithNewton.restype = None

    #------------------------------------------------------------
    # Given material, calculate Z according to Abbema parameterization
    #------------------------------------------------------------
    def CalculateZ(self, mat):
        # assume density = 1.0
        # this is fine, self.CalculateGPU() the CT number will be converted back to mu, which is numerically
        # mu/rho, and in the subsequent Abbema method, density will be canceled out and have no influence
        # on the result

        muHigh = mat.CalculateMacAtE(self.com.Ehigh)
        muLow  = mat.CalculateMacAtE(self.com.Elow)

        huHigh = self.com.ConvertMuToCTNumber(muHigh, self.com.muWater_Ehigh, self.com.muAir_Ehigh)
        huLow  = self.com.ConvertMuToCTNumber(muLow , self.com.muWater_Elow , self.com.muAir_Elow )

        imageHighKVP = np.zeros((1, 1), dtype = np.float64)
        imageLowKVP = np.zeros((1, 1), dtype = np.float64)
        imageZeff = np.zeros((1, 1), dtype = np.float64)

        imageHighKVP[0][0] = huHigh
        imageLowKVP[0][0]  = huLow

        imageZeff = self.CalculateGPU(imageHighKVP, imageLowKVP)

        return imageZeff[0][0]

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateGivenCTNumber(self, imageHighKVP, imageLowKVP):
        return self.CalculateGPU(imageHighKVP, imageLowKVP)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateGivenCTNumberGPU(self, imageHighKVP, imageLowKVP):
        return self.CalculateGPU(imageHighKVP, imageLowKVP)
