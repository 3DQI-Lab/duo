#------------------------------------------------------------
# Generate 4096 x 4096 table for lookup
# Plot the table as 2D image
#------------------------------------------------------------

import os
import sys
import duo.zeff.bourque as bourque
import duo.zeff.taylor as taylor
import duo.zeff.abbema as abbema
import duo.zeff.abbema_coeff_calculator as abbema_coeff_calculator
import duo.zeff.cai as cai
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import duo.zeff.common as common
import duo.core.material as material
import matplotlib.lines as lines
from matplotlib import cm
import matplotlib.colors as mc

#------------------------------------------------------------
#------------------------------------------------------------
class ZeffTableGenerator:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self,
                 dataPath,
                 step = 1,
                 controlPointType = "colonECExtra",
                 use140kVpWithSn = True,
                 use100kVp = False,
                 originalBourque = False):
        self.dataPath = dataPath
        self.fileName = "zeff_table.bin"
        self.imageName = "zeff_table.pdf"
        self.caseSuffix = None
        self.lowHU  = -1000 # start: lowest node
        self.highHU = 3095 # stop (inclusive): highest node
        self.step = step
        self.num = int((self.highHU - self.lowHU + 1) / self.step) # actual number of nodes along one dimension

        self.controlPointType = controlPointType
        self.use140kVpWithSn = use140kVpWithSn
        self.use100kVp = use100kVp

        self.caseSuffix = self.controlPointType
        if self.use100kVp:
            self.caseSuffix += "_100"
        else:
            self.caseSuffix += "_80"

        if self.use140kVpWithSn:
            self.caseSuffix += "_140_Sn"
        else:
            self.caseSuffix += "_140"


        self.zeffOriginal = None
        self.zeffRbf = None

        self.huLowKVPImage = None
        self.huHighKVPImage = None

        self.zeffMin = 1.0
        self.zeffMax = 36.0

        self.com = None

        self.originalBourque = originalBourque

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GenerateZeffTable(self):
        self.InitializeHuLowHighImages()

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # original bourque
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        com_bourque_original = common.Common(self.dataPath,
                                             use140kVpWithSn = self.use140kVpWithSn,
                                             use100kVp = self.use100kVp)

        if self.originalBourque:
            bq = bourque.Bourque(com_bourque_original, method = 'original')
        else:
            bq = bourque.Bourque(com_bourque_original, method = 'bspline')

        bq.ParameterizeDualEnergyRatioAndZ(com_bourque_original.Ehigh,
                                           com_bourque_original.Elow)
        self.imageZeff_bq_original = com_bourque_original.CalculateGivenCTNumber(self.huHighKVPImage,
                                                                    self.huLowKVPImage,
                                                                    bq)

        print("--> imageZeff_bq_original")
        print("    max = {:f} min = {:f}".format(np.amax(self.imageZeff_bq_original), np.amin(self.imageZeff_bq_original)))

        # this class will borrow some useful data from external Common object
        self.com = com_bourque_original

        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ## cai-bourque
        ##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #com_bourque_cai = common.Common(self.dataPath,
        #                                use140kVpWithSn = self.use140kVpWithSn,
        #                                use100kVp = self.use100kVp)
        #cai_bourque = cai.Cai(com_bourque_cai, method = "bourque",
        #                      controlPointType = self.controlPointType)
        #self.imageZeff_bq_cai = com_bourque_cai.CalculateGivenCTNumber(self.huHighKVPImage,
        #                                                  self.huLowKVPImage,
        #                                                  cai_bourque)
        ## clamp zeff within [1, 36]
        #self.imageZeff_bq_cai = self.ClampZeff(self.imageZeff_bq_cai)

        #print("--> imageZeff_bq_cai")
        #print("    max = {:f} min = {:f}".format(np.amax(self.imageZeff_bq_cai), np.amin(self.imageZeff_bq_cai)))

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # von abbema
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        a = abbema.Abbema(self.com)
        libPath = "../../../../duo_build_release/bin/Release/zeff-calculator-abbema-lib.dll"
        a.LoadDll(libPath)
        self.imageZeff_abbema = a.CalculateGivenCTNumberGPU(self.huHighKVPImage, self.huLowKVPImage)

        # clamp zeff within [1, 36]
        self.imageZeff_abbema = self.ClampZeff(self.imageZeff_abbema)

        print("--> imageZeff_abbema")
        print("    max = {:f} min = {:f}".format(np.amax(self.imageZeff_abbema), np.amin(self.imageZeff_abbema)))


        self.diffBqAbbema = self.imageZeff_bq_original - self.imageZeff_abbema

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # cai-taylor
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        com_taylor_cai = common.Common(self.dataPath,
                                       use140kVpWithSn = self.use140kVpWithSn,
                                       use100kVp = self.use100kVp)
        cai_taylor = cai.Cai(com_taylor_cai, method = "taylor",
                             controlPointType = self.controlPointType)
        self.imageZeff_tl_cai = com_taylor_cai.CalculateGivenCTNumber(self.huHighKVPImage,
                                                         self.huLowKVPImage,
                                                         cai_taylor)
        # clamp zeff within [1, 36]
        self.imageZeff_tl_cai = self.ClampZeff(self.imageZeff_tl_cai)

        print("--> imageZeff_tl_cai")
        print("    max = {:f} min = {:f}".format(np.amax(self.imageZeff_tl_cai), np.amin(self.imageZeff_tl_cai)))

        self.PlotZeff(self.imageZeff_bq_original, "zeff_table_bourque_original_" + self.caseSuffix + ".pdf")
        self.PlotZeff(self.imageZeff_abbema, "zeff_table_abbema_" + self.caseSuffix + ".pdf")
        self.PlotZeff(self.imageZeff_tl_cai, "zeff_table_cai_taylor_" + self.caseSuffix + ".pdf", addControlPoints = True, caiObj = cai_taylor)
        self.PlotDiff(self.diffBqAbbema, "zeff_table_diff_" + self.caseSuffix + ".pdf")

        #self.SaveDataToDisk()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ClampZeff(self, zeffMap):
        zeffMap = np.where(zeffMap > self.zeffMax, self.zeffMax, zeffMap)
        zeffMap = np.where(zeffMap < self.zeffMin, self.zeffMin, zeffMap)

        return zeffMap

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetImageCoordinateFromHU(self, HU):
        imageCoord = (HU + 1000.0) / self.step

        return imageCoord

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddControlPointsToImage(self, cpManager, ax):
        # add nist data
        for idx in range(len(cpManager.cpListNist.xList)):
            CTNumber_Ehigh = cpManager.cpListNist.xList[idx]
            CTNumber_Elow  = cpManager.cpListNist.yList[idx]
            ax.plot(self.GetImageCoordinateFromHU(CTNumber_Ehigh),
                    self.GetImageCoordinateFromHU(CTNumber_Elow),
                    marker = 'o',
                    fillstyle = 'none',
                    markersize = 6,
                    color = "#ff0000")

        # add non-nist data
        for idx in range(len(cpManager.cpListCai.xList)):
            CTNumber_Ehigh = cpManager.cpListCai.xList[idx]
            CTNumber_Elow  = cpManager.cpListCai.yList[idx]
            ax.plot(self.GetImageCoordinateFromHU(CTNumber_Ehigh),
                    self.GetImageCoordinateFromHU(CTNumber_Elow),
                    marker = 'x',
                    markersize = 6,
                    color = "#ffff00")

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddDelimitingLinesToImage(self, ax):
        CTNumber_Ehigh_Z1_d1  , CTNumber_Elow_Z1_d1  = self.CalculateCTNumbers(ax, 1 , 1.0)
        CTNumber_Ehigh_Z1_d2  , CTNumber_Elow_Z1_d2  = self.CalculateCTNumbers(ax, 1 , 2.0)
        CTNumber_Ehigh_Z36_d1 , CTNumber_Elow_Z36_d1 = self.CalculateCTNumbers(ax, 30, 0.1)
        CTNumber_Ehigh_Z36_d2 , CTNumber_Elow_Z36_d2 = self.CalculateCTNumbers(ax, 30, 0.2)

        ax.axline((self.GetImageCoordinateFromHU(CTNumber_Ehigh_Z1_d1 ), self.GetImageCoordinateFromHU(CTNumber_Elow_Z1_d1 )),
                  (self.GetImageCoordinateFromHU(CTNumber_Ehigh_Z1_d2 ), self.GetImageCoordinateFromHU(CTNumber_Elow_Z1_d2 )),
                  color = "#000000", linewidth = 1, linestyle = "--")
        ax.axline((self.GetImageCoordinateFromHU(CTNumber_Ehigh_Z36_d1), self.GetImageCoordinateFromHU(CTNumber_Elow_Z36_d1)),
                  (self.GetImageCoordinateFromHU(CTNumber_Ehigh_Z36_d2), self.GetImageCoordinateFromHU(CTNumber_Elow_Z36_d2)),
                  color = "#000000", linewidth = 1, linestyle = "--")

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateCTNumbers(self, ax, Z, density):
        result = self.CalculateCTNumberForOneMaterial(Z, density)

        return result

    #------------------------------------------------------------
    #------------------------------------------------------------
    def InitializeHuLowHighImages(self):
        huLowKVPNodes = np.arange(0, self.num) * self.step + self.lowHU # given start, step and number
        huHighKVPNodes = np.arange(0, self.num) * self.step + self.lowHU
        self.huLowKVPImage, self.huHighKVPImage = np.meshgrid(huLowKVPNodes, huHighKVPNodes, indexing='ij')

        print("    num = ", self.num)
        print("    step = ", self.step)

        fig = plt.figure(figsize = (16, 6))
        plt.subplots_adjust(wspace = 0.1, hspace = 0.1, right = 0.9)

        ax = fig.add_subplot(1, 2, 1)
        m = ax.imshow(self.huLowKVPImage, interpolation = 'none', cmap = "Greys_r")
        ax.set_title("HU low kVp")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        cbar = fig.colorbar(m, ax = ax)
        cbar.ax.locator_params(nbins = 8) # nbins is the max number of bins

        ax = fig.add_subplot(1, 2, 2)
        m = ax.imshow(self.huHighKVPImage, interpolation = 'none', cmap = "Greys_r")
        ax.set_title("HU high kVp")
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        cbar = fig.colorbar(m, ax = ax)
        cbar.ax.locator_params(nbins = 8) # nbins is the max number of bins

        plt.savefig("hu_low_high.pdf", bbox_inches = 'tight')

    #------------------------------------------------------------
    #------------------------------------------------------------
    def LoadAndPlotData(self, fileName, imageName):
        result1D = np.fromfile(fileName, dtype = np.float64, sep='')
        print("--> Shape : ", result1D.shape)

        result = np.reshape(result1D, (self.num, self.num), order='C')

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111)
        m = ax.imshow(result, interpolation = 'none', origin='lower')
        cbar = fig.colorbar(m, ax = ax)

        ax.set_xlabel("CT number at low keV")
        ax.set_xticks([-0.5, self.num - 0.5])
        ax.set_xticklabels(["-1000", "3095"])

        ax.set_ylabel("CT number at high keV")
        ax.set_yticks([-0.5, self.num - 0.5])
        ax.set_yticklabels(["-1000", "3095"])

        plt.savefig(imageName)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def PlotZeff(self, imageZeff, title, addControlPoints = False, caiObj = None):
        # plot
        # intended tick labels: -1000, 0, 3095
        # note that the endpoints -1000 and 3095 are located in the center of the starting and ending pixels
        # not the tips of the axis!!!
        # intended tick locations in terms of pixel index, if step = 1: 0, 1000, 4095
        # actual pixel index along each dimension: 0, 1, 2, ... , self.num - 1
        myTicks = [0, self.GetImageCoordinateFromHU(0.0), self.num - 1]
        myTickLabels = ["-1000", "0", "3095"]

        #fig = plt.figure(figsize = (20, 6))
        fig = plt.figure(figsize = (7, 6))
        plt.subplots_adjust(wspace = 0.2, hspace = 0.1, right = 0.8)

        ax = fig.add_subplot(1, 1, 1)
        mappable = ax.imshow(imageZeff,
                      origin = 'lower',
                      interpolation = 'none',
                      cmap = 'viridis',
                      vmin = self.zeffMin,
                      vmax = self.zeffMax)
        #ax.set_title("Original Bourque")
        ax.set_xticks(myTicks)
        ax.set_yticks(myTicks)
        ax.set_xticklabels(myTickLabels)
        ax.set_yticklabels(myTickLabels)
        ax.axhline(y = myTicks[1], color = "#000000", linestyle = "-", linewidth = 0.5)
        ax.axvline(x = myTicks[1], color = "#000000", linestyle = "-", linewidth = 0.5)
        ax.set_xlabel("CT number at high kVp", labelpad = -10)
        ax.set_ylabel("CT number at low kVp" , labelpad = -20)

        self.AddDelimitingLinesToImage(ax)

        if addControlPoints == True:
            self.AddControlPointsToImage(caiObj.cp, ax)

        cax = plt.axes([0.82, 0.12, 0.04, 0.76]) # left, bottom, width, height
        cbar = fig.colorbar(mappable, cax = cax)
        cbar.ax.locator_params(nbins = 8) # nbins is the max number of bins

        plt.savefig(title, bbox_inches = 'tight')

    #------------------------------------------------------------
    #------------------------------------------------------------
    def PlotDiff(self, imageDiff, title):
        myTicks = [0, self.GetImageCoordinateFromHU(0.0), self.num - 1]
        myTickLabels = ["-1000", "0", "3095"]

        #fig = plt.figure(figsize = (20, 6))
        fig = plt.figure(figsize = (7, 6))
        plt.subplots_adjust(wspace = 0.2, hspace = 0.1, right = 0.8)

        ax = fig.add_subplot(1, 1, 1)
        mappable = ax.imshow(imageDiff,
                      origin = 'lower',
                      interpolation = 'none',
                      cmap = self.CreateNewColormap(imageDiff))

        ax.set_xticks(myTicks)
        ax.set_yticks(myTicks)
        ax.set_xticklabels(myTickLabels)
        ax.set_yticklabels(myTickLabels)
        ax.axhline(y = myTicks[1], color = "#000000", linestyle = "-", linewidth = 0.5)
        ax.axvline(x = myTicks[1], color = "#000000", linestyle = "-", linewidth = 0.5)
        ax.set_xlabel("CT number at high kVp", labelpad = -10)
        ax.set_ylabel("CT number at low kVp" , labelpad = -20)

        self.AddDelimitingLinesToImage(ax)

        cax = plt.axes([0.82, 0.12, 0.04, 0.76]) # left, bottom, width, height
        cbar = fig.colorbar(mappable, cax = cax)
        cbar.ax.locator_params(nbins = 8) # nbins is the max number of bins

        plt.savefig(title, bbox_inches = 'tight')


    #------------------------------------------------------------
    #------------------------------------------------------------
    def SaveDataToDisk(self):
        # imageZeff is a 2-D array
        #         high kVp
        #         |----->
        # low kVp |----->
        #         |----->
        #         .
        #
        # x axis (column direction where j increments) is high kVp
        # y axis (row direction where i increments) is low kVp
        # imageZeff memory layout follows row major where j index changes fastest
        #
        # we wish to have a target data
        #          low kVp
        #          |----->
        # high kVp |----->
        #          |----->
        #          .
        # so just apply transpose

        temp = self.imageZeff_bq_original
        temp = temp.T
        temp = temp.flatten(order = 'C')
        filename = "bourque_original_" + self.caseSuffix + ".bin"
        temp.tofile(filename, sep = "")

        temp = self.imageZeff_bq_cai
        temp = temp.T
        temp = temp.flatten(order = 'C')
        filename = "bourque_RBF_" + self.caseSuffix + ".bin"
        temp.tofile(filename, sep = "")

        temp = self.imageZeff_tl_cai
        temp = temp.T
        temp = temp.flatten(order = 'C')
        filename = "taylor_RBF_" + self.caseSuffix + ".bin"
        temp.tofile(filename, sep = "")

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateCTNumberForOneMaterial(self, Z, density):
        myMaterial = material.Material("material {:d} {:f}".format(Z, density), self.com.elementTable)
        myMaterial.AddElement(Z, atomicFraction = 1)
        myMaterial.Commit()
        myMaterial.density = density

        mu_Elow  = myMaterial.CalculateMacAtE(self.com.Elow) * myMaterial.density
        CTNumber_Elow = self.com.ConvertMuToCTNumber(mu_Elow, self.com.muWater_Elow, self.com.muAir_Elow)

        mu_Ehigh  = myMaterial.CalculateMacAtE(self.com.Ehigh) * myMaterial.density
        CTNumber_Ehigh = self.com.ConvertMuToCTNumber(mu_Ehigh, self.com.muWater_Ehigh, self.com.muAir_Ehigh)

        return (CTNumber_Ehigh, CTNumber_Elow)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CreateNewColormap(self, image):
        imageMaxValue = np.amax(image)
        imageMinValue = np.amin(image)

        # image value (diff value) exceeding upper or lower bound will be specially colored
        upperBoundValue = 0.2
        lowerBoundValue = -0.2
        numBin = 512

        delta = (imageMaxValue - imageMinValue) / numBin

        # new colormap used to mark problematic pixels
        viridis = cm.get_cmap('viridis', numBin)
        newcolors = viridis(np.linspace(0, 1, numBin))
        upperBoundColor = np.array([1.0, 1.0, 0.0, 1])
        lowerBoundColor = np.array([1.0, 0.0, 0.0, 1])

        temp = (upperBoundValue - imageMinValue) / delta
        newcolors[int(temp) : numBin, :] = upperBoundColor

        temp = (lowerBoundValue - imageMinValue) / delta
        newcolors[0 : int(temp) + 1, :] = lowerBoundColor

        newcmp = mc.ListedColormap(newcolors)

        return newcmp