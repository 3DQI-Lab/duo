import duo.zeff.common as common
import duo.application.application as appl
import os
import sys
import scipy
import pydicom
import duo.core.material as material
import duo.core.element_table as element_table
import numpy as np
import duo.core.timer as timer
import duo.core.dicom_decoder as dicom_decoder
import duo.zeff.enhance as enhance
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mc
import duo.zeff.blend as blend
import pydicom.pixel_data_handlers.util as pdutil
import cv2

#------------------------------------------------------------
#------------------------------------------------------------
class ApplicationKidneyStone(appl.Application):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        super(ApplicationKidneyStone, self).__init__(com)

        self.CTNumberThreshold = 200
        self.CTNumberAir = -1000

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ApplyModel(self, lowHighDicomDir, fileName, outputImageName, methodObj, **kwargs):
        applyDenoise = kwargs.get("applyDenoise", True) # by default, apply denoise to Zeff map after it is calculated
        zeffThreshold = kwargs.get("zeffThreshold", 7.0)

        enh = enhance.Enhance()

        lowKVPFileName  = ""
        highKVPFileName = ""
        dd = dicom_decoder.dect_dicom_decoder()

        dd.ProcessDirectory(lowHighDicomDir)

        lowOrHighKVP = dd.CheckLowOrHighKVP(fileName)

        if lowOrHighKVP == "low":
            lowKVPFileName  = fileName
            highKVPFileName = dd.FindPairingFile(fileName)
        elif lowOrHighKVP == "high":
            highKVPFileName = fileName
            lowKVPFileName  = dd.FindPairingFile(fileName)

        df = pydicom.dcmread(lowKVPFileName)
        img = df.pixel_array
        plt.imsave("low_image.png", img, cmap = "Greys_r")

        df = pydicom.dcmread(highKVPFileName)
        img = df.pixel_array
        plt.imsave("high_image.png", img, cmap = "Greys_r")

        imageLowKVP  = self.com.ProcessCTImage(lowKVPFileName)
        imageHighKVP = self.com.ProcessCTImage(highKVPFileName)

        imageZeff = self.com.CalculateGivenCTNumber(imageHighKVP, imageLowKVP, methodObj)

        if applyDenoise:
            imageZeff = enh.Denoise(imageZeff, filterSize = 3)

        imageZeff = self.ApplyMaskToZeff(imageZeff, imageHighKVP)
        imageZeff = self.RemoveSmallIsland(imageZeff)

        np.save("image_high_kvp", imageHighKVP)
        np.save("image_low_kvp", imageLowKVP)
        np.save("image_zeff", imageZeff)

        self.imageRow = 1
        self.imageCol = 1
        self.imageIdx = 0
        fig = plt.figure(figsize = (7, 6))
        #plt.subplots_adjust(wspace = 0.01, hspace = 0.1, right = 0.9)

        #self.AddImageRegular(fig, imageLowKVP, "Low kVp CT image")

        zeffMin = np.amin(imageZeff)
        zeffMax = np.amax(imageZeff)
        msg = "--> after masking and small island removal: Zeff min = {:f}, max = {:f}".format(zeffMin, zeffMax)
        print(msg)
        self.AddImageZeff(fig, imageZeff, imageLowKVP, zeffMin, zeffMax, zeffThreshold, "Zeff")

        plt.savefig(outputImageName, bbox_inches = 'tight')

        plt.savefig(outputImageName.replace(".pdf", ".png"), bbox_inches='tight', pad_inches = 0)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ApplyMaskToZeff(self, myZeff, myImage):
        zeffMin = np.amin(myZeff)
        zeffMax = np.amax(myZeff)
        msg = "--> before masking: Zeff min = {:f}, max = {:f}".format(zeffMin, zeffMax)
        print(msg)

        myZeff = np.where(myImage < self.CTNumberThreshold, 0, myZeff)

        zeffMin = np.amin(myZeff)
        zeffMax = np.amax(myZeff)
        msg = "--> after masking: Zeff min = {:f}, max = {:f}".format(zeffMin, zeffMax)
        print(msg)

        return myZeff

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CreateNewColormap(self, zeffMinValue, zeffMaxValue, zeffThreshold):
        numBin = 512

        delta = (zeffMaxValue - zeffMinValue) / numBin
        temp = int((zeffThreshold - zeffMinValue) / delta)

        actualZeffThreshold = zeffMinValue + delta * temp
        msg = ">>> actual Zeff threshold = {:f}".format(actualZeffThreshold)
        print(msg)

        colorMapTemplate = cm.get_cmap('bwr_r', numBin)
        newcolors = colorMapTemplate(np.linspace(0, 1, numBin))

        newcolors[0 : 1, :] = np.array([0.0, 0.0, 0.0, 0.6])
        newcolors[1 : temp, :] = np.array([1.0, 0.0, 0.0, 0.6])
        newcolors[temp : -1, :] = np.array([0.0, 0.0, 1.0, 0.6])

        newcmp = mc.ListedColormap(newcolors)

        return newcmp

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddImageZeff(self, fig, imageZeff, imageBase, zeffMin, zeffMax, zeffThreshold, title):
        self.imageIdx += 1

        newcmp = self.CreateNewColormap(zeffMin, zeffMax, zeffThreshold)
        ax = fig.add_subplot(self.imageRow, self.imageCol, self.imageIdx)

        ax.imshow(imageBase, interpolation = 'none', cmap = cm.gray, vmin = -150, vmax = 400)
        ax.imshow(imageZeff, interpolation = 'none', cmap = newcmp)
        #ax.set_title("({0:d}) {1:s}".format(self.imageIdx, title), fontsize = 12)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def PlotReferenceColorImage(self, filePath):
        df = pydicom.dcmread(filePath)
        img = df.pixel_array
        plt.imsave("reference_color_image.png", img)

        print("--> reference dicom size: ", df.pixel_array.shape)

        fig = plt.figure(figsize = (7, 6))
        ax = fig.add_subplot(1, 1, 1)
        ax.imshow(df.pixel_array, interpolation = 'none')
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        plt.savefig("reference_color_image.pdf", bbox_inches = 'tight')

    #------------------------------------------------------------
    #------------------------------------------------------------
    def RemoveSmallIsland(self, image):
        # zeff is float64
        # returned oldMask is int32
        oldMask = np.where(image != 0.0, 1, 0)
        oldMask = oldMask.astype(np.int8)

        # find all your connected components (white blobs in your image)
        nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(oldMask, connectivity = 8)

        # connectedComponentswithStats yields every seperated component with information on each of them, such as size
        # the following part is just taking out the background which is also considered a component, but most of the time we don't want that.
        sizes = stats[1:, -1]
        nb_components = nb_components - 1

        # minimum size of particles we want to keep (number of pixels)
        # here, it's a fixed value, but you can set it as you want, eg the mean of the sizes or whatever
        min_size = 2

        # your answer image
        newMask = np.zeros((output.shape))
        # for every component in the image, you keep it only if it's above min_size
        for i in range(0, nb_components):
            if sizes[i] >= min_size:
                newMask[output == i + 1] = 1

        myResult = np.where(newMask == 1, image, 0.0)

        return myResult