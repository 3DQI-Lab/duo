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
import cv2 as cv
import string

#------------------------------------------------------------
#------------------------------------------------------------
class Common:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self,
                 dataPath,
                 method = "original",
                 use140kVpWithSn = True,
                 use100kVp = False):
        self.method = method

        self.dataPath = dataPath
        self.elementTable = None

        self.CTNumber_Elow  = -1000
        self.CTNumber_Ehigh = 3095

        # mu/rho
        self.muWater_Ehigh = 0.0
        self.muWater_Elow  = 0.0

        self.muAir_Ehigh = 0.0
        self.muAir_Elow  = 0.0

        # Van Abbema, J.K., Van der Schaaf, A., Kristanto, W.,
        # Groen, J.M. and Greuter, M.J., 2012. Feasibility and
        # accuracy of tissue characterization with dual source
        # computed tomography. Physica Medica, 28(1), pp.25-32.

        # according to https://radiologykey.com/technical-aspects-of-dual-energy-ct-with-dual-source-ct-systems/
        # The SOMATOM Definition Flash uses an additional tin filter (Sn) with a thickness of 0.4 mm to shift the
        # mean energy of the 140 kV spectrum from 69 to 89 keV, see Fig. 2.2. The mean energy of the 80 kV spectrum is 52 keV.
        if use140kVpWithSn: # 140 kVp with Sn (Tin) filter
            self.Ehigh = 89
        else: # without the filter
            self.Ehigh = 69.28

        if use100kVp: # temporary fix
            self.Elow  = 63.292 # keV, 100 kVp
        else: # normal case: 80 kVp
            self.Elow  = 51.93 # keV, 80 kVp

        self.Eave = (self.Ehigh + self.Elow) / 2.0

        self.Initialize()

    #------------------------------------------------------------
    # For some reason, NIST material has two versions
    # original: https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=276
    # alt: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html
    #------------------------------------------------------------
    def Initialize(self):
        print("--> Python version = ", sys.version)
        print("    matplotlib version = ", matplotlib.__version__)
        print("    numpy version = ", np.version.version)
        print("    scipy version = ", scipy.__version__)

        # load xs data
        self.elementTable = element_table.ElementTable(self.dataPath)

        water = material.Material("Water", self.elementTable)
        if self.method == "original":
            water.AddElement(1 , weightFraction =  0.111894)
            water.AddElement(8 , weightFraction =  0.888106)
            water.Commit()
            water.density = 1.0
        elif self.method == "alt":
            water.AddElement(1 , weightFraction =  0.111898)
            water.AddElement(8 , weightFraction =  0.888102)
            water.Commit()
            water.density = 1.0

        air = material.Material("Air", self.elementTable)
        if self.method == "original":
            air.AddElement(6  , weightFraction = 0.000124)
            air.AddElement(7  , weightFraction = 0.755267)
            air.AddElement(8  , weightFraction = 0.231781)
            air.AddElement(18 , weightFraction = 0.012827)
            air.Commit()
            air.density = 1.20479e-03
        elif self.method == "alt":
            air.AddElement(6  , weightFraction = 0.000124)
            air.AddElement(7  , weightFraction = 0.755268)
            air.AddElement(8  , weightFraction = 0.231781)
            air.AddElement(18 , weightFraction = 0.012827)
            air.Commit()
            air.density = 1.205E-03

        self.muWater_Ehigh = water.CalculateMacAtE(self.Ehigh) * water.density
        self.muWater_Elow  = water.CalculateMacAtE(self.Elow) * water.density
        self.muAir_Ehigh   = air.CalculateMacAtE(self.Ehigh) * air.density
        self.muAir_Elow    = air.CalculateMacAtE(self.Elow) * air.density

        print("    muWater_Ehigh = {0:.16e}".format(self.muWater_Ehigh))
        print("    muWater_Elow  = {0:.16e}".format(self.muWater_Elow ))
        print("    muAir_Ehigh   = {0:.16e}".format(self.muAir_Ehigh  ))
        print("    muAir_Elow    = {0:.16e}".format(self.muAir_Elow   ))

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ProcessCTImage(self, fileName):
        df = pydicom.dcmread(fileName)

        print("--> process CT image")
        print("    RescaleSlope = ", df.RescaleSlope)
        print("    RescaleIntercept = ", df.RescaleIntercept)

        # converts raw pixel data values to a specific (possibly unitless) physical quantity,
        # such as Hounsfield units for CT
        # applies Modality LUT transformation to raw data
        # df.pixel_array may be uint16
        # image is float64
        image = pdutil.apply_modality_lut(df.pixel_array, df)
        print("    Before adjustment: pixel data range = [{0:f} {1:f}]".format(np.amin(image), np.amax(image)))

        print("    df.pixel_array type = {0:s}, shape = {1:s}".format(
            str(df.pixel_array.dtype),
            str(df.pixel_array.shape)))
        print("    image type = {0:s}, shape = {1:s}".format(
            str(image.dtype),
            str(image.shape)))

        # clamp
        # for out-of-range data
        image = np.where(image < self.CTNumber_Elow,  self.CTNumber_Elow, image)
        image = np.where(image > self.CTNumber_Ehigh, self.CTNumber_Ehigh, image)
        print("    After adjustment: pixel data range = [{0:f} {1:f}]".format(np.amin(image), np.amax(image)))

        # clamp
        # for pixels with low HU, set them to air
        image = np.where(image < -800, -1000, image)

        return image


    #------------------------------------------------------------
    #------------------------------------------------------------
    def ConvertCTNumberToMu(self, CTNumber, mu_water, mu_air):
        mu = CTNumber / 1000.0 * (mu_water - mu_air) + mu_water
        return mu

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ConvertMuToCTNumber(self, mu, mu_water, mu_air):
        CTNumber = (mu - mu_water) / (mu_water - mu_air) * 1000.0
        return CTNumber

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateGivenCTNumber(self, imageHighKVP, imageLowKVP, methodObj, useGPU = False):
        imageZeff = None

        if useGPU:
            imageZeff = methodObj.CalculateGivenCTNumberGPU(imageHighKVP, imageLowKVP)
        else:
            imageZeff = methodObj.CalculateGivenCTNumber(imageHighKVP, imageLowKVP)

        return imageZeff

    #------------------------------------------------------------
    # full plot:
    #     highKVP
    #     lowKVP
    #     mixed
    #     Zeff
    #     reference
    #     mixed with EC
    #------------------------------------------------------------
    def ApplyModel(self, lowHighDicomDir, fileName, outputImageName, methodObj, **kwargs):
        msg = "--> Mean high energy = {:f}, mean low energy = {:f}, average = {:f}".format(self.Ehigh, self.Elow, self.Eave)
        print(msg)

        # whether to apply denoising to the images
        # the original highKVP, lowKVP, mixed and reference are always plotted, whereas
        # the denoised Zeff and mixed with EC can be optionally plotted when applyDenoise is True
        applyDenoise = kwargs.get("applyDenoise", True)

        imageMask = kwargs.get("imageMask", None)

        # whether to process the mixed image in the specified dicom dir
        # the mixed image shall be stored in a directory separate from low/high kVp images
        mixedImageDicomDir = kwargs.get("mixedImageDicomDir", None)
        showMixedImage = kwargs.get("showMixedImage", True)

        # whether to apply EC to mixed image
        applyEC = kwargs.get("applyEC", False)

        saveECImage = kwargs.get("saveECImage", False)

        referenceImagePath = kwargs.get("referenceImagePath", None)

        useGPU = kwargs.get("useGPU", False)
        zeffWindowLower = kwargs.get("zeffWindowLower", 7)
        zeffWindowUpper = kwargs.get("zeffWindowUpper", 12)


        tmg = timer.TimerManager()
        tmg.StartOrResume("total")

        lowKVPFileName  = ""
        highKVPFileName = ""
        mixedFileName   = ""
        dd = dicom_decoder.dect_dicom_decoder()

        dd.ProcessDirectory(lowHighDicomDir)

        if mixedImageDicomDir:
            dd.ProcessDirectoryWithMixedImage(mixedImageDicomDir)

        lowOrHighKVP = dd.CheckLowOrHighKVP(fileName)

        if lowOrHighKVP == "low":
            lowKVPFileName  = fileName
            highKVPFileName = dd.FindPairingFile(fileName)
            if mixedImageDicomDir:
                mixedFileName = dd.FindMixedFile(fileName)
        elif lowOrHighKVP == "high":
            highKVPFileName = fileName
            lowKVPFileName  = dd.FindPairingFile(fileName)
            if mixedImageDicomDir:
                mixedFileName = dd.FindMixedFile(fileName)

        enh = enhance.Enhance()

        imageLowKVP  = self.ProcessCTImage(lowKVPFileName)
        imageLowKVPDenoised = enh.Denoise(imageLowKVP)

        imageHighKVP = self.ProcessCTImage(highKVPFileName)
        imageHighKVPDenoised = enh.Denoise(imageHighKVP)

        imageMixed = None
        if mixedImageDicomDir:
            imageMixed = self.ProcessCTImage(mixedFileName)
            imageMixedDenoised = enh.Denoise(imageMixed)

        imageZeff = None
        if applyDenoise:
            imageZeff = self.CalculateGivenCTNumber(imageHighKVPDenoised, imageLowKVPDenoised, methodObj, useGPU)
        else:
            imageZeff = self.CalculateGivenCTNumber(imageHighKVP, imageLowKVP, methodObj, useGPU)

        imageMixed_ec = None
        if applyEC:
            if applyDenoise:
                imageMixed_ec = self.ElectronicCleanse(imageLowKVPDenoised,
                                                       imageHighKVPDenoised,
                                                       imageMixedDenoised,
                                                       imageZeff,
                                                       zeffWindowLower,
                                                       zeffWindowUpper,
                                                       imageMask)
            else:
                imageMixed_ec = self.ElectronicCleanse(imageLowKVP,
                                                       imageHighKVP,
                                                       imageMixed,
                                                       imageZeff,
                                                       zeffWindowLower,
                                                       zeffWindowUpper,
                                                       imageMask)

            print("    ec mixed image type = {0:s}, shape = {1:s}".format(
                str(imageMixed_ec.dtype),
                str(imageMixed_ec.shape)))

            if saveECImage:
                outputName = os.path.basename(mixedFileName)
                outputFullPath = "ec_" + outputName

                self.SaveDicom(mixedFileName,
                               imageMixed_ec,
                               outputFullPath)

        if applyDenoise:
            np.save("image_high_kvp", imageHighKVPDenoised)
            np.save("image_low_kvp", imageLowKVPDenoised)
        else:
            np.save("image_high_kvp", imageHighKVP)
            np.save("image_low_kvp", imageLowKVP)
        np.save("image_zeff", imageZeff)

        referenceImage = None
        if referenceImagePath:
            referenceImage = self.ProcessCTImage(referenceImagePath)

        tmg.Stop("total")
        totalTime = tmg.GetElapsedTimeInSecond("total")
        print("--> total time = {0:f} [s]".format(totalTime))

        totalNumImages = 3

        if mixedImageDicomDir:
            totalNumImages += 1

        if applyEC:
            totalNumImages += 1

        if referenceImagePath:
            totalNumImages += 1



        self.AddImageRegular(imageHighKVP, "high_kvp_{:s}".format(outputImageName))
        self.AddImageRegular(imageLowKVP, "low_kvp_{:s}".format(outputImageName))

        if mixedImageDicomDir and showMixedImage:
            self.AddImageRegular(imageMixed, "mixed_image_{:s}".format(outputImageName))

        self.AddImageZeff(imageZeff, "zeff_image_{:s}".format(outputImageName))

        if referenceImagePath:
            self.AddImageRegular(referenceImage, "reference_{:s}".format(outputImageName))

        if applyEC:
            self.AddImageRegular(imageMixed_ec, "mixed_image_with_ec_{:s}".format(outputImageName))

        imageDict = {"imageZeff" : imageZeff,
                     "imageMixed" : imageMixed,
                     "imageMixed_ec" : imageMixed_ec,
                     "referenceImage" : referenceImage,
                     "imageHighKVP" : imageHighKVP,
                     "imageLowKVP" : imageLowKVP}

        return imageDict

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ElectronicCleanse(self, imageLowKVP,
                                imageHighKVP,
                                imageMixed,
                                imageZeff,
                                zeffWindowLower,
                                zeffWindowUpper,
                                imageMask):

        def AlphaTransform(pixelMixed,
                           pixelZeff,
                           pixelGradY,
                           thresholdHuGrad,
                           HUAir,
                           pixelMask):

            newPixelMixed = pixelMixed

            # do not apply alpha transform outside of mask
            #if pixelMask == 0:
            #    return newPixelMixed

            # if Zeff is too large, set the pixel value to air
            if pixelZeff >= zeffWindowUpper:
                newPixelMixed = HUAir
            # if Zeff is not large, do nothing
            elif pixelZeff <= zeffWindowLower:
                newPixelMixed = pixelMixed
            # if Zeff is between the window, transform it
            else:
                # method 1: linear
                # temp = (pixelZeff - zeffWindowUpper) / (zeffWindowLower - zeffWindowUpper)
                # newPixelMixed = (pixelMixed - HUAir) * temp + HUAir

                # method 2: non-linear
                temp = (pixelZeff - zeffWindowUpper) / (zeffWindowLower - zeffWindowUpper)
                temp = np.power(temp, 2.0)
                newPixelMixed = (pixelMixed - HUAir) * temp + HUAir

            #if pixelGradY >= thresholdHuGrad:
            #    newPixelMixed = HUAir

            return newPixelMixed

        # reference: https://docs.opencv.org/2.4/modules/imgproc/doc/filtering.html?highlight=sobel#sobel
        gradY = cv.Sobel(imageMixed, cv.CV_64F, 0, 1, ksize = 5)
        thresholdHuGrad = np.quantile(gradY, 0.98)

        HUAir = -1000

        AlphaTransform_v = np.vectorize(AlphaTransform)
        imageMixed_ec = AlphaTransform_v(imageMixed,
                                         imageZeff,
                                         gradY,
                                         thresholdHuGrad,
                                         HUAir,
                                         imageMask)

        return imageMixed_ec


    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddImageRegular(self, imageRegular, title):
        fig = plt.figure(figsize = (6, 6))

        ax = fig.add_subplot(1,1,1)
        m = ax.imshow(imageRegular,
                      interpolation = 'none',
                      cmap = cm.gray
                      #vmin = -1000,
                      #vmax = 1800
                      )

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        plt.savefig(title, bbox_inches = 'tight')

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddImageZeff(self, imageZeff, title):
        fig = plt.figure(figsize = (6, 4))
        plt.subplots_adjust(right = 0.9)

        ax = fig.add_subplot(1,1,1)

        m = ax.imshow(imageZeff,
                      interpolation = 'none',
                      cmap = cm.viridis,
                      vmin = 0,
                      vmax = 20
                      )

        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        cbar = fig.colorbar(m, ax = ax)
        cbar.ax.locator_params(nbins = 12) # nbins is the max number of bins

        plt.savefig(title, bbox_inches = 'tight')

    #------------------------------------------------------------
    #------------------------------------------------------------
    def PlotSingleImage(self, imageData, outputImageName, markWrongData = True):
        imageMaxValue = 10.0
        imageMinValue = 0.0

        if markWrongData:
            upperBoundValue = 2.2
            lowerBoundValue = 1.0
            numBin = 512

            delta = (imageMaxValue - imageMinValue) / numBin

            # new colormap used to mark problematic pixels
            viridis = cm.get_cmap('viridis', numBin)
            newcolors = viridis(np.linspace(0, 1, numBin))
            upperBoundColor = np.array([0.5, 1.0, 1.0, 1])
            lowerBoundColor = np.array([1.0, 0.0, 0.0, 1])

            temp = (upperBoundValue - imageMinValue) / delta
            newcolors[int(temp) : numBin, :] = upperBoundColor

            temp = (lowerBoundValue - imageMinValue) / delta
            newcolors[0 : int(temp) + 1, :] = lowerBoundColor

            newcmp = mc.ListedColormap(newcolors)
        else:
            newcmp = cm.viridis

        fig = plt.figure(figsize = (10, 10))
        ax = fig.add_subplot(1, 1, 1)
        m = ax.imshow(imageData, interpolation = 'none', cmap = newcmp, vmin = imageMinValue, vmax = imageMaxValue)
        cbar = fig.colorbar(m, ax = ax)
        cbar.ax.locator_params(nbins = 12) # nbins is the max number of bins

        plt.savefig(outputImageName)

        np.save("image_data", imageData)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ShowImageInfo(self, image):
        print("    Data type = ", image.dtype)
        print("    Max = ", np.amax(image))
        print("    Min = ", np.amin(image))
        print("    Shape = ", image.shape)
        print("    Order = ", np.isfortran(image))

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ApplyModelIteratively(self, lowHighDicomDir, mixedImageDicomDir, outputECDicomDir, methodObj, **kwargs):
        applyDenoise = kwargs.get("applyDenoise", True) # by default, apply denoise
        useGPU = kwargs.get("useGPU", False)
        zeffWindowLower = kwargs.get("zeffWindowLower", 7)
        zeffWindowUpper = kwargs.get("zeffWindowUpper", 12)

        enh = enhance.Enhance()

        if not os.path.exists(outputECDicomDir):
            os.mkdir(outputECDicomDir)

        dd = dicom_decoder.dect_dicom_decoder()
        dd.ProcessDirectoryWithLowHighAndMixedImages(lowHighDicomDir, mixedImageDicomDir)

        tmg = timer.TimerManager()
        tmg.StartOrResume("total")

        # dectDicomInfo is a list of DECTDicomPair objects
        dectDicomInfo = list(dd.dectDicomDict.values())[0]
        for idx in range(len(dectDicomInfo)):
            item = dectDicomInfo[idx]

            print("\n\n\n--> idx = ", idx)
            # print("    {0:s}".format(item.dicomFileLowKVP.fullFileName))
            # print("    {0:s}".format(item.dicomFileHighKVP.fullFileName))
            # print("    {0:s}".format(item.dicomFileMixed.fullFileName))

            imageLowKVP  = self.ProcessCTImage(item.dicomFileLowKVP.fullFileName)
            if applyDenoise:
                imageLowKVP = enh.Denoise(imageLowKVP)

            imageHighKVP = self.ProcessCTImage(item.dicomFileHighKVP.fullFileName)
            if applyDenoise:
                imageHighKVP = enh.Denoise(imageHighKVP)

            imageMixed = self.ProcessCTImage(item.dicomFileMixed.fullFileName)
            if applyDenoise:
                imageMixed = enh.Denoise(imageMixed)

            imageZeff = self.CalculateGivenCTNumber(imageHighKVP,
                                                    imageLowKVP,
                                                    methodObj,
                                                    useGPU)

            imageMixed_ec = self.ElectronicCleanse(imageLowKVP,
                                                   imageHighKVP,
                                                   imageMixed,
                                                   imageZeff,
                                                   zeffWindowLower,
                                                   zeffWindowUpper)

            outputName = os.path.basename(item.dicomFileMixed.fullFileName)
            outputName = "ec_{0:03d}_{1:s}".format(idx, outputName)
            outputFullPath = os.path.join(outputECDicomDir, outputName)

            self.SaveDicom(item.dicomFileMixed.fullFileName,
                           imageMixed_ec,
                           outputFullPath)

        tmg.Stop("total")
        totalTime = tmg.GetElapsedTimeInSecond("total")
        print("--> total time = {0:f} [s]".format(totalTime))

    #------------------------------------------------------------
    #------------------------------------------------------------
    def SaveDicom(self, inputFileName,
                        imageMixed_ec,
                        outputFullPath):
        df = pydicom.dcmread(inputFileName)

        tempImage = (imageMixed_ec - df.RescaleIntercept) / df.RescaleSlope

        tempImage = tempImage.astype(np.uint16)
        df.PixelData = tempImage.tobytes()

        df.save_as(outputFullPath)


    #------------------------------------------------------------
    # Given a known material, calculate its CT numbers at high and low energies
    # to form a one-pixel image. Then apply given method to the one-pixel image
    # to calculate its Z values.
    #------------------------------------------------------------
    def VerifyModel(self, knownMaterial, methodObj):
        muLow  = knownMaterial.CalculateMacAtE(self.Elow) * knownMaterial.density
        muHigh = knownMaterial.CalculateMacAtE(self.Ehigh) * knownMaterial.density

        imageLowKVP = self.ConvertMuToCTNumber(muLow, self.muWater_Elow, self.muAir_Elow)
        imageHighKVP = self.ConvertMuToCTNumber(muHigh, self.muWater_Ehigh, self.muAir_Ehigh)

        imageZeff = self.CalculateGivenCTNumber(imageHighKVP, imageLowKVP, methodObj)

        msg = ">>> material = {:s}, CT low = {:f}, CT high = {:f}, Zeff = {:f}".format(knownMaterial.name,
                                                                                       imageLowKVP,
                                                                                       imageHighKVP,
                                                                                       imageZeff)
        print(msg)

