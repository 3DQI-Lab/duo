# ideally
#     HU
#         air (-999 at 80kVp or 120 kVp)
#         soft tissue (57 at 80kVp, 53 at 120 kVp)
#         20 mg/cc iodine water (997 at 80kVp, 528 at 120 kVp)
#     Zeff
#         soft tissue (3.7)
#         20 mg/cc iodine water (6)
#         air (7.4)
#
# in out test
#     HU is in agreement with ideal case
#     Zeff is not in agreement with ideal case
#         soft tissue (3~4)
#         air (7~8)
#         iodine water (>11)
#
# we want
#     reduce iodine water's HU to a point way below soft tissue's
#
# intraluminal remains
#     refers to stool that needs to be separate from soft tissue


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
import duo.zeff.zeff_calculator as zeff_calculator
import duo.zeff.blend as blend
import duo.zeff.enhance as enhance
import duo.core.surface_fitting as sf
import scipy.interpolate as sip
from matplotlib import cm
import scipy
import scipy.ndimage
import pydicom
import duo.core.timer as timer
import matplotlib.colors
import cv2 as cv
import duo.zeff.common as common
import duo.zeff.control_point as control_point


#------------------------------------------------------------
# Given high and low CT images, calculate Zeff map
#------------------------------------------------------------
class Cai:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method = 'bourque', controlPointType = "colonECExtra"):
        self.com = com
        self.method = method
        self.elementTable = com.elementTable

        self.cp = None
        if controlPointType == "colonEC":
            self.cp = control_point.ControlPointManagerColonEC(com, self.method, isToDumpIntermediateData = True)
        elif controlPointType == "kidneyStone":
            self.cp = control_point.ControlPointManagerKidneyStone(com, self.method, isToDumpIntermediateData = True)
        elif controlPointType == "special1":
            self.cp = control_point.ControlPointManagerSpecial1(com, self.method, isToDumpIntermediateData = True)
        elif controlPointType == "special2":
            self.cp = control_point.ControlPointManagerSpecial2(com, self.method, isToDumpIntermediateData = True)
        elif controlPointType == "colonECExtra":
            self.cp = control_point.ControlPointManagerColonECExtra(com, self.method, isToDumpIntermediateData = True)
        elif controlPointType == "kidneyStoneExtra":
            self.cp = control_point.ControlPointManagerKidneyStoneExtra(com, self.method, isToDumpIntermediateData = True)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateGivenCTNumber(self, imageHighKVP, imageLowKVP):
        return self.Calculate(imageHighKVP, imageLowKVP)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Calculate(self, imageHighKVP, imageLowKVP):
        imageZeff = self.cp.rbf(imageHighKVP, imageLowKVP)
        return imageZeff

    #------------------------------------------------------------
    # first yield totalNumImage
    # then yield imageIdx
    #------------------------------------------------------------
    def ApplyModelRecursively(self, dicomDir, outputImageDir):
        tmg = timer.TimerManager()
        tmg.StartOrResume("total")

        dd = dicom_decoder.dect_dicom_decoder()
        dd.ProcessDirectory(dicomDir)

        # get total number of dicom files
        totalNumImage = 0
        for key, value in dd.dectDicomDict.items():
            totalNumImage += len(value)
        yield totalNumImage

        if not os.path.exists(outputImageDir):
            os.makedirs(outputImageDir)

        # generate log file
        fileName = os.path.join(outputImageDir, "log.txt")
        with open(fileName, "w") as outFile:
            msg = "--> This file is generated by duo utility. It lists the parent files (high and low kVp DICOM files) for each figure.\n\n\n\n"
            outFile.write(msg)

            imageIdx = 0
            for key, value in dd.dectDicomDict.items():
                msg = "--> Low kVp series: {0:d}\n    High kVp series: {1:d}\n".format(key[0], key[1])
                outFile.write(msg)
                print(msg)

                msg = "{0:<20s}{1:<120s}{2:<120s}{3:>40s}{4:>40s}\n".format("Image name",
                                                                            "Low kVp DICOM file",
                                                                            "High kVp DICOM file",
                                                                            "Number of pixels out of range",
                                                                            "Percent of pixels out of range")
                outFile.write(msg)
                counterSeries = 0
                for item in value:
                    # calculate image
                    imageLowKVP  = self.ProcessCTImage(item.dicomFileLowKVP.fullFileName)
                    imageHighKVP = self.ProcessCTImage(item.dicomFileHighKVP.fullFileName)

                    imageZeff = self.rbf(imageHighKVP, imageLowKVP)
                    stat = np.where(np.logical_or(imageZeff < self.allowedZMin, imageZeff > self.allowedZMax), 1, 0)
                    sum = np.sum(stat)
                    totalNumPixel = imageZeff.shape[0] * imageZeff.shape[1]
                    ratio = sum / totalNumPixel

                    imageZeff = np.where(imageZeff < self.allowedZMin, self.allowedZMin, imageZeff)
                    imageZeff = np.where(imageZeff > self.allowedZMax, self.allowedZMax, imageZeff)


                    # increment counters
                    counterSeries += 1
                    imageIdx += 1
                    yield imageIdx

                    # write log file
                    imageName = str(key[0]) + "_" + str(key[1]) + "_" + str(counterSeries).zfill(4) + ".png"
                    msg = "{0:<20s}{1:<120s}{2:<120s}{3:>40d}{4:>40.2f}%\n".format(imageName,
                                                                                os.path.basename(item.dicomFileLowKVP.fullFileName),
                                                                                os.path.basename(item.dicomFileHighKVP.fullFileName),
                                                                                sum,
                                                                                ratio * 100.0)
                    outFile.write(msg)
                    print("    Generating image " + imageName)

                    # create image
                    fullImageName = os.path.join(outputImageDir, imageName)

                    # # use imageio package
                    # imageio.imwrite(fullImageName, imageZeff)

                    # use matplotlib
                    fig = plt.figure(figsize = (12, 10))
                    ax3 = fig.add_subplot(111)
                    m3 = ax3.imshow(imageZeff, interpolation = 'none', cmap = cm.viridis)
                    ax3.set_title("Average effective atomic number map", fontsize = 10)
                    cbar3 = fig.colorbar(m3, ax = ax3)
                    cbar3.ax.locator_params(nbins = 12) # nbins is the max number of bins
                    plt.savefig(fullImageName)

                outFile.write("\n\n\n")

        tmg.Stop("total")
        totalTime = tmg.GetElapsedTimeInSecond("total")
        print("--> total time = {0:f} [s]".format(totalTime))





