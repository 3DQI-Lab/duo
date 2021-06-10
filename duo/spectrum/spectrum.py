import os
import os.path
import numpy as np

#------------------------------------------------------------
#------------------------------------------------------------
class SpectrumManager:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, dataPath):
        self.spectrumPath = os.path.join(dataPath, "spectrum")
        self.energyList = []
        self.intensityList = []
        self.effectiveEnergy = 0.0

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ImportSpectrum(self, filePath):
        self.energyList = []
        self.intensityList = []

        fullFilePath = os.path.join(self.spectrumPath, filePath)
        with open(fullFilePath, 'r') as myFile:
            for line in myFile:
                result = line.split()
                if len(result) == 0: # empty line
                    continue
                elif result[0] == '#': # line starting with # is a comment line
                    continue
                else: # data line
                    energy = float(result[0])
                    intensity = float(result[1])
                    self.energyList.append(energy)
                    self.intensityList.append(intensity)

        # convert to numpy
        self.energyList = np.array(self.energyList)
        self.intensityList = np.array(self.intensityList)

        # statistics
        up = np.sum(np.multiply(self.energyList, self.intensityList))
        down = np.sum(self.intensityList)
        self.effectiveEnergy = up / down
        print("--> effective energy = ", self.effectiveEnergy)









