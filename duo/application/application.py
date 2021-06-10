import os
import sys
result = os.path.join("..", "..", "..")
sys.path.append(result)
print("--> sys.path = ", sys.path)
import duo.zeff.bourque as bourque
import duo.zeff.taylor as taylor
import duo.zeff.cai as cai
import duo.core.dicom_decoder as dicom_decoder
import duo.core.element_table as element_table
import duo.core.material as material
import numpy as np
import matplotlib
matplotlib.use('PDF')
from matplotlib import cm
import matplotlib.pyplot as plt
import duo.zeff.common as common
import duo.application.kidney_stone as kidney_stone

#------------------------------------------------------------
#------------------------------------------------------------
class Application:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        self.com = com

        self.imageRow = 0
        self.imageCol = 0
        self.imageIdx = 0

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddImageRegular(self, fig, imageRegular, title):
        self.imageIdx += 1

        ax = fig.add_subplot(self.imageRow, self.imageCol, self.imageIdx)
        m = ax.imshow(imageRegular, interpolation = 'none', cmap = cm.gray, vmin = -1000, vmax = 1800)
        ax.set_title("({0:d}) {1:s}".format(self.imageIdx, title), fontsize = 12)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddImageZeff(self, fig, imageZeff, title):
        self.imageIdx += 1

        ax = fig.add_subplot(self.imageRow, self.imageCol, self.imageIdx)
        m = ax.imshow(imageZeff, interpolation = 'none', cmap = cm.viridis, vmin = 0, vmax = 20)
        ax.set_title("({0:d}) {1:s}".format(self.imageIdx, title), fontsize = 12)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        cax = None
        if self.imageCol == 1:
            cax = fig.add_axes([0.9, 0.1, 0.02, 0.8]) # left, bottom, width, height
        elif self.imageCol == 2:
            cax = fig.add_axes([0.9, 0.5, 0.02, 0.8]) # left, bottom, width, height

        cbar = fig.colorbar(m)
        cbar.ax.locator_params(nbins = 12) # nbins is the max number of bins