import numpy as np
import scipy
import scipy.ndimage
import matplotlib
import matplotlib.colors
import cv2 as cv
import duo.core.duo_exception as de

#------------------------------------------------------------
#------------------------------------------------------------
class Enhance:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        pass

    #------------------------------------------------------------
    #------------------------------------------------------------
    def RemoveVerticalPartialVolumeEffect(self, image):
        imageCopy = np.copy(image)

        # gradY = cv.Laplacian(imageCopy, cv.CV_64F)
        gradY = cv.Sobel(imageCopy, cv.CV_64F, 0, 1, ksize = 5)

        threshold = (np.amax(gradY) - np.amin(gradY)) * 0.7 + np.amin(gradY)

        imageCopy = np.where(gradY > threshold, np.amin(imageCopy), imageCopy)

        return imageCopy

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Denoise(self, image, filterSize = 5):
        imageCopy = np.copy(image)

        imageCopy = scipy.ndimage.median_filter(imageCopy, size = filterSize)

        return imageCopy


