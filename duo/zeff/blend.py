import numpy as np
import scipy
import scipy.ndimage
import matplotlib
import matplotlib.colors
import duo.core.duo_exception as de
import cv2 as cv

#------------------------------------------------------------
#------------------------------------------------------------
class Blend:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, allowedZMin, allowedZMax, CTNumber_Elow, CTNumber_Ehigh):
        self.allowedZMin = allowedZMin
        self.allowedZMax = allowedZMax
        self.CTNumber_Elow  = CTNumber_Elow
        self.CTNumber_Ehigh = CTNumber_Ehigh

    #------------------------------------------------------------
    # blend images without using alpha
    # on the mixed HU image, replace iodine solution pixels (from
    # the Zeff image) with air pixels
    #------------------------------------------------------------
    def BlendWithoutAlpha(self, imageMixed, imageZeff):
        threshold = 6.0
        imageZeffCopy = np.copy(imageZeff)
        imageBlend = np.copy(imageMixed)

        # clamp imageZeffCopy
        imageZeffCopy = np.where(imageZeffCopy < self.allowedZMin, self.allowedZMin, imageZeffCopy)
        imageZeffCopy = np.where(imageZeffCopy > self.allowedZMax, self.allowedZMax, imageZeffCopy)

        # denoise
        imageZeffCopy = scipy.ndimage.median_filter(imageZeffCopy, size = 5)

        # modify imageMixed
        idx = np.where(imageZeffCopy > threshold)
        imageBlend[idx] = self.CTNumber_Elow

        return imageBlend

    #------------------------------------------------------------
    # blend images using a constant alpha value
    #------------------------------------------------------------
    def BlendWithConstantAlpha(self, imageMixed, imageZeff):
        threshold = 6.0
        imageZeffCopy = np.copy(imageZeff)

        # initialize alpha
        alpha = np.ones_like(imageZeffCopy)

        # clamp imageZeffCopy
        imageZeffCopy = np.where(imageZeffCopy < self.allowedZMin, self.allowedZMin, imageZeffCopy)
        imageZeffCopy = np.where(imageZeffCopy > self.allowedZMax, self.allowedZMax, imageZeffCopy)

        # denoise
        enh = enhance.Enhance()
        imageZeffCopy = enh.Denoise(imageZeffCopy)

        # modify imageMixed
        idx = np.where(imageZeffCopy > threshold)

        # set alpha
        alpha[idx] = 0.0

        # normalize imageMixed to (0.0, 1.0)
        obj = matplotlib.colors.Normalize(vmin = np.amin(imageMixed), vmax = np.amax(imageMixed), clip = False)
        temp = obj(imageMixed)

        # set up rgba image
        imageBlend = np.zeros((imageMixed.shape[0], imageMixed.shape[1], 4))
        imageBlend[:, :, 0] = temp # r
        imageBlend[:, :, 1] = temp # g
        imageBlend[:, :, 2] = temp # b
        imageBlend[:, :, 3] = alpha # a

        return imageBlend

    #------------------------------------------------------------
    # blend images using a varying alpha value
    #------------------------------------------------------------
    def BlendWithAlpha(self, imageHighKVP, imageLowKVP, imageMixed, imageZeff, method = 0):
        #------------------------------
        # Zeff
        #------------------------------
        thresholdZeff = 4.32
        imageZeffCopy = np.copy(imageZeff)

        # clamp imageZeffCopy
        imageZeffCopy = np.where(imageZeffCopy < self.allowedZMin, self.allowedZMin, imageZeffCopy)
        imageZeffCopy = np.where(imageZeffCopy > self.allowedZMax, self.allowedZMax, imageZeffCopy)

        # denoise
        # imageZeffCopy = scipy.ndimage.median_filter(imageZeffCopy, size = 4)




        #------------------------------
        # HU difference
        #------------------------------
        thresholdHuDiff = 120.0
        imageHuDiff = imageLowKVP - imageHighKVP
        allowedHuDiffMin = 0
        allowedHuDiffMax = np.amax(imageHuDiff)







        #------------------------------
        # HU gradient
        #------------------------------
        gradY = cv.Sobel(imageMixed, cv.CV_64F, 0, 1, ksize = 5)
        thresholdHuGrad = (np.amax(gradY) - np.amin(gradY)) * 0.8 + np.amin(gradY)





        #------------------------------
        # alpha
        #------------------------------
        # initialize alpha
        alpha = np.ones_like(imageZeffCopy)

        # modify imageMixed
        temp = np.logical_and(imageHuDiff > thresholdHuDiff, imageZeffCopy > thresholdZeff)
        idx = np.where(temp)

        # set alpha if Zeff is [threshold, self.allowedZMax]
        # alphaMin for Zeff = self.allowedZMax
        # alphaMax for Zeff = threshold
        alphaMin = 0.0
        alphaMax = 1.0

        fZeff = (alphaMin - alphaMax) / (self.allowedZMax - thresholdZeff) * (imageZeffCopy[idx] - thresholdZeff) + alphaMax
        fHuDiff = (alphaMin - alphaMax) / (allowedHuDiffMax - thresholdHuDiff) * (imageHuDiff[idx] - thresholdHuDiff) + alphaMax

        print("fZeff: max = {0:f}, min = {1:f}".format(np.amax(fZeff), np.amin(fZeff)))
        print("fHuDiff: max = {0:f}, min = {1:f}".format(np.amax(fHuDiff), np.amin(fHuDiff)))

        if method == 0:
            alpha[idx] = 0.0
        elif method == 1:
            alpha[idx] = (fZeff + fHuDiff) / 2
        elif method == 2:
            alpha[idx] = (np.power(fZeff, 4.0) + np.power(fHuDiff, 4.0)) / 2

        idx2 = np.where(gradY > thresholdHuGrad)
        alpha[idx2] = 0.0


        # normalize imageMixed to (0.0, 1.0)
        obj = matplotlib.colors.Normalize(vmin = np.amin(imageMixed), vmax = np.amax(imageMixed), clip = False)
        temp = obj(imageMixed)

        # set up rgba image
        imageBlend = np.zeros((imageMixed.shape[0], imageMixed.shape[1], 4))
        imageBlend[:, :, 0] = temp # r
        imageBlend[:, :, 1] = temp # g
        imageBlend[:, :, 2] = temp # b
        imageBlend[:, :, 3] = alpha # a

        return imageBlend

