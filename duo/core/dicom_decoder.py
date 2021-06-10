import pydicom
import os
import duo.core.duo_exception as de
import numpy as np
import duo.core.timer as timer

#------------------------------------------------------------
#------------------------------------------------------------
class DicomFile:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        self.fullFileName  = None
        self.PatientID     = None
        self.StudyID       = None
        self.SeriesNumber  = None
        self.KVP           = None
        self.SliceLocation = None
        self.ImageType     = None

        self.dicomData     = None

#------------------------------------------------------------
#------------------------------------------------------------
class DECTDicomPair:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        self.dicomFileLowKVP  = None
        self.dicomFileHighKVP = None
        self.dicomFileMixed   = None

#------------------------------------------------------------
# dual energy CT DICOM decoder
#------------------------------------------------------------
class dect_dicom_decoder:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        self.dicomDir = ""
        self.dicomFileList = []
        self.dectDicomDict = {} # (series number low kVp, series number high kVp) : list of DECTDicomPair data

        self.lowKVP = 80.0
        self.lowKVP100 = 100.0
        self.highKVP = 140.0

    #------------------------------------------------------------
    # Initialize self.dectDicomDict member
    # This function assumes that in the given directory dicomDir,
    # --- the data is for one patient only
    # --- there are even number of series
    # --- for a low kVp image, there exists a pairing high kVp image, and vice versa
    # If either one does not hold, exception will be thrown.
    # This function does not assume the existence of mixed images nor handle them.
    #------------------------------------------------------------
    def ProcessDirectory(self, dicomDir):
        print("--> Process directory")

        tmg = timer.TimerManager()
        tmg.StartOrResume("ProcessDirectory")

        self.dicomDir = dicomDir
        fileNameList = sorted(os.listdir(self.dicomDir))

        # get basic info
        for i in range(len(fileNameList)):
            fullFileName = os.path.join(self.dicomDir, fileNameList[i])
            df = pydicom.dcmread(fullFileName, specific_tags = ["PatientID", "StudyID", "SeriesNumber", "SliceLocation", "KVP", "ImageType"])

            dif = DicomFile()
            dif.fullFileName   = fullFileName
            dif.PatientID      = df.PatientID
            dif.StudyID        = df.StudyID
            dif.SeriesNumber   = df.SeriesNumber
            dif.KVP            = df.KVP
            dif.SliceLocation  = df.SliceLocation
            dif.ImageType      = df.ImageType

            # skip mixed images and only consider original images
            # but among the test data, some mixed images are somehow also marked as ORIGINAL
            if dif.ImageType[0] == 'ORIGINAL':
                self.dicomFileList.append(dif)

        # check if patient is unique
        temp = self.dicomFileList[0].PatientID
        for item in self.dicomFileList:
            if item.PatientID != temp:
                raise de.DuoException("--> More than one patient is found.")

        # in some test data, study id is empty, so it is not checked

        # get series number
        seriesNumberList = []
        for item in self.dicomFileList:
            if item.SeriesNumber not in seriesNumberList:
                seriesNumberList.append(item.SeriesNumber)
        if len(seriesNumberList) % 2 != 0:
            raise de.DuoException("--> Number of series (excluding mixed images) is not even.")
        print("    Number of series (excluding mixed images) : {0:d}".format(len(seriesNumberList)))

        seriesNumberFileListDictLowKVP = {} # series number : file list
        seriesNumberFileListDictHighKVP = {} # series number : file list
        for seriesNumber in seriesNumberList:
            fileList = []
            for item in self.dicomFileList:
                if item.SeriesNumber == seriesNumber:
                    fileList.append(item)

            # sort according to slice location
            fileList.sort(key = lambda file : file.SliceLocation)

            if isinstance(fileList[0].KVP, float):
                # dirty fix: both 80 and 100 kVp are allowed for low kVp
                if np.fabs(fileList[0].KVP - self.lowKVP) < 1e-8 or\
                   np.fabs(fileList[0].KVP - self.lowKVP100) < 1e-8:
                    seriesNumberFileListDictLowKVP[seriesNumber] = fileList
                elif np.fabs(fileList[0].KVP - self.highKVP) < 1e-8:
                    seriesNumberFileListDictHighKVP[seriesNumber] = fileList
                else:
                    raise de.DuoException("--> Unknown kVp.")

        # pair low kVp and high kVp dicom files
        # pretty messy method
        pairList = []
        for keyLow, valueLow in seriesNumberFileListDictLowKVP.items():
            paired = False

            for keyHigh, valueHigh in seriesNumberFileListDictHighKVP.items():
                if len(valueLow) != len(valueHigh):
                    paired = False
                    continue
                else:
                    paired = True
                    for i in range(len(valueLow)):
                        if valueLow[i].SliceLocation != valueHigh[i].SliceLocation:
                            paired = False
                            break

                if paired == True:
                    pairList.append((keyLow, keyHigh))
                    break

            if paired == False:
                raise de.DuoException("--> Pairing files not found.")

        # generate result: self.dectDicomDict
        for pair in pairList:
            valueLow  = seriesNumberFileListDictLowKVP[pair[0]]
            valueHigh = seriesNumberFileListDictHighKVP[pair[1]]

            dectDicomPairList = []
            for i in range(len(valueLow)):
                dectDicomPair = DECTDicomPair()
                dectDicomPair.dicomFileLowKVP  = valueLow[i]
                dectDicomPair.dicomFileHighKVP = valueHigh[i]
                dectDicomPairList.append(dectDicomPair)

            self.dectDicomDict[pair] = dectDicomPairList

        # validate dicom files
        for key, value in self.dectDicomDict.items():
            for item in value:
                if (item.dicomFileLowKVP.KVP != self.lowKVP and item.dicomFileLowKVP.KVP != self.lowKVP100) or\
                   item.dicomFileHighKVP.KVP != self.highKVP or\
                   item.dicomFileLowKVP.SliceLocation != item.dicomFileHighKVP.SliceLocation:
                    raise de.DuoException("--> Invalid data.")

        tmg.Stop("ProcessDirectory")
        result = tmg.GetElapsedTimeInSecond("ProcessDirectory")
        print("    time = {0:f} [s]".format(result))


    #------------------------------------------------------------
    # Given a file name, find its pairing file.
    # self.ProcessDirectory() must be called beforehand.
    #------------------------------------------------------------
    def FindPairingFile(self, givenFileName):
        for key, value in self.dectDicomDict.items():
            for item in value:
                if item.dicomFileLowKVP.fullFileName == givenFileName:
                    return item.dicomFileHighKVP.fullFileName
                elif item.dicomFileHighKVP.fullFileName == givenFileName:
                    return item.dicomFileLowKVP.fullFileName

        # at this point, nothing is found, throw an error
        raise de.DuoException("--> File not found.")

    #------------------------------------------------------------
    # Given a file name, check if it is low kVp or high kVp.
    # self.ProcessDirectory() must be called beforehand.
    #------------------------------------------------------------
    def CheckLowOrHighKVP(self, givenFileName):
        for key, value in self.dectDicomDict.items():
            for item in value:
                if item.dicomFileLowKVP.fullFileName == givenFileName:
                    return "low"
                elif item.dicomFileHighKVP.fullFileName == givenFileName:
                    return "high"

        # at this point, nothing is found, throw an error
        raise de.DuoException("--> File not found.")

    #------------------------------------------------------------
    # Given mixedImageDicomDir, find the mixed images and add the file objects to member self.dectDicomDict.
    # self.ProcessDirectory() must be called beforehand.
    # This function assumes that the mixed images are stored in a directory separate from low/high kVp images.
    #------------------------------------------------------------
    def ProcessDirectoryWithMixedImage(self, mixedImageDicomDir):
        fileNameList = sorted(os.listdir(mixedImageDicomDir))

        # get basic info
        dicomFileList = []
        for i in range(len(fileNameList)):
            fullFileName = os.path.join(mixedImageDicomDir, fileNameList[i])
            df = pydicom.dcmread(fullFileName, specific_tags = ["SeriesNumber", "SliceLocation"])

            dif = DicomFile()
            dif.fullFileName   = fullFileName
            dif.SeriesNumber   = df.SeriesNumber
            dif.SliceLocation  = df.SliceLocation

            dicomFileList.append(dif)

        # get series number
        seriesNumberList = []
        for item in dicomFileList:
            if item.SeriesNumber not in seriesNumberList:
                seriesNumberList.append(item.SeriesNumber)

        seriesNumberFileListDict = {} # series number : file list
        for seriesNumber in seriesNumberList:
            mixedImageFileList = []
            for item in dicomFileList:
                if item.SeriesNumber == seriesNumber:
                    mixedImageFileList.append(item)
            seriesNumberFileListDict[seriesNumber] = mixedImageFileList

        # find mixed image for self.dectDicomDict
        for seriesNumber, mixedImageFileList in seriesNumberFileListDict.items():
            for seriesPair, dectDicomPairList in self.dectDicomDict.items():
                if len(mixedImageFileList) == len(dectDicomPairList):
                    for itemMixed in mixedImageFileList:
                        for i in range(len(dectDicomPairList)):
                            if itemMixed.SliceLocation == dectDicomPairList[i].dicomFileLowKVP.SliceLocation:
                                self.dectDicomDict[seriesPair][i].dicomFileMixed = itemMixed

        # sanity check
        for key, value in self.dectDicomDict.items():
            for item in value:
                if item.dicomFileMixed == None:
                    raise de.DuoException("--> Mixed image not found.")

    #------------------------------------------------------------
    # Given a file name, find its mixed image file.
    # self.ProcessDirectory() and self.ProcessDirectoryWithMixedImage() must be called beforehand.
    #------------------------------------------------------------
    def FindMixedFile(self, givenFileName):
        for key, value in self.dectDicomDict.items():
            for item in value:
                if item.dicomFileLowKVP.fullFileName == givenFileName or\
                   item.dicomFileHighKVP.fullFileName == givenFileName:
                    return item.dicomFileMixed.fullFileName

        # at this point, nothing is found, throw an error
        raise de.DuoException("--> File not found.")

    #------------------------------------------------------------
    # Give 2 directories (one with low and high images, the other with mixed images)
    # pair them (i.e. initializing the DECTDicomPair object)
    #------------------------------------------------------------
    def ProcessDirectoryWithLowHighAndMixedImages(self, dicomDirLowAndHigh, mixedImageDicomDir):
        self.ProcessDirectory(dicomDirLowAndHigh)

        self.ProcessDirectoryWithMixedImage(mixedImageDicomDir)


