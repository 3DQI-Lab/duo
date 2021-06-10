import os
import xml.etree.ElementTree as ET
import duo.core.duo_exception as de
import duo.core.element as duoelement
import duo.core.photoatomic_xs as duophotoatomic_xs

#------------------------------------------------------------
# Two formats of the ENDF-B-VIII.0 photoatomic libraries are used
# ---- endfb format: for microscopic cross-section
# ---- gnd format: for atomic weight ratio (AWR)
#------------------------------------------------------------
class PhotoAtomicXSIOManager:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, endfbDir="", gndDir=""):
        self.endfbDir = endfbDir
        self.endfbFileList = []

        self.gndDir = gndDir
        self.gndFileList = []

        self.elementList = {}

    #------------------------------------------------------------
    #------------------------------------------------------------
    def InputData(self):
        self.InputAWRFromEndfb()
        self.InputXSFromGnd()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def InputAWRFromEndfb(self):
        print("--> Input atomic weight ratio")

        # use sorted() to enforce order that is consistent across different file systems
        self.endfbFileList = sorted(os.listdir(self.endfbDir))

        for fileName in self.endfbFileList:
            fullPath = os.path.join(self.endfbDir, fileName)
            with open(fullPath, "r") as infile:
                lineNumber = 0
                for line in infile:
                    lineNumber += 1

                    if lineNumber == 2:
                        lineString = line.split()
                        element = duoelement.Element()
                        element.Z = int(float(lineString[0]) / 1000.0)
                        element.AWR = float(lineString[1])

                        element.GetAFromAWR()

                        self.elementList[element.Z] = element



    #------------------------------------------------------------
    #------------------------------------------------------------
    def OutputAWRToFile(self, outputAWRPath = "awr.txt"):
        print("--> Output atomic weight ratio")
        with open(outputAWRPath, "w") as outfile:
            for Z, element in self.elementList.items():
                result = "{0:10d}, {1:10f}, {2:10f}\n".format(element.Z, element.AWR, element.A)
                outfile.write(result)


    #------------------------------------------------------------
    #------------------------------------------------------------
    def InputXSFromGnd(self):
        print("--> Input total microscopic cross-section")

        self.gndFileList = os.listdir(self.gndDir)

        counter = 0 # for debugging purpose
        # iterate files (elements)
        for fileName in self.gndFileList:
            # counter += 1
            # if counter > 1:
                # break

            fullPath = os.path.join(self.gndDir, fileName)
            tree = ET.parse(fullPath)
            root = tree.getroot()

            element = self.GetBasicInfo(root)

            self.GetTotalXS(root, element)

            self.GetAllOtherXS(root, element)

    #------------------------------------------------------------
    # Add element symbol
    #------------------------------------------------------------
    def GetBasicInfo(self, root):
        # use xpath feature to find specific node
        result = root.findall(".//chemicalElement")

        if len(result) != 1:
            raise de.DuoException("--> More than 1 chemical info.")

        info = result[0]

        Z = int(info.attrib["Z"])
        if Z not in self.elementList.keys():
            raise de.DuoException("--> Element not found.")
        element = self.elementList[Z]
        element.symbol = info.attrib["symbol"]

        return element

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetTotalXS(self, root, element):
        # use xpath feature to find specific node
        result = root.findall(".//crossSectionSum[@ENDF_MT='501']//crossSection//regions1d")

        if len(result) != 1:
            raise de.DuoException("--> More than 1 interesting section found. Need further handling.")

        xsXML = result[0].findall(".//XYs1d/values")

        xsList = []
        for subXML in xsXML:
            stringList = subXML.text.split()

            for i in range(len(stringList) // 2):
                energy = float(stringList[2 * i]) / 1000.0 # eV to keV
                microXS = float(stringList[2 * i + 1]) # barn

                xs = duophotoatomic_xs.PhotoAtomicXS(energy, microXS)
                xsList.append(xs)

        element.xsTable["total_ref"] = xsList

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetAllOtherXS(self, root, element):
        result = root.findall(".//crossSectionSum[@ENDF_MT='501']/summands/*")

        for item in result:
            xpath = item.attrib["href"]
            # replace substring "/reactionSuite" with "./" to make a legitimate xpath
            xpath = xpath.replace("/reactionSuite", "./")

            # get xs name
            result = xpath.replace("/crossSection", "")
            result = root.findall(result)
            xsName = result[0].attrib["label"]

            block = root.findall(xpath)

            # use xpath feature to find specific node
            xsXML = block[0].findall(".//XYs1d/values")

            xsList = []
            for subXML in xsXML:
                stringList = subXML.text.split()

                for i in range(len(stringList) // 2):
                    energy = float(stringList[2 * i]) / 1000.0 # eV to keV
                    microXS = float(stringList[2 * i + 1]) # barn

                    xs = duophotoatomic_xs.PhotoAtomicXS(energy, microXS)
                    xsList.append(xs)

            element.xsTable[xsName] = xsList

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Show(self):
        for key, value in self.elementList.items():
            value.Show()
