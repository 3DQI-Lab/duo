import os
import duo.nuclear_data.photoatomic.photoatomic_xs_io as xs

#------------------------------------------------------------
#------------------------------------------------------------
class ElementTable:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, dataPath):
        self.dataPath = dataPath
        self.elementList = {}

        self.Initialize()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        mg = xs.PhotoAtomicXSIOManager(endfbDir=os.path.join(self.dataPath, "photoatomic_endfb"),
                                       gndDir=os.path.join(self.dataPath, "photoatomic_gnd"))

        mg.InputAWRFromEndfb()
        mg.InputXSFromGnd()

        self.elementList = mg.elementList

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetElementByZ(self, Z):
        return self.elementList[Z]
