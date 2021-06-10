import os
import sys
import numpy as np
import duo.core.element_table as element_table
import duo.core.material as material
import duo.core.duo_exception as de
import duo.zeff.common as common

#------------------------------------------------------------
#------------------------------------------------------------
class Nist:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com):
        self.dataPath = com.dataPath
        self.elementTable = com.elementTable
        self.materialList = []

        self.ImportNistMaterial()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ImportNistMaterial(self):
        path = os.path.join(self.dataPath, "nist", "material.txt")
        with open(path, "r") as file:
            mat = None
            for line in file:
                lineString = line.split()

                if line.find("name:") != -1: # if line contains "name:"
                    name = line.replace("name:", "")
                    name = name.lstrip()
                    name = name.rstrip("\n")

                    mat = material.Material(name, self.elementTable)

                elif line.find(":") != -1 and line.find("name") == -1: # if line contains ":" but not "name"
                    lineString = line.split(":")
                    Z = int(lineString[0])
                    wf = float(lineString[1])
                    mat.AddElement(Z, weightFraction = wf)

                elif line.find("end") != -1:
                    self.materialList.append(mat)

                elif len(lineString) == 3:
                    mat.density = float(lineString[2])

            for mat in self.materialList:
                mat.Commit()

        print("--> {0:d} nist materials have been imported.".format(len(self.materialList)))

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Show(self):
        for mat in self.materialList:
            mat.Show()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def SearchMaterialFromNist(self, name):
        material = None
        for mat in self.materialList:
            if mat.name.find(name) != -1:
                material = mat
                break
        if None == material:
            raise de.DuoException("--> Material not found.")
        return material


