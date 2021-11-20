import os
import sys

# you need to set duoModulePath to the actual location of duo
duoModulePath = ".."
sys.path.append(duoModulePath)

import duo.zeff.common as common
import duo.zeff.abbema as abbema
import duo.core.element_table as element_table
import duo.core.material as material
import duo.zeff.power_law as power_law
import duo.zeff.taylor as taylor
import duo.zeff.bourque as bourque
import duo.zeff.torikoshi as torikoshi
import numpy as np

if __name__ == '__main__':
    # you need to set dataPath to the actual location of duo's data directory
    dataPath = os.path.join("..", "data")
    et = element_table.ElementTable(dataPath)

    # create materials
    waterMat = material.Material("Water", et) # material name is "Water"
    waterMat.AddElement(1, weightFraction = 0.111894) # Z = 1, i.e. hydrogen
    waterMat.AddElement(8, weightFraction = 0.888106)
    waterMat.Commit()
    waterMat.density = 1.0

    airMat = material.Material("Dry Air", et)
    airMat.AddElement(6  , weightFraction = 0.000124)
    airMat.AddElement(7  , weightFraction = 0.755267)
    airMat.AddElement(8  , weightFraction = 0.231781)
    airMat.AddElement(18 , weightFraction = 0.012827)
    airMat.Commit()
    airMat.density = 1.20479e-03

    materialList = [waterMat, airMat]

    # choose an energy value
    Ehigh = 69.28 # keV
    Elow  = 51.93 # keV
    Eave = (Ehigh + Elow) / 2.0

    # calculate Zeff using different formalisms under Eave
    com = common.Common(dataPath)
    for material in materialList:
        print("\n--> {:s} at {:f} keV".format(material.name, Eave))

        # Mayneord
        md = power_law.Mayneord()
        result = md.CalculateZeff(material)
        print("    {:30s}: {:.2f}".format("Mayneord", result))

        # Bourque
        bq = bourque.Bourque(com)
        result = bq.CalculateZeffAtE(material, Eave)
        print("    {:30s}: {:.2f}".format("Bourque", result))

        # Taylor
        tl = taylor.Taylor(com)
        result = tl.CalculateZeffAtE(material, Eave)
        print("    {:30s}: {:.2f}".format("Taylor", result))

        # Torikoshi
        tk = torikoshi.Torikoshi(com)
        result = tk.CalculateZeffAtE(material, Eave)
        print("    {:30s}: {:.2f}".format("Torikoshi", result))

