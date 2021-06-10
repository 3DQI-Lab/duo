import duo.core.element as duoelement
import duo.core.duo_exception as de
import duo.core.material as material
import numpy as np

#------------------------------------------------------------
#------------------------------------------------------------
class MaterialComponent:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, material):
        self.material = material
        self.materialWeightFraction = 0.0
        self.materiaMolarFraction = 0.0

#------------------------------------------------------------
# a mixture must be composed of materials
#------------------------------------------------------------
class Mixture(material.Material):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, name, elementTable):
        # call super class method
        super().__init__(name, elementTable)

        self.materialList = []

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddMaterial(self, otherMaterial, **kwargs):
        if not otherMaterial.isCommitted:
            raise de.DuoException("--> Material not committed.")

        for key, value in kwargs.items():
            if "materialWeightFraction" == key:
                materialWeightFraction = kwargs["materialWeightFraction"]

                for Z, ec in otherMaterial.elementList.items():
                    temp = materialWeightFraction * ec.weightFraction
                    self.AddElement(Z, weightFraction = temp)

                mc = MaterialComponent(otherMaterial)
                mc.materialWeightFraction = materialWeightFraction
            elif "materiaMolarFraction" == key:
                materiaMolarFraction = kwargs["materiaMolarFraction"]

                for Z, ec in otherMaterial.elementList.items():
                    temp = materiaMolarFraction * ec.atomicFraction
                    self.AddElement(Z, atomicFraction = temp)

                mc = MaterialComponent(otherMaterial)
                mc.materiaMolarFraction = materiaMolarFraction
            else:
                raise de.DuoException("--> Invalid argument.")


        self.materialList.append(mc)

    #------------------------------------------------------------
    # note that when two material is mixed together
    # Vmix < V1 + V2, so this function is not very accurate
    #------------------------------------------------------------
    def CalculateMixtureDensity(self):
        mass = 0.0
        volume = 0.0
        for mc in self.materialList:
            mass += mc.materialWeightFraction
            volume += mc.materialWeightFraction / mc.material.density
        self.density = mass / volume

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Commit(self):
        # call super class method
        super().Commit()

        # double-check that for each material, material and atomic weight fractions are normalized
        for mc in self.materialList:
            weightFractionSum = 0.0
            atomicFractionSum = 0.0

            for Z, ec in mc.material.elementList.items():
                weightFractionSum += ec.weightFraction
                atomicFractionSum += ec.atomicFraction

            if np.fabs(weightFractionSum - 1.0) > 1e-6 or\
               np.fabs(atomicFractionSum - 1.0) > 1e-6:
                raise de.DuoException("--> Unnormalized fraction values found.")

        # validate fraction values
        # the fraction values should be either all weight fractions or all molar fraction
        # mix of both is not allowed
        useMaterialWeightFraction = False;
        usemateriaMolarFraction = False;
        for mc in self.materialList:
            if (not useMaterialWeightFraction) and (mc.materialWeightFraction != 0.0):
                useMaterialWeightFraction = True

            if (not usemateriaMolarFraction) and (mc.materiaMolarFraction != 0.0):
                usemateriaMolarFraction = True

        if useMaterialWeightFraction and usemateriaMolarFraction:
            raise de.DuoException("--> Mixing material weight fraction and material molar fraction is not supported.")

        if (not useMaterialWeightFraction) and (not usemateriaMolarFraction):
            raise de.DuoException("--> Wrong material weight fraction or material molar fraction.")

        # convert between weight and molar fraction
        if useMaterialWeightFraction:
            for mc in self.materialList:
                A_all = 0.0
                for Z, ec in mc.material.elementList.items():
                    element = self.elementTable.GetElementByZ(Z)
                    A_all += element.A
                mc.materiaMolarFraction = mc.materialWeightFraction / A_all

        elif usemateriaMolarFraction:
            for mc in self.materialList:
                A_all = 0.0
                for Z, ec in mc.material.elementList.items():
                    element = self.elementTable.GetElementByZ(Z)
                    A_all += element.A
                mc.materialWeightFraction = A_all * mc.materiaMolarFraction

        # normalize fraction values
        weightFractionSum = 0.0
        molarFractionSum = 0.0

        for mc in self.materialList:
            weightFractionSum += mc.materialWeightFraction
            molarFractionSum += mc.materiaMolarFraction

        for mc in self.materialList:
            mc.materialWeightFraction /= weightFractionSum
            mc.materiaMolarFraction /= molarFractionSum






