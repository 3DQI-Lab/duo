import duo.core.element as duoelement
import duo.core.duo_exception as de
import duo.core.constant as constant

#------------------------------------------------------------
#------------------------------------------------------------
class ElementComponent:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, Z):
        self.Z = Z
        self.weightFraction = 0.0
        self.atomicFraction = 0.0

#------------------------------------------------------------
#------------------------------------------------------------
class Material:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, name, elementTable):
        self.name = name
        self.elementList = {}
        self.elementTable = elementTable
        self.density = 0.0
        self.isCommitted = False

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddElement(self, Z, **kwargs):
        # if Z exists in the material, accumulate
        if Z in self.elementList.keys():
            for key, value in kwargs.items():
                if "weightFraction" == key:
                    self.elementList[Z].weightFraction += kwargs["weightFraction"]
                elif "atomicFraction" == key:
                    self.elementList[Z].atomicFraction += kwargs["atomicFraction"]
                else:
                    raise de.DuoException("--> Invalid argument.")
        else:
            ec = ElementComponent(Z)

            for key, value in kwargs.items():
                if "weightFraction" == key:
                    ec.weightFraction = kwargs["weightFraction"]
                elif "atomicFraction" == key:
                    ec.atomicFraction = kwargs["atomicFraction"]
                else:
                    raise de.DuoException("--> Invalid argument.")

            self.elementList[Z] = ec

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Commit(self):
        # validate fraction values
        # the fraction values should be either all weight fractions or all atomic fraction
        # mix of both is not allowed
        useWeightFraction = False;
        useAtomicFraction = False;
        for Z, ec in self.elementList.items():
            if (not useWeightFraction) and (ec.weightFraction != 0.0):
                useWeightFraction = True

            if (not useAtomicFraction) and (ec.atomicFraction != 0.0):
                useAtomicFraction = True

        if useWeightFraction and useAtomicFraction:
            raise de.DuoException("--> Mixing weight fraction and atomic fraction is not supported.")

        if (not useWeightFraction) and (not useAtomicFraction):
            raise de.DuoException("--> Wrong weight fraction or atomic fraction.")

        # convert between weight and atomic fraction
        # theory: w_i = n_i * A_i / A
        if useWeightFraction:
            for Z, ec in self.elementList.items():
                element = self.elementTable.GetElementByZ(Z)
                ec.atomicFraction = ec.weightFraction / element.A

        elif useAtomicFraction:
            for Z, ec in self.elementList.items():
                element = self.elementTable.GetElementByZ(Z)
                ec.weightFraction = ec.atomicFraction * element.A

        # normalize fraction values
        weightFractionSum = 0.0
        atomicFractionSum = 0.0

        for Z, ec in self.elementList.items():
            weightFractionSum += ec.weightFraction
            atomicFractionSum += ec.atomicFraction

        for Z, ec in self.elementList.items():
            ec.weightFraction /= weightFractionSum
            ec.atomicFraction /= atomicFractionSum


        self.isCommitted = True

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Show(self):
        message = "--> material : {0:s}".format(self.name)
        print(message)

        for Z, ec in self.elementList.items():
            message = "    Z = {0:5d}, w = {1:10f}, n = {2:10f}".format(Z, ec.weightFraction, ec.atomicFraction)
            print(message)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateMacAtE(self, energy):
        if not self.isCommitted:
            raise de.DuoException("--> Material not committed.")

        mac = 0.0

        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)
            mac += ec.weightFraction / element.A * element.CalculateTotalMicroXSAtE(energy)

        mac *= constant.Endfb.Avogadro * 1e-24

        return mac # unit: cm ^ 2 / g

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculatePEMacAtE(self, energy):
        if not self.isCommitted:
            raise de.DuoException("--> Material not committed.")

        mac = 0.0

        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)
            mac += ec.weightFraction / element.A * element.CalculatePEMicroXSAtE(energy)

        mac *= constant.Endfb.Avogadro * 1e-24

        return mac # unit: cm ^ 2 / g

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateCSMacAtE(self, energy):
        if not self.isCommitted:
            raise de.DuoException("--> Material not committed.")

        mac = 0.0

        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)
            mac += ec.weightFraction / element.A * element.CalculateCSMicroXSAtE(energy)

        mac *= constant.Endfb.Avogadro * 1e-24

        return mac # unit: cm ^ 2 / g

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateRLMacAtE(self, energy):
        if not self.isCommitted:
            raise de.DuoException("--> Material not committed.")

        mac = 0.0

        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)
            mac += ec.weightFraction / element.A * element.CalculateRLMicroXSAtE(energy)

        mac *= constant.Endfb.Avogadro * 1e-24

        return mac # unit: cm ^ 2 / g

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateZeffAtE(self, energy):
        if not self.isCommitted:
            raise de.DuoException("--> Material not committed.")

        Zeff = 0.0
        up = 0.0
        down = 0.0

        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)
            sigma = element.CalculateTotalMicroXSAtE(energy)
            temp = ec.atomicFraction * sigma
            up   += temp
            down += temp / Z
        Zeff = up / down

        return Zeff

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateElectronXSAtE(self, energy):
        if not self.isCommitted:
            raise de.DuoException("--> Material not committed.")

        up = 0.0
        down = 0.0
        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)
            sigma_e = element.CalculateElectronXSAtE(energy)
            up   += ec.weightFraction * element.Z * sigma_e / element.A
            down += ec.weightFraction * element.Z / element.A

        result = up / down
        return result

    #------------------------------------------------------------
    #------------------------------------------------------------
    def ShowInfo(self, energy):
        print("--> Material ShowInfo(): energy = {:f}".format(energy))

        totalElectrons = 0
        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)
            totalElectrons += element.Z * ec.atomicFraction

        up = 0.0
        down = 0.0
        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)
            sigma_e = element.CalculateElectronXSAtE(energy)
            up   += ec.weightFraction * element.Z * sigma_e / element.A
            down += ec.weightFraction * element.Z / element.A

            currentElectrons = element.Z * ec.atomicFraction
            print("    Z = {:d}, lambda = {:f}, sigma_e = {:f}".format(element.Z, currentElectrons / totalElectrons, sigma_e))

        result = up / down
        print("    electron xs = {:f}".format(result))

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateTotalMicroXSAtEPerAtom(self, energy):
        up = 0.0
        down = 0.0

        for Z, ec in self.elementList.items():
            element = self.elementTable.GetElementByZ(Z)

            temp = ec.weightFraction / element.A

            up += temp * element.CalculateTotalMicroXSAtE(energy)
            down += temp

        return up / down

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateNumElectronPerAtom(self):
        result = 0.0

        for Z, ec in self.elementList.items():
            result += Z * ec.atomicFraction

        return result

