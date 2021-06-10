import os
import sys
import numpy as np
import duo.zeff.nist as nist
import duo.core.element_table as element_table
import duo.core.material as material
import duo.zeff.common as common
import duo.zeff.bourque as bourque
import duo.zeff.taylor as taylor
import duo.core.duo_exception as de

#------------------------------------------------------------
# Given a known material, calculate Z, mu, CTNumber at high
# and low energies, and ZeffAve
#------------------------------------------------------------
class ZeffCalculator:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, material, method = 'bourque', isNist = False):
        self.com            = com
        self.method         = method

        self.Ehigh          = com.Ehigh
        self.Z_Ehigh        = 0.0
        self.mu_Ehigh       = 0.0
        self.CTNumber_Ehigh = 0.0

        self.Elow           = com.Elow
        self.Z_Elow         = 0.0
        self.mu_Elow        = 0.0
        self.CTNumber_Elow  = 0.0

        self.ZeffAve = 0.0
        self.material = material

        self.isChosen = False
        self.isNist = isNist

        if self.method == 'bourque':
            self.bq = bourque.Bourque(com, method = 'bspline')
        elif self.method == 'taylor':
            self.tl = taylor.Taylor(com)


    #------------------------------------------------------------
    #------------------------------------------------------------
    def Calculate(self):
        # high
        if self.method == 'direct':
            self.Z_Ehigh  = self.material.CalculateZeffAtE(self.Ehigh)
        elif self.method == 'bourque':
            try:
                self.Z_Ehigh  = self.bq.CalculateZeffAtE(self.material, self.Ehigh)
            except de.DuoException as err:
                msg = "--> Error when calculating Z_Ehigh for material {:s} using {:s} : {:s}".format(self.material.name, self.method, str(err))
                print(msg)
                self.Z_Ehigh = -1.0
        elif self.method == 'taylor':
            try:
                self.Z_Ehigh  = self.tl.CalculateZeffAtE(self.material, self.Ehigh)
            except de.DuoException as err:
                msg = "--> Error when calculating Z_Ehigh for material {:s} using {:s} : {:s}".format(self.material.name, self.method, str(err))
                print(msg)
                self.Z_Ehigh = -1.0

        self.mu_Ehigh = self.material.CalculateMacAtE(self.Ehigh) * self.material.density
        self.CTNumber_Ehigh = self.com.ConvertMuToCTNumber(self.mu_Ehigh, self.com.muWater_Ehigh, self.com.muAir_Ehigh)

        # low
        if self.method == 'direct':
            self.Z_Elow   = self.material.CalculateZeffAtE(self.Elow)
        elif self.method == 'bourque':
            try:
                self.Z_Elow  = self.bq.CalculateZeffAtE(self.material, self.Elow)
            except de.DuoException as err:
                msg = "--> Error when calculating Z_Elow for material {:s} using {:s} : {:s}".format(self.material.name, self.method, str(err))
                print(msg)
                self.Z_Elow  = -1.0
        elif self.method == 'taylor':
            try:
                self.Z_Elow  = self.tl.CalculateZeffAtE(self.material, self.Elow)
            except de.DuoException as err:
                msg = "--> Error when calculating Z_Elow for material {:s} using {:s} : {:s}".format(self.material.name, self.method, str(err))
                print(msg)
                self.Z_Elow  = -1.0

        self.mu_Elow  = self.material.CalculateMacAtE(self.Elow) * self.material.density
        self.CTNumber_Elow = self.com.ConvertMuToCTNumber(self.mu_Elow, self.com.muWater_Elow, self.com.muAir_Elow)

        # ave
        self.ZeffAve = (self.Z_Elow + self.Z_Ehigh) / 2.0

        # choose according to CT number
        if self.CTNumber_Elow  >= self.com.CTNumber_Elow and \
           self.CTNumber_Elow  <= self.com.CTNumber_Ehigh and \
           self.CTNumber_Ehigh >= self.com.CTNumber_Elow and \
           self.CTNumber_Ehigh <= self.com.CTNumber_Ehigh:
            self.isChosen = True

        # choose manually
        # remove these material according to Cai's paper
        if self.material.name.find("Tissue-Equivalent Gas, Methane Based") != -1 or\
           self.material.name.find("Tissue-Equivalent Gas, Propane Based") != -1:
            self.isChosen = False

        if self.material.name.find("Polytetrafluoroethylene, (Teflon)") != -1:
            self.isChosen = False

        # retain these material according to Cai's paper
        if self.isNist:
            if self.material.name.find("A-150 Tissue-Equivalent Plastic") != -1 or\
               self.material.name.find("Adipose Tissue (ICRU-44)") != -1 or\
               self.material.name.find("Air, Dry (near sea level)") != -1 or\
               self.material.name.find("B-100 Bone-Equivalent Plastic") != -1 or\
               self.material.name.find("Blood, Whole (ICRU-44)") != -1 or\
               self.material.name.find("Bone, Cortical (ICRU-44)") != -1 or\
               self.material.name.find("Brain, Grey/White Matter (ICRU-44)") != -1 or\
               self.material.name.find("Breast Tissue (ICRU-44)") != -1 or\
               self.material.name.find("Eye Lens (ICRU-44)") != -1 or\
               self.material.name.find("Lung Tissue (ICRU-44)") != -1 or\
               self.material.name.find("Muscle, Skeletal (ICRU-44)") != -1 or\
               self.material.name.find("Ovary (ICRU-44)") != -1 or\
               self.material.name.find("Testis (ICRU-44)") != -1 or\
               self.material.name.find("Tissue, Soft (ICRU-44)") != -1 or\
               self.material.name.find("Tissue, Soft (ICRU Four-Component)") != -1 or\
               self.material.name.find("Water, Liquid") != -1:
                self.isChosen = True
            else:
                self.isChosen = False
        else:
            self.isChosen = True


