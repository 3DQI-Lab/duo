import os
import sys
import numpy as np
import duo.zeff.nist as nist
import duo.core.element_table as element_table
import duo.core.material as material
import duo.core.mixture as mixture
import duo.core.dicom_decoder as dicom_decoder
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import duo.zeff.zeff_calculator as zeff_calculator
import duo.zeff.blend as blend
import duo.zeff.enhance as enhance
import duo.core.surface_fitting as sf
import scipy.interpolate as sip
from matplotlib import cm
import scipy
import scipy.ndimage
import pydicom
import duo.core.timer as timer
import matplotlib.colors
import cv2 as cv
import duo.zeff.common as common
import copy

#------------------------------------------------------------
#------------------------------------------------------------
class ControlPointContainer:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        self.xList = []
        self.yList = []
        self.zList = []
        self.textList = []

#------------------------------------------------------------
#------------------------------------------------------------
class ControlPointManager:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method, isToDumpIntermediateData = True):
        self.com = com
        self.method = method
        self.isToDumpIntermediateData = isToDumpIntermediateData
        self.elementTable = com.elementTable

        self.nist = None # 16 nist data
        self.extraMaterialList = [] # 8 mixture data
        self.fullMaterialList = [] # 16 nist data + 8 mixture data

        self.zeffCalculatorList = []
        self.cpListTotal = ControlPointContainer()
        self.cpListNist = ControlPointContainer()
        self.cpListCai = ControlPointContainer()

        self.rbf = None
        self.Initialize()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        # instantiate nist object
        self.nist = nist.Nist(self.com)

        # instantiate extra materials
        self.AddExtraMaterial2()

        # concatenate material lists
        self.fullMaterialList = self.nist.materialList + self.extraMaterialList

        # control points consist of 16 nist materials + a bunch of extra materials
        self.CalculateZeffForControlPoint()

        self.InitializeControlPointContainer()

        self.InitializeRBF()

        if self.isToDumpIntermediateData:
            self.PlotDiscretePoint()

            self.PlotFullSurface()

            # self.PlotSelectedData()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetMaterialFromNist(self, name):
        material = None
        for mat in self.nist.materialList:
            if mat.name.find(name) != -1:
                material = mat
                break
        return material

    #------------------------------------------------------------
    #------------------------------------------------------------
    def InitializeControlPointContainer(self):
        # # manual adjustment
        # for item in self.zeffCalculatorList:
            # if item.material.name == 'Air, Dry (near sea level)':
                # item.ZeffAve = 30.0

        for item in self.zeffCalculatorList:
            if item.isChosen:
                self.cpListTotal.xList.append(item.CTNumber_Ehigh)
                self.cpListTotal.yList.append(item.CTNumber_Elow)
                self.cpListTotal.zList.append(item.ZeffAve)
                self.cpListTotal.textList.append(item.material.name)

        for item in self.zeffCalculatorList:
            if item.isChosen and item.isNist:
                self.cpListNist.xList.append(item.CTNumber_Ehigh)
                self.cpListNist.yList.append(item.CTNumber_Elow)
                self.cpListNist.zList.append(item.ZeffAve)
                self.cpListNist.textList.append(item.material.name)
        print("    number of NIST data = {0:d}".format(len(self.cpListNist.xList)))

        for item in self.zeffCalculatorList:
            if item.isChosen and not item.isNist:
                self.cpListCai.xList.append(item.CTNumber_Ehigh)
                self.cpListCai.yList.append(item.CTNumber_Elow)
                self.cpListCai.zList.append(item.ZeffAve)
                self.cpListCai.textList.append(item.material.name)
        print("    number of non-NIST data = {0:d}".format(len(self.cpListCai.xList)))

    #------------------------------------------------------------
    #------------------------------------------------------------
    def InitializeRBF(self):
        # build model

        # method 1: surface fit. not used because the generated polynomial surface has extremely wide range
        # self.model = sf.PolyFit2D(np.array(self.collocationPointList.xList),
                                  # np.array(self.collocationPointList.yList),
                                  # np.array(self.collocationPointList.zList),
                                  # order = 3)

        # x = np.linspace(-1000, 3095, 4096)
        # y = np.linspace(-1000, 3095, 4096)
        # xx, yy = np.meshgrid(x, y)
        # zz = sf.PolyVal2D(xx, yy, self.model)

        # method 2: interpolation
        # according to matlab doc: https://www.mathworks.com/help/curvefit/interpolation-methods.html
        # "Biharmonic for surfaces which is a radial basis function interpolant"
        # "Biharmonic spline interpolation"
        #
        # scipy doc
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.Rbf.html
        #
        # according to this projecthttps://github.com/mattfoster/matlab-interpolation-toolkit/blob/master/rbf/rbf.m
        # and matlab source code in griddata.m, biharmonic (aka v4) is equivalent to thin plate
        self.rbf = sip.Rbf(self.cpListTotal.xList,
                           self.cpListTotal.yList,
                           self.cpListTotal.zList,
                           function = 'thin_plate')

        # calculate rbf error
        sum = 0.0
        for item in self.zeffCalculatorList:
            ZRbf = self.rbf(item.CTNumber_Ehigh, item.CTNumber_Elow)
            diff = np.absolute(ZRbf - item.ZeffAve) / item.ZeffAve
            sum += diff
        sum /= len(self.zeffCalculatorList)
        print("    rbf error: average diff = {0:.6e}%".format(sum * 100.0))

    #------------------------------------------------------------
    #------------------------------------------------------------
    def PlotDiscretePoint(self):
        print("--> Plot discrete point")

        fig = plt.figure(figsize=(16, 10))
        ax = fig.add_subplot(111, projection = '3d')
        counter = 1
        for i in range(len(self.cpListNist.xList)):
            ax.scatter(self.cpListNist.xList[i],
                       self.cpListNist.yList[i],
                       self.cpListNist.zList[i],
                       c = '#ff0000',
                       depthshade = False,
                       label = str(counter) + " : "+ self.cpListNist.textList[i])

            ax.text(self.cpListNist.xList[i],
                    self.cpListNist.yList[i],
                    self.cpListNist.zList[i] + 1,
                    str(counter), fontsize = 10)

            counter += 1

        for i in range(len(self.cpListCai.xList)):
            ax.scatter(self.cpListCai.xList[i],
                       self.cpListCai.yList[i],
                       self.cpListCai.zList[i],
                       c = '#0000ff',
                       depthshade = False,
                       label = str(counter) + " : "+ self.cpListCai.textList[i])

            ax.text(self.cpListCai.xList[i],
                    self.cpListCai.yList[i],
                    self.cpListCai.zList[i] + 1,
                    str(counter), fontsize = 10)

            counter += 1

        ax.set_xlim([-1000, 3095])
        ax.set_xlabel("CT number at 69.28 keV")

        ax.set_ylim([-1000, 3095])
        ax.set_ylabel("CT number at 51.93 keV")

        # upperZ = np.amax([np.amax(self.cpListNist.zList), np.amax(self.cpListCai.zList)]) + 1
        ax.set_zlim(bottom = 0.0)
        ax.set_zlabel("$Z_{eff}$")

        # add vertical lines
        for i in range(len(self.cpListTotal.xList)):
            ax.plot([self.cpListTotal.xList[i], self.cpListTotal.xList[i]],
                    [self.cpListTotal.yList[i], self.cpListTotal.yList[i]],
                    [0, self.cpListTotal.zList[i]],
                    color = '#000000', linestyle = 'dotted', linewidth = 1)

        # no simple way to customize 3d plot, so according to stackoverflow:
        # set pane color (background color)
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

        ax.tick_params(bottom = True, top = True, left = True, right = True)

        plt.legend(loc='center left', bbox_to_anchor=(1.07, 0.5), edgecolor = '#000000', shadow = True)
        fig.subplots_adjust(right = 0.7)
        plt.savefig("discrete_point_" + self.method +".pdf")


    #------------------------------------------------------------
    #------------------------------------------------------------
    def PlotFullSurface(self):
        """
        Known issue: Matplotlib 3.4.2 still has a bug where 3D plot zorder is oftentimes ignored.
        """
        print("--> Plot full surface")

        ZMax = 0.0
        for item in self.zeffCalculatorList:
            if item.ZeffAve > ZMax:
                ZMax = item.ZeffAve

        ZMin = 100.0
        for item in self.zeffCalculatorList:
            if item.ZeffAve < ZMin:
                ZMin = item.ZeffAve

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection = '3d')

        #------------------------------------------------------------
        # plot surface
        #------------------------------------------------------------
        x = np.linspace(-1000, 3095, 32)
        y = np.linspace(-1000, 3095, 32)
        xx, yy = np.meshgrid(x, y)
        zz = self.rbf(xx, yy)

        # clamp
        zeffMax = 11
        zeffMin = 1
        zz = np.where(zz > zeffMax, zeffMax, zz)
        zz = np.where(zz < zeffMin, zeffMin, zz)

        ax.plot_surface(xx,
                        yy,
                        zz,
                        cmap = 'jet',
                        antialiased = True,
                        vmin = 1,
                        vmax = 10,
                        linewidth = 0.1,
                        edgecolor='#00000088',
                        alpha = 0.6,
                        zorder = 1)

        #------------------------------------------------------------
        # plot points
        #------------------------------------------------------------
        counter = 1
        for i in range(len(self.cpListNist.xList)):
            ax.scatter(self.cpListNist.xList[i],
                       self.cpListNist.yList[i],
                       self.cpListNist.zList[i],
                       c = '#ff0000',
                       depthshade = False,
                       label = str(counter) + " : "+ self.cpListNist.textList[i],
                       zorder = 5)

            ax.text(self.cpListNist.xList[i],
                    self.cpListNist.yList[i],
                    self.cpListNist.zList[i] + 1,
                    str(counter), fontsize = 5)

            counter += 1

        for i in range(len(self.cpListCai.xList)):
            ax.scatter(self.cpListCai.xList[i],
                       self.cpListCai.yList[i],
                       self.cpListCai.zList[i],
                       c = '#0000ff',
                       depthshade = False,
                       label = str(counter) + " : "+ self.cpListCai.textList[i],
                       zorder = 5)

            ax.text(self.cpListCai.xList[i],
                    self.cpListCai.yList[i],
                    self.cpListCai.zList[i] + 1,
                    str(counter), fontsize = 5)

            counter += 1

        ax.set_xlim([-1000, 3095])
        ax.set_xlabel("CT number at $E_{high}$ keV")

        ax.set_ylim([-1000, 3095])
        ax.set_ylabel("CT number at $E_{low}$ keV")

        ax.set_zlim(zmin = 1, zmax = 11)
        ax.set_zlabel("Effective atomic number")
        ax.set_zticks(np.arange(1, 13, 2)) # 1, 3, 5 ... 15

        # add vertical lines
        for i in range(len(self.cpListTotal.xList)):
            ax.plot([self.cpListTotal.xList[i], self.cpListTotal.xList[i]],
                    [self.cpListTotal.yList[i], self.cpListTotal.yList[i]],
                    [1, self.cpListTotal.zList[i]],
                    color = '#000000', linestyle = 'dotted', linewidth = 1)

        # no simple way to customize 3d plot, so according to stackoverflow:
        # set pane color (background color)
        ax.w_xaxis.set_pane_color((0.9, 0.9, 0.9, 1.0))
        ax.w_yaxis.set_pane_color((0.9, 0.9, 0.9, 1.0))
        ax.w_zaxis.set_pane_color((0.9, 0.9, 0.9, 1.0))

        ax.tick_params(bottom = True, top = True, left = True, right = True)
        plt.legend(loc='center left', bbox_to_anchor=(1.07, 0.5), edgecolor = '#000000', shadow = True)
        fig.subplots_adjust(right = 0.7)

        plt.savefig("full_surface_" + self.method +".pdf", bbox_inches='tight')

    # #------------------------------------------------------------
    # #------------------------------------------------------------
    # def PlotSelectedData(self):
        # print("--> Plot selected data")


        # fig = plt.figure(figsize = (8, 8))
        # ax = fig.add_subplot(111, projection = '3d')

        # #------------------------------------------------------------
        # # plot surface
        # #------------------------------------------------------------
        # ax.plot_trisurf(self.predictedData.xList,
                        # self.predictedData.yList,
                        # self.predictedData.zList,
                        # cmap = cm.viridis,
                        # antialiased = True,
                        # linewidth = 0,
                        # alpha = 0.8,
                        # zorder = 1)

        # ax.set_xlim(-1000, 3095)
        # ax.set_xlabel("CT number at 69.28 keV")

        # ax.set_ylim(-1000, 3095)
        # ax.set_ylabel("CT number at 51.93 keV")

        # ax.set_zlabel("Average effective atomic number")
        # ax.locator_params(axis = 'z', nbins = 12) # nbins is the max number of bins

        # #------------------------------------------------------------
        # # plot scatter
        # #------------------------------------------------------------
        # lineNist, = ax.plot(self.ControlPointListNist.xList,
                            # self.ControlPointListNist.yList,
                            # self.ControlPointListNist.zList,
                            # markerfacecolor = '#ff0000',
                            # marker = 'o',
                            # markeredgewidth = 0,
                            # linestyle = 'none',
                            # zorder = 5)
        # lineNoneNist, = ax.plot(self.ControlPointListNoneNist.xList,
                                # self.ControlPointListNoneNist.yList,
                                # self.ControlPointListNoneNist.zList,
                                # markerfacecolor = '#0000ff',
                                # marker = 'o',
                                # markeredgewidth = 0,
                                # linestyle = 'none',
                                # zorder = 5)

        # # add vertical lines
        # for i in range(len(self.cpListTotal.xList)):
            # ax.plot([self.cpListTotal.xList[i], self.cpListTotal.xList[i]],
                    # [self.cpListTotal.yList[i], self.cpListTotal.yList[i]],
                    # [0.0, self.cpListTotal.zList[i]],
                    # color = '#000000', linestyle = 'dotted', linewidth = 1,
                    # zorder = 5)

        # lgd = plt.legend([lineNist, lineNoneNist], ['16 NIST materials', '8 mixture materials'], shadow = True)
        # lgd.get_frame().set_edgecolor('#000000')

        # plt.savefig("selected_data.png", bbox_inches='tight', dpi=200)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def CalculateZeffForControlPoint(self):
        print("--> Calculate Zeff for Control points")

        #------------------------------------------------------------
        # calculate 48 nist materials
        # choose 16 from them
        #------------------------------------------------------------
        for mat in self.nist.materialList:
            dm = zeff_calculator.ZeffCalculator(self.com, mat, self.method, True)
            dm.Calculate()
            if dm.isChosen:
                self.zeffCalculatorList.append(dm)

        #------------------------------------------------------------
        # calculate extra material
        #------------------------------------------------------------
        for mat in self.extraMaterialList:
            dm = zeff_calculator.ZeffCalculator(self.com, mat, self.method)
            dm.Calculate()
            self.zeffCalculatorList.append(dm)

        # sort data list
        self.zeffCalculatorList = sorted(self.zeffCalculatorList, key = lambda dm : dm.ZeffAve)

        if self.isToDumpIntermediateData:
            self.WriteControlPointToFile()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def WriteControlPointToFile(self):
        with open("control_points_" + self.method + ".txt", "w") as outfile:
            msg = "Mean energy for high kVp: {:f}\n".format(self.com.Ehigh)
            outfile.write(msg)

            msg = "Mean energy for low kVp: {:f}\n\n\n\n".format(self.com.Elow)
            outfile.write(msg)

            msg = "{:50s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}{:>20s}\n".format("material",
                                                                                            "ZeffAve",
                                                                                            "Z_Elow",
                                                                                            "Z_Ehigh",
                                                                                            "Z variation [%]",
                                                                                            "density",
                                                                                            "CTNumber_Elow",
                                                                                            "CTNumber_Ehigh",
                                                                                            "Chosen")
            outfile.write(msg)

            msg = "\n\n\n--> NIST material\n"
            outfile.write(msg)

            for item in self.zeffCalculatorList:
                if item.isNist:
                    isChosenStr = ""
                    if item.isChosen:
                        isChosenStr = "yes"

                    zVariation = np.fabs(item.Z_Ehigh - item.Z_Elow) / (item.Z_Ehigh + item.Z_Elow) * 2.0 * 100.0
                    msg = "{:50s}{:20.6f}{:20.6f}{:20.6f}{:20.6f}{:20.6f}{:20.6f}{:20.6f}{:>20s}\n".format(item.material.name,
                                                                                                           item.ZeffAve,
                                                                                                           item.Z_Elow,
                                                                                                           item.Z_Ehigh,
                                                                                                           zVariation,
                                                                                                           item.material.density,
                                                                                                           item.CTNumber_Elow,
                                                                                                           item.CTNumber_Ehigh,
                                                                                                           isChosenStr)
                    outfile.write(msg)

            msg = "\n\n\n--> non-NIST material\n"
            outfile.write(msg)

            for item in self.zeffCalculatorList:
                if not item.isNist:
                    isChosenStr = ""
                    if item.isChosen:
                        isChosenStr = "yes"

                    zVariation = np.fabs(item.Z_Ehigh - item.Z_Elow) / (item.Z_Ehigh + item.Z_Elow) * 2.0 * 100.0
                    msg = "{:50s}{:20.6f}{:20.6f}{:20.6f}{:20.6f}{:20.6f}{:20.6f}{:20.6f}{:>20s}\n".format(item.material.name,
                                                                                                           item.ZeffAve,
                                                                                                           item.Z_Elow,
                                                                                                           item.Z_Ehigh,
                                                                                                           zVariation,
                                                                                                           item.material.density,
                                                                                                           item.CTNumber_Elow,
                                                                                                           item.CTNumber_Ehigh,
                                                                                                           isChosenStr)
                    outfile.write(msg)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddManuallyTunedExtraMaterial(self):
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # add adipose tissue variation
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        scaleFactorList = [0.9, 1.5, 2.1, 2.7, 3.3, 3.9]
        for idx in range(len(scaleFactorList)):
            tempMaterial = material.Material("Adipose Tissue (ICRU-44) (Var {:d})".format(idx), self.elementTable)
            tempMaterial.AddElement(1 , weightFraction = 0.114000)
            tempMaterial.AddElement(6 , weightFraction = 0.598000)
            tempMaterial.AddElement(7 , weightFraction = 0.007000)
            tempMaterial.AddElement(8 , weightFraction = 0.278000)
            tempMaterial.AddElement(11, weightFraction = 0.001000)
            tempMaterial.AddElement(16, weightFraction = 0.001000)
            tempMaterial.AddElement(17, weightFraction = 0.001000)
            tempMaterial.Commit()
            tempMaterial.density = 0.95 * scaleFactorList[idx]
            self.extraMaterialList.append(tempMaterial)

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # add breast tissue variation
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        scaleFactorList = [0.9, 1.5, 2.1, 2.7, 3.3, 3.9]
        for idx in range(len(scaleFactorList)):
            tempMaterial = material.Material("Breast Tissue (ICRU-44) (Var {:d})".format(idx), self.elementTable)
            tempMaterial.AddElement(1 , weightFraction = 0.106000)
            tempMaterial.AddElement(6 , weightFraction = 0.332000)
            tempMaterial.AddElement(7 , weightFraction = 0.030000)
            tempMaterial.AddElement(8 , weightFraction = 0.527000)
            tempMaterial.AddElement(11, weightFraction = 0.001000)
            tempMaterial.AddElement(15, weightFraction = 0.001000)
            tempMaterial.AddElement(16, weightFraction = 0.002000)
            tempMaterial.AddElement(17, weightFraction = 0.001000)
            tempMaterial.Commit()
            tempMaterial.density = 1.02 * scaleFactorList[idx]
            self.extraMaterialList.append(tempMaterial)

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # add soft tissue variation
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        scaleFactorList = [0.9, 1.5, 2.1, 2.7, 3.3]
        for idx in range(len(scaleFactorList)):
            tempMaterial = material.Material("Tissue, Soft (ICRU-44) (Var {:d})".format(idx), self.elementTable)
            tempMaterial.AddElement(1 , weightFraction = 0.102000)
            tempMaterial.AddElement(6 , weightFraction = 0.143000)
            tempMaterial.AddElement(7 , weightFraction = 0.034000)
            tempMaterial.AddElement(8 , weightFraction = 0.708000)
            tempMaterial.AddElement(11, weightFraction = 0.002000)
            tempMaterial.AddElement(15, weightFraction = 0.003000)
            tempMaterial.AddElement(16, weightFraction = 0.003000)
            tempMaterial.AddElement(17, weightFraction = 0.002000)
            tempMaterial.AddElement(19, weightFraction = 0.003000)
            tempMaterial.Commit()
            tempMaterial.density = 1.06 * scaleFactorList[idx]
            self.extraMaterialList.append(tempMaterial)

        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # add air variation
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #scaleFactorList = [5.0, 10.0, 20.0]
        #for idx in range(len(scaleFactorList)):
        #    tempMaterial = material.Material("Air, Dry (near sea level)  (Var {:d})".format(idx), self.elementTable)
        #    tempMaterial.AddElement(6 , weightFraction = 0.000124)
        #    tempMaterial.AddElement(7 , weightFraction = 0.755268)
        #    tempMaterial.AddElement(8 , weightFraction = 0.231781)
        #    tempMaterial.AddElement(18, weightFraction = 0.012827)
        #    tempMaterial.Commit()
        #    tempMaterial.density = 1.205e-03 * scaleFactorList[idx]
        #    self.extraMaterialList.append(tempMaterial)




        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        # add bone variation
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        scaleFactorList = [0.9, 1.3, 1.7]
        for idx in range(len(scaleFactorList)):
            tempMaterial = material.Material("B-100 Bone-Equivalent Plastic (Var {:d})".format(idx), self.elementTable)
            tempMaterial.AddElement(1 , weightFraction = 0.065473)
            tempMaterial.AddElement(6 , weightFraction = 0.536942)
            tempMaterial.AddElement(7 , weightFraction = 0.021500)
            tempMaterial.AddElement(8 , weightFraction = 0.032084)
            tempMaterial.AddElement(9 , weightFraction = 0.167415)
            tempMaterial.AddElement(20, weightFraction = 0.176585)
            tempMaterial.Commit()
            tempMaterial.density = 1.45 * scaleFactorList[idx]
            self.extraMaterialList.append(tempMaterial)

        scaleFactorList = [0.9, 1.15]
        for idx in range(len(scaleFactorList)):
            tempMaterial = material.Material("Bone, Cortical (ICRU-44) (Var {:d})".format(idx), self.elementTable)
            tempMaterial.AddElement(1   , weightFraction = 0.034000)
            tempMaterial.AddElement(6   , weightFraction = 0.155000)
            tempMaterial.AddElement(7   , weightFraction = 0.042000)
            tempMaterial.AddElement(8   , weightFraction = 0.435000)
            tempMaterial.AddElement(11  , weightFraction = 0.001000)
            tempMaterial.AddElement(12  , weightFraction = 0.002000)
            tempMaterial.AddElement(15  , weightFraction = 0.103000)
            tempMaterial.AddElement(16  , weightFraction = 0.003000)
            tempMaterial.AddElement(20  , weightFraction = 0.225000)
            tempMaterial.Commit()
            tempMaterial.density = 1.92 * scaleFactorList[idx]
            self.extraMaterialList.append(tempMaterial)



#------------------------------------------------------------
#------------------------------------------------------------
class ControlPointManagerColonEC(ControlPointManager):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method, isToDumpIntermediateData = True):
        super(ControlPointManagerColonEC, self).__init__(com, method, isToDumpIntermediateData)

    #------------------------------------------------------------
    # This function is not actually used
    #------------------------------------------------------------
    def AddExtraMaterial(self):
        # dry air
        dryAir = self.nist.SearchMaterialFromNist("Air, Dry (near sea level)")

        # water
        water = self.nist.SearchMaterialFromNist("Water, Liquid")

        # fat
        fat = self.nist.SearchMaterialFromNist("Adipose Tissue (ICRU-44)")

        # soft tissue
        softTissue = self.nist.SearchMaterialFromNist("Tissue, Soft (ICRU-44)")

        # iodine + water
        iodine = material.Material("iodine", self.elementTable)
        iodine.AddElement(53, weightFraction = 1.0)
        iodine.Commit()
        iodine.density = 4.93

        iodineWater = mixture.Mixture("iodine-water", self.elementTable)
        iodineWater.AddMaterial(iodine, materialWeightFraction = 20.0 / 1000.0) # 20 mg iodine
        iodineWater.AddMaterial(water, materialWeightFraction = 1.0) # 1 cc water, i.e. 1 g
        iodineWater.Commit()
        iodineWater.CalculateMixtureDensity()
        self.extraMaterialList.append(iodineWater)

        # build non-nist material
        # assume the fraction given in Cai's draft is the volume fraction
        # instead of weight (mass) fraction or molar fraction
        # 9% dry air + 91% (iodine + water)
        mix1 = mixture.Mixture("air(9%)-iodine-water", self.elementTable)
        ratio = 0.09 * dryAir.density / iodineWater.density
        mix1.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix1.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix1.Commit()
        mix1.CalculateMixtureDensity()
        self.extraMaterialList.append(mix1)

        # 72% dry air + 28% fat
        mix2 = mixture.Mixture("air(72%)-fat", self.elementTable)
        ratio = 0.72 * dryAir.density / fat.density
        mix2.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix2.AddMaterial(fat, materialWeightFraction = 1.0 - ratio)
        mix2.Commit()
        mix2.CalculateMixtureDensity()
        self.extraMaterialList.append(mix2)

        # 72% dry air + 28% soft tissue
        mix3 = mixture.Mixture("air(72%)-softTissue", self.elementTable)
        ratio = 0.72 * dryAir.density / softTissue.density
        mix3.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix3.AddMaterial(softTissue, materialWeightFraction = 1.0 - ratio)
        mix3.Commit()
        mix3.CalculateMixtureDensity()
        self.extraMaterialList.append(mix3)

        # 23% dry air + 77% (iodine + water)
        mix4 = mixture.Mixture("air(23%)-iodine-water", self.elementTable)
        ratio = 0.23 * dryAir.density / iodineWater.density
        mix4.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix4.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix4.Commit()
        mix4.CalculateMixtureDensity()
        self.extraMaterialList.append(mix4)

        # 3% dry air + 97% (iodine + water)
        mix5 = mixture.Mixture("air(3%)-iodine-water", self.elementTable)
        ratio = 0.03 * dryAir.density / iodineWater.density
        mix5.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix5.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix5.Commit()
        mix5.CalculateMixtureDensity()
        self.extraMaterialList.append(mix5)

        # 30% dry air + 70% (iodine + water)
        mix6 = mixture.Mixture("air(30%)-iodine-water", self.elementTable)
        ratio = 0.3 * dryAir.density / iodineWater.density
        mix6.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix6.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix6.Commit()
        mix6.CalculateMixtureDensity()
        self.extraMaterialList.append(mix6)

        # 72% soft tissue + 28% (iodine + water)
        mix7 = mixture.Mixture("softTissue(72%)-iodine-water", self.elementTable)
        ratio = 0.72 * softTissue.density / iodineWater.density
        mix7.AddMaterial(softTissue, materialWeightFraction = ratio)
        mix7.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix7.Commit()
        mix7.CalculateMixtureDensity()
        self.extraMaterialList.append(mix7)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddExtraMaterial2(self):
        # dry air
        dryAir = self.nist.SearchMaterialFromNist("Air, Dry (near sea level)")

        # water
        water = self.nist.SearchMaterialFromNist("Water, Liquid")

        # fat
        fat = self.nist.SearchMaterialFromNist("Adipose Tissue (ICRU-44)")

        # soft tissue
        softTissue = self.nist.SearchMaterialFromNist("Tissue, Soft (ICRU-44)")

        # iodine + water
        iodine = material.Material("iodine", self.elementTable)
        iodine.AddElement(53, weightFraction = 1.0)
        iodine.Commit()
        iodine.density = 4.93

        iodineWater = mixture.Mixture("iodine-water (20 mg/ml)", self.elementTable)
        iodineWater.AddMaterial(iodine, materialWeightFraction = 20.0 / 1000.0) # 20 mg iodine
        iodineWater.AddMaterial(water, materialWeightFraction = 1.0) # 1 cc water, i.e. 1 g
        iodineWater.Commit()
        iodineWater.CalculateMixtureDensity()
        self.extraMaterialList.append(iodineWater)

        # iodine + water
        iodineWater2 = mixture.Mixture("iodine-water (50 mg/ml)", self.elementTable)
        iodineWater2.AddMaterial(iodine, materialWeightFraction = 50.0 / 1000.0) # 50 mg iodine
        iodineWater2.AddMaterial(water, materialWeightFraction = 1.0) # 1 cc water, i.e. 1 g
        iodineWater2.Commit()
        iodineWater2.CalculateMixtureDensity()
        self.extraMaterialList.append(iodineWater2)

        # build non-nist material
        # assume the fraction given in Cai's draft is the volume fraction
        # instead of weight (mass) fraction or molar fraction

        # 10% dry air + 90% (iodine + water)
        mix = mixture.Mixture("air(10%)-iodine-water", self.elementTable)
        ratio = 0.1 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        # 20% dry air + 80% (iodine + water)
        mix = mixture.Mixture("air(20%)-iodine-water", self.elementTable)
        ratio = 0.2 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        # 30% dry air + 70% (iodine + water)
        mix = mixture.Mixture("air(30%)-iodine-water", self.elementTable)
        ratio = 0.3 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        # 40% dry air + 60% (iodine + water)
        mix = mixture.Mixture("air(40%)-iodine-water", self.elementTable)
        ratio = 0.4 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        # 50% dry air + 50% (iodine + water)
        mix = mixture.Mixture("air(50%)-iodine-water", self.elementTable)
        ratio = 0.5 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        # 60% dry air + 40% (iodine + water)
        mix = mixture.Mixture("air(60%)-iodine-water", self.elementTable)
        ratio = 0.6 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        # 70% dry air + 30% (iodine + water)
        mix = mixture.Mixture("air(70%)-iodine-water", self.elementTable)
        ratio = 0.7 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        # 80% dry air + 20% (iodine + water)
        mix = mixture.Mixture("air(80%)-iodine-water", self.elementTable)
        ratio = 0.8 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        # 90% dry air + 10% (iodine + water)
        mix = mixture.Mixture("air(90%)-iodine-water", self.elementTable)
        ratio = 0.9 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)



#------------------------------------------------------------
#------------------------------------------------------------
class ControlPointManagerKidneyStone(ControlPointManager):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method, isToDumpIntermediateData = True):
        super(ControlPointManagerKidneyStone, self).__init__(com, method, isToDumpIntermediateData)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddExtraMaterial2(self):
        # uric acid
        # https://byjus.com/uric-acid-formula/
        # C5 H4 N4 O3
        uricAcid = material.Material("uric acid", self.elementTable)
        uricAcid.AddElement(6, atomicFraction = 5)
        uricAcid.AddElement(1, atomicFraction = 4)
        uricAcid.AddElement(7, atomicFraction = 4)
        uricAcid.AddElement(8, atomicFraction = 3)
        uricAcid.Commit()
        uricAcid.density = 1.87
        self.extraMaterialList.append(uricAcid)

        # calcium oxalate monohydrate
        # https://www.americanelements.com/calcium-oxalate-monohydrate-5794-28-5
        # Ca C2 H2 O5
        calciumOxMono = material.Material("calcium oxalate monohydrate", self.elementTable)
        calciumOxMono.AddElement(20, atomicFraction = 1)
        calciumOxMono.AddElement(6 , atomicFraction = 2)
        calciumOxMono.AddElement(1 , atomicFraction = 2)
        calciumOxMono.AddElement(8 , atomicFraction = 5)
        calciumOxMono.Commit()
        calciumOxMono.density = 2.2
        self.extraMaterialList.append(calciumOxMono)

        # calcium oxalate dihydrate
        # https://www.americanelements.com/calcium-oxalate-monohydrate-5794-28-5
        # Ca C2 H4 O6
        calciumOxDi = material.Material("calcium oxalate dihydrate", self.elementTable)
        calciumOxDi.AddElement(20, atomicFraction = 1)
        calciumOxDi.AddElement(6 , atomicFraction = 2)
        calciumOxDi.AddElement(1 , atomicFraction = 4)
        calciumOxDi.AddElement(8 , atomicFraction = 6)
        calciumOxDi.Commit()
        calciumOxDi.density = 2.2
        self.extraMaterialList.append(calciumOxDi)

        # hydroxyapatite
        # http://cameo.mfa.org/wiki/Calcium_hydroxyapatite
        # https://pubmed.ncbi.nlm.nih.gov/23263603/
        # Ca10 H2 O26 P6
        hydroxyapatite = material.Material("hydroxyapatite", self.elementTable)
        hydroxyapatite.AddElement(20 , atomicFraction = 10)
        hydroxyapatite.AddElement(1  , atomicFraction = 2)
        hydroxyapatite.AddElement(8  , atomicFraction = 26)
        hydroxyapatite.AddElement(15 , atomicFraction = 6)
        hydroxyapatite.Commit()
        hydroxyapatite.density = 3.18
        self.extraMaterialList.append(hydroxyapatite)

        # cystine
        # https://www.sigmaaldrich.com/catalog/product/mm/102836?lang=en&region=US
        # C6 H12 N2 O4 S2
        cystine = material.Material("cystine", self.elementTable)
        cystine.AddElement(6  , atomicFraction = 6 )
        cystine.AddElement(1  , atomicFraction = 12)
        cystine.AddElement(7  , atomicFraction = 2 )
        cystine.AddElement(8  , atomicFraction = 4 )
        cystine.AddElement(16 , atomicFraction = 2 )
        cystine.Commit()
        cystine.density = 1.677
        self.extraMaterialList.append(cystine)

        # struvite
        # https://pubmed.ncbi.nlm.nih.gov/23263603/
        # Ca Mg N H16 P O10
        struvite = material.Material("struvite", self.elementTable)
        struvite.AddElement(20 , atomicFraction = 1 )
        struvite.AddElement(12 , atomicFraction = 1 )
        struvite.AddElement(7  , atomicFraction = 1 )
        struvite.AddElement(1  , atomicFraction = 16)
        struvite.AddElement(15 , atomicFraction = 1 )
        struvite.AddElement(8  , atomicFraction = 10)
        struvite.Commit()
        struvite.density = 1.711
        self.extraMaterialList.append(struvite)


#------------------------------------------------------------
# same density, different Z
#------------------------------------------------------------
class ControlPointManagerSpecial1(ControlPointManager):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method, isToDumpIntermediateData = True):
        super(ControlPointManagerSpecial1, self).__init__(com, method, isToDumpIntermediateData)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddExtraMaterial2(self):
        ZList = np.arange(1, 36 + 1)

        for Z in ZList:
            tempMaterial = material.Material("material {:d}".format(Z), self.elementTable)
            tempMaterial.AddElement(Z, atomicFraction = 1)
            tempMaterial.Commit()
            tempMaterial.density = 0.2
            self.extraMaterialList.append(tempMaterial)

#------------------------------------------------------------
# same Z, different density
#------------------------------------------------------------
class ControlPointManagerSpecial2(ControlPointManager):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method, isToDumpIntermediateData = True):
        super(ControlPointManagerSpecial2, self).__init__(com, method, isToDumpIntermediateData)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddExtraMaterial2(self):
        Z = 26

        densityList = np.linspace(0.1, 1, 20)

        for density in densityList:
            tempMaterial = material.Material("material {:f}".format(density), self.elementTable)
            tempMaterial.AddElement(Z, atomicFraction = 1)
            tempMaterial.Commit()
            tempMaterial.density = density
            self.extraMaterialList.append(tempMaterial)


#------------------------------------------------------------
#------------------------------------------------------------
class ControlPointManagerColonECExtra(ControlPointManager):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method, isToDumpIntermediateData = True):
        super(ControlPointManagerColonECExtra, self).__init__(com, method, isToDumpIntermediateData)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddExtraMaterial2(self):
        self.AddManuallyTunedExtraMaterial()

        # dry air
        dryAir = self.nist.SearchMaterialFromNist("Air, Dry (near sea level)")

        # water
        water = self.nist.SearchMaterialFromNist("Water, Liquid")

        # fat
        fat = self.nist.SearchMaterialFromNist("Adipose Tissue (ICRU-44)")

        # soft tissue
        softTissue = self.nist.SearchMaterialFromNist("Tissue, Soft (ICRU-44)")

        # iodine + water
        iodine = material.Material("iodine", self.elementTable)
        iodine.AddElement(53, weightFraction = 1.0)
        iodine.Commit()
        iodine.density = 4.93

        iodineWater = mixture.Mixture("iodine-water (20 mg/ml)", self.elementTable)
        iodineWater.AddMaterial(iodine, materialWeightFraction = 20.0 / 1000.0) # 20 mg iodine
        iodineWater.AddMaterial(water, materialWeightFraction = 1.0) # 1 cc water, i.e. 1 g
        iodineWater.Commit()
        iodineWater.CalculateMixtureDensity()
        self.extraMaterialList.append(iodineWater)

        # iodine + water
        iodineWater2 = mixture.Mixture("iodine-water (50 mg/ml)", self.elementTable)
        iodineWater2.AddMaterial(iodine, materialWeightFraction = 50.0 / 1000.0) # 50 mg iodine
        iodineWater2.AddMaterial(water, materialWeightFraction = 1.0) # 1 cc water, i.e. 1 g
        iodineWater2.Commit()
        iodineWater2.CalculateMixtureDensity()
        self.extraMaterialList.append(iodineWater2)

        # build non-nist material
        # assume the fraction given in Cai's draft is the volume fraction
        # instead of weight (mass) fraction or molar fraction

        # 10% dry air + 90% (iodine + water)
        mix = mixture.Mixture("air(10%)-iodine-water", self.elementTable)
        ratio = 0.1 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        ## 20% dry air + 80% (iodine + water)
        #mix = mixture.Mixture("air(20%)-iodine-water", self.elementTable)
        #ratio = 0.2 * dryAir.density / iodineWater.density
        #mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        #mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        #mix.Commit()
        #mix.CalculateMixtureDensity()
        #self.extraMaterialList.append(mix)

        # 30% dry air + 70% (iodine + water)
        mix = mixture.Mixture("air(30%)-iodine-water", self.elementTable)
        ratio = 0.3 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        ## 40% dry air + 60% (iodine + water)
        #mix = mixture.Mixture("air(40%)-iodine-water", self.elementTable)
        #ratio = 0.4 * dryAir.density / iodineWater.density
        #mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        #mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        #mix.Commit()
        #mix.CalculateMixtureDensity()
        #self.extraMaterialList.append(mix)

        # 50% dry air + 50% (iodine + water)
        mix = mixture.Mixture("air(50%)-iodine-water", self.elementTable)
        ratio = 0.5 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        ## 60% dry air + 40% (iodine + water)
        #mix = mixture.Mixture("air(60%)-iodine-water", self.elementTable)
        #ratio = 0.6 * dryAir.density / iodineWater.density
        #mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        #mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        #mix.Commit()
        #mix.CalculateMixtureDensity()
        #self.extraMaterialList.append(mix)

        # 70% dry air + 30% (iodine + water)
        mix = mixture.Mixture("air(70%)-iodine-water", self.elementTable)
        ratio = 0.7 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

        ## 80% dry air + 20% (iodine + water)
        #mix = mixture.Mixture("air(80%)-iodine-water", self.elementTable)
        #ratio = 0.8 * dryAir.density / iodineWater.density
        #mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        #mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        #mix.Commit()
        #mix.CalculateMixtureDensity()
        #self.extraMaterialList.append(mix)

        # 90% dry air + 10% (iodine + water)
        mix = mixture.Mixture("air(90%)-iodine-water", self.elementTable)
        ratio = 0.9 * dryAir.density / iodineWater.density
        mix.AddMaterial(dryAir, materialWeightFraction = ratio)
        mix.AddMaterial(iodineWater, materialWeightFraction = 1.0 - ratio)
        mix.Commit()
        mix.CalculateMixtureDensity()
        self.extraMaterialList.append(mix)

#------------------------------------------------------------
#------------------------------------------------------------
class ControlPointManagerKidneyStoneExtra(ControlPointManager):
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, com, method, isToDumpIntermediateData = True):
        super(ControlPointManagerKidneyStoneExtra, self).__init__(com, method, isToDumpIntermediateData)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def AddExtraMaterial2(self):
        self.AddManuallyTunedExtraMaterial()

        # uric acid
        # https://byjus.com/uric-acid-formula/
        # C5 H4 N4 O3
        uricAcid = material.Material("uric acid", self.elementTable)
        uricAcid.AddElement(6, atomicFraction = 5)
        uricAcid.AddElement(1, atomicFraction = 4)
        uricAcid.AddElement(7, atomicFraction = 4)
        uricAcid.AddElement(8, atomicFraction = 3)
        uricAcid.Commit()
        uricAcid.density = 1.87
        self.extraMaterialList.append(uricAcid)

        # calcium oxalate monohydrate
        # https://www.americanelements.com/calcium-oxalate-monohydrate-5794-28-5
        # Ca C2 H2 O5
        calciumOxMono = material.Material("calcium oxalate monohydrate", self.elementTable)
        calciumOxMono.AddElement(20, atomicFraction = 1)
        calciumOxMono.AddElement(6 , atomicFraction = 2)
        calciumOxMono.AddElement(1 , atomicFraction = 2)
        calciumOxMono.AddElement(8 , atomicFraction = 5)
        calciumOxMono.Commit()
        calciumOxMono.density = 2.2
        self.extraMaterialList.append(calciumOxMono)

        # calcium oxalate dihydrate
        # https://www.americanelements.com/calcium-oxalate-monohydrate-5794-28-5
        # Ca C2 H4 O6
        calciumOxDi = material.Material("calcium oxalate dihydrate", self.elementTable)
        calciumOxDi.AddElement(20, atomicFraction = 1)
        calciumOxDi.AddElement(6 , atomicFraction = 2)
        calciumOxDi.AddElement(1 , atomicFraction = 4)
        calciumOxDi.AddElement(8 , atomicFraction = 6)
        calciumOxDi.Commit()
        calciumOxDi.density = 2.2
        self.extraMaterialList.append(calciumOxDi)

        # hydroxyapatite
        # http://cameo.mfa.org/wiki/Calcium_hydroxyapatite
        # https://pubmed.ncbi.nlm.nih.gov/23263603/
        # Ca10 H2 O26 P6
        # CT number exceeds 3095 with density 3.18
        # manually tune it down
        scaleFactorList = [0.5]
        for idx in range(len(scaleFactorList)):
            hydroxyapatite = material.Material("hydroxyapatite (Var {:d})".format(idx), self.elementTable)
            hydroxyapatite.AddElement(20 , atomicFraction = 10)
            hydroxyapatite.AddElement(1  , atomicFraction = 2)
            hydroxyapatite.AddElement(8  , atomicFraction = 26)
            hydroxyapatite.AddElement(15 , atomicFraction = 6)
            hydroxyapatite.Commit()
            hydroxyapatite.density = 3.18
            hydroxyapatite.density *= scaleFactorList[idx]
            self.extraMaterialList.append(hydroxyapatite)

        # cystine
        # https://www.sigmaaldrich.com/catalog/product/mm/102836?lang=en&region=US
        # C6 H12 N2 O4 S2
        cystine = material.Material("cystine", self.elementTable)
        cystine.AddElement(6  , atomicFraction = 6 )
        cystine.AddElement(1  , atomicFraction = 12)
        cystine.AddElement(7  , atomicFraction = 2 )
        cystine.AddElement(8  , atomicFraction = 4 )
        cystine.AddElement(16 , atomicFraction = 2 )
        cystine.Commit()
        cystine.density = 1.677
        self.extraMaterialList.append(cystine)

        # struvite
        # https://pubmed.ncbi.nlm.nih.gov/23263603/
        # Ca Mg N H16 P O10
        struvite = material.Material("struvite", self.elementTable)
        struvite.AddElement(20 , atomicFraction = 1 )
        struvite.AddElement(12 , atomicFraction = 1 )
        struvite.AddElement(7  , atomicFraction = 1 )
        struvite.AddElement(1  , atomicFraction = 16)
        struvite.AddElement(15 , atomicFraction = 1 )
        struvite.AddElement(8  , atomicFraction = 10)
        struvite.Commit()
        struvite.density = 1.711
        self.extraMaterialList.append(struvite)

