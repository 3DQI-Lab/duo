import sys
import os
result = os.path.join("..", "..")
sys.path.append(result)

from PySide2 import QtCore, QtGui, QtWidgets
import duo.core.duo_exception as de
import duo.zeff.cai as cai
import duo.gui.progress as progress

#------------------------------------------------------------
#------------------------------------------------------------
class Worker(QtCore.QThread):
    #------------------------------------------------------------
    #------------------------------------------------------------
    signalTotalNumImage = QtCore.Signal(int, int)
    signalImageIdx      = QtCore.Signal(int)

    #------------------------------------------------------------
    # run on old thread
    #------------------------------------------------------------
    def __init__(self, inputDir, outputDir):
        QtCore.QThread.__init__(self)

        self.inputDir       = inputDir
        self.outputDir      = outputDir



    #------------------------------------------------------------
    # run on new thread
    #------------------------------------------------------------
    def run(self):
        dataPath = os.path.join("..", "..", "data")
        c = cai.Cai(dataPath, isToDumpIntermediateData = False)
        it = c.ApplyModelRecursively(self.inputDir, self.outputDir)
        totalNumImage = next(it)
        self.signalTotalNumImage.emit(0, totalNumImage)

        # consume the rest of generator
        while True:
            try:
                imageIdx = next(it)
                self.signalImageIdx.emit(imageIdx)
            except StopIteration:
                break

#------------------------------------------------------------
#------------------------------------------------------------
class DuoGUI:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        self.app = None
        self.mainWindow = None
        self.inputDir = None
        self.outputDir = None

        self.inputDirLineEdit = None
        self.outputDirLineEdit = None

        self.inputDirButton = None
        self.outputDirButton = None

        self.runButton = None
        # self.abortButton = None

        self.worker = None

        self.Initialize()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Initialize(self):
        self.app = QtWidgets.QApplication(sys.argv)
        self.desktopGeometry = QtWidgets.QApplication.desktop().screenGeometry()
        self.screenCenter = QtWidgets.QDesktopWidget().availableGeometry().center()

        self.mainWindow = QtWidgets.QMainWindow()
        self.mainWindow.setWindowTitle("duo")
        self.mainWindow.setMinimumSize(480, 300)
        self.mainWindow.resize(480, 300)

        vbox = QtWidgets.QVBoxLayout()

        hbox = QtWidgets.QHBoxLayout()
        self.inputDirLineEdit = QtWidgets.QLineEdit()
        hbox.addWidget(self.inputDirLineEdit)
        self.inputDirButton = QtWidgets.QPushButton("Specify directory")
        hbox.addWidget(self.inputDirButton)
        gbox = QtWidgets.QGroupBox("DICOM directory (Input. Must exist.)")
        gbox.setLayout(hbox)
        vbox.addWidget(gbox)

        hbox = QtWidgets.QHBoxLayout()
        self.outputDirLineEdit = QtWidgets.QLineEdit()
        hbox.addWidget(self.outputDirLineEdit)
        self.outputDirButton = QtWidgets.QPushButton("Specify directory")
        hbox.addWidget(self.outputDirButton)
        gbox = QtWidgets.QGroupBox("Zeff image directory (Output. Will be created if not existing)")
        gbox.setLayout(hbox)
        vbox.addWidget(gbox)

        hbox = QtWidgets.QHBoxLayout()
        self.runButton = QtWidgets.QPushButton("Run")
        # self.abortButton = QtWidgets.QPushButton("Abort")
        # self.abortButton.setEnabled(False)
        hbox.addWidget(self.runButton)
        # hbox.addWidget(self.abortButton)
        temp = QtWidgets.QWidget()
        temp.setLayout(hbox)
        vbox.addWidget(temp)

        centralWidget = QtWidgets.QWidget()
        centralWidget.setLayout(vbox)
        self.mainWindow.setCentralWidget(centralWidget)
        self.mainWindow.show()

        self.inputDirButton.clicked.connect(lambda: self.SetUpFileDialog(self.inputDirLineEdit))
        self.outputDirButton.clicked.connect(lambda: self.SetUpFileDialog(self.outputDirLineEdit))
        self.runButton.clicked.connect(lambda: self.Run())
        # self.abortButton.clicked.connect(lambda: self.Abort())

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Loop(self):
        # enter main event loop and wait until exit() is called
        sys.exit(self.app.exec_())

    #------------------------------------------------------------
    #------------------------------------------------------------
    def SetUpFileDialog(self, lineEdit):
        dialog = QtWidgets.QFileDialog()
        dialog.setFileMode(QtWidgets.QFileDialog.Directory)

        if dialog.exec_():
            dirPath = dialog.selectedFiles()
            if len(dirPath) > 0:
                lineEdit.setText(dirPath[0])

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Run(self):
        try:
            self.inputDir = self.inputDirLineEdit.text()
            # the specified path does not exist or is not a directory
            if (not os.path.exists(self.inputDir)) or (not os.path.isdir(self.inputDir)):
                raise de.DuoException("Invalid DICOM directory.")

            self.outputDir = self.outputDirLineEdit.text()
            # the specified path exists but is a file instead of a directory
            if os.path.exists(self.outputDir) and (not os.path.isdir(self.inputDir)):
                raise de.DuoException("Invalid Zeff image directory.")

            if not os.path.exists(self.outputDir):
                os.mkdir(self.outputDir)

            self.progressWindow = progress.ProgressWindow()
            self.worker = Worker(self.inputDir, self.outputDir)
            self.worker.started.connect(lambda : self.OnThreadStart())
            self.worker.finished.connect(lambda : self.OnThreadFinish())
            self.worker.signalTotalNumImage.connect(self.progressWindow.SetMinMax)
            self.worker.signalImageIdx.connect(self.progressWindow.SetProgress)
            self.worker.start() # non-blocking execution !!!

        except de.DuoException as err:
            print(err)
            QtWidgets.QMessageBox.critical(self.mainWindow, "Error", str(err))

        except:
            errMsg = str(sys.exc_info()[1])
            print("Unexpected error:", errMsg)
            QtWidgets.QMessageBox.critical(self.mainWindow, "Error", errMsg)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def OnThreadStart(self):
        self.runButton.setEnabled(False)
        # self.abortButton.setEnabled(True)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def OnThreadFinish(self):
        self.runButton.setEnabled(True)
        # self.abortButton.setEnabled(False)
        if self.progressWindow:
            self.progressWindow.close()
            self.progressWindow = None
        QtWidgets.QMessageBox.information(self.mainWindow, "Info", "Calculation has been completed successfully.")

    # #------------------------------------------------------------
    # #------------------------------------------------------------
    # def Abort(self):
        # self.runButton.setEnabled(True)
        # self.abortButton.setEnabled(False)
