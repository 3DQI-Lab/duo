from PySide2 import QtCore, QtGui, QtWidgets

#------------------------------------------------------------
#------------------------------------------------------------
class ProgressWindow(QtWidgets.QMainWindow):
    #------------------------------------------------------------
    # minValue, maxValue must be int
    #------------------------------------------------------------
    def __init__(self):
        super().__init__()
        self.progressBar = None

        self.setWindowTitle("Progress")
        self.setMinimumSize(480, 50)
        self.resize(480, 50)

        self.progressBar = QtWidgets.QProgressBar()
        self.setCentralWidget(self.progressBar)
        self.show()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def SetMinMax(self, minValue, maxValue):
        self.progressBar.setMinimum(minValue)
        self.progressBar.setMaximum(maxValue)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def SetProgress(self, value):
        self.progressBar.setValue(value)