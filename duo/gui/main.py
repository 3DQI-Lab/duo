from PySide2 import QtCore, QtGui, QtWidgets

import sys
import os
result = os.path.join("..", "..")
sys.path.append(result)

import gui
import duo.core.duo_exception as de

#------------------------------------------------------------
#------------------------------------------------------------
if __name__ == '__main__':
    try:
        duoGUI = gui.DuoGUI()
        duoGUI.Loop()
    except de.DuoException as e:
        print(e)
        input()

