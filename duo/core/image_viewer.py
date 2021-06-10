import sys
import os
import sys
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
plt.ion()

from matplotlib.patches import Circle
from mpl_toolkits.mplot3d import Axes3D
import scipy
import scipy.optimize
import pydicom
from matplotlib import cm


#------------------------------------------------------------
# This image viewer plots 2D data stored in .npy format.
# Usage:
# --- Create a python script that instantiates image_viewer object
# --- Pass the full path of the .npy file to the script
#
# import os
# import sys
# result = os.path.join("..", "..", "..")
# sys.path.append(result)
# print("--> sys.path = ", sys.path)
# import duo.core.image_viewer as iv
#
# #------------------------------------------------------------
# #------------------------------------------------------------
# if __name__ == '__main__':
#     iv = ImageViewer(sys.argv)
#     iv.Plot()
#------------------------------------------------------------
class ImageViewer:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self, argv): # psargs is positional arguments (a tuple)
        print(argv)
        print(">>>>>")
        self.filename = argv[1]

        self.vmin = None
        self.vmax = None
        if len(argv) == 4:
            self.vmin = argv[2]
            self.vmax = argv[3]

        self.title = None
        if len(argv) == 5:
            self.title = argv[4]

        self.data = np.load(self.filename)
        self.fig = []

        self.circle = []
        self.text = []

        self.ax = []

        self.col_i = 0
        self.row_i = 0
        self.z = 0.0

        print("    Data type = ", self.data.dtype)
        print("    Max = ", np.amax(self.data))
        print("    Min = ", np.amin(self.data))
        print("    Shape = ", self.data.shape)

        self.DisableDefaultKeyMap()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Plot(self):
        self.fig = plt.figure(figsize = (8, 8))
        self.ax = self.fig.add_subplot(1, 1, 1)

        if self.vmin != None and self.vmax != None:
            m = self.ax.imshow(self.data, interpolation = 'none', cmap = cm.gray,
                               vmin = self.vmin, vmax = self.vmax)
        else:
            m = self.ax.imshow(self.data, interpolation = 'none', cmap = cm.gray)

        self.circle = Circle(xy = (10.0, 10.0), radius = 4.0, fill = False, color = '#ff0000')
        self.ax.add_artist(self.circle)

        self.text = plt.text(x = 10.0, y = 10.0, s = "text here", color = '#ffff00', fontsize = 12,
                             bbox = dict(facecolor = '#000000', alpha = 0.5))

        self.fig.canvas.mpl_connect('button_press_event', self.OnClick)
        self.fig.canvas.mpl_connect('key_press_event', self.OnKeyPress)

        if self.title != None:
            self.ax.set_title(self.title)

        plt.show(block = True)

    #------------------------------------------------------------
    #------------------------------------------------------------
    def OnClick(self, event):
        if event.inaxes is None:
                    return

        self.col_i = int(np.rint(event.xdata))
        self.row_i = int(np.rint(event.ydata))

        self.UpdateInfo()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def OnKeyPress(self, event):
        if event.inaxes is None:
                    return

        if (event.key == "up" or event.key == "w") and self.row_i > 0:
            self.row_i -= 1

        if (event.key == "down" or event.key == "s") and self.row_i < self.data.shape[0] - 1:
            self.row_i += 1

        if (event.key == "left" or event.key == "a") and self.col_i > 0:
            self.col_i -= 1

        if (event.key == "right" or event.key == "d") and self.col_i < self.data.shape[1] - 1:
            self.col_i += 1

        self.UpdateInfo()

    #------------------------------------------------------------
    #keymap.fullscreen : f, ctrl+f       # toggling
    #keymap.home : h, r, home            # home or reset mnemonic
    #keymap.back : left, c, backspace    # forward / backward keys to enable
    #keymap.forward : right, v           #   left handed quick navigation
    #keymap.pan : p                      # pan mnemonic
    #keymap.zoom : o                     # zoom mnemonic
    #keymap.save : s                     # saving current figure
    #keymap.quit : ctrl+w, cmd+w         # close the current figure
    #keymap.grid : g                     # switching on/off a grid in current axes
    #keymap.yscale : l                   # toggle scaling of y-axes ('log'/'linear')
    #keymap.xscale : L, k                # toggle scaling of x-axes ('log'/'linear')
    #------------------------------------------------------------
    def DisableDefaultKeyMap(self):
        plt.rcParams['keymap.fullscreen'] = ''
        plt.rcParams['keymap.home'      ] = ''
        plt.rcParams['keymap.back'      ] = ''
        plt.rcParams['keymap.forward'   ] = ''
        plt.rcParams['keymap.pan'       ] = ''
        plt.rcParams['keymap.zoom'      ] = ''
        plt.rcParams['keymap.save'      ] = ''
        plt.rcParams['keymap.quit'      ] = ''
        plt.rcParams['keymap.grid'      ] = ''
        plt.rcParams['keymap.yscale'    ] = ''
        plt.rcParams['keymap.xscale'    ] = ''

    #------------------------------------------------------------
    #------------------------------------------------------------
    def UpdateInfo(self):
        self.z = self.data[self.row_i, self.col_i]

        msg = "col_i = {0:d}, row_i = {1:d}, z = {2:f}".format(self.col_i, self.row_i, self.z)
        print(msg)

        self.circle.center = self.col_i, self.row_i

        msg = "col index = {0:d}\nrow index = {1:d}\nvalue = {2:f}".format(self.col_i, self.row_i, self.z)
        self.text.set_text(msg)
        self.text.set_position((self.col_i + 10, self.row_i - 10))

        self.fig.canvas.draw()






#------------------------------------------------------------
#------------------------------------------------------------
if __name__ == '__main__':
    iv = ImageViewer(sys.argv)
    iv.Plot()





