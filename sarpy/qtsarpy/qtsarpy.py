# -*- coding: utf-8 -*-
"""
qt-based image data viewer with pyqtgrpah to do the heavy lifting
in the display department
"""

from qtsarpy_mainwindow_auto import Ui_MainWindow
# importing PyQt modules from ui generated file is a round-about way 
# of trying to follow the choice of PyQt vs PySide that's been made there... 
from qtsarpy_mainwindow_auto import QtGui, QtCore
import os

# This code should be in the auto-generated file but since
# pyside-uic does not understand switch -w it can only be generated
# if pyuic4 is present. Hence we move this in here to become 
# independent of pyside-uic / pyuic4 compilatino of the *ui file.
class MainWindow(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None, f=QtCore.Qt.WindowFlags()):
        QtGui.QMainWindow.__init__(self, parent, f)
        self.setupUi(self)


class viewer(MainWindow):
    def __init__(self, parent=None, name=None, fl=0):
        MainWindow.__init__(self)

    def on_actionChange_dataroot_triggered(self, checked=None):
        if checked is None: return
        fname = QtGui.QFileDialog.getExistingDirectory(
                    self, 
                    'Root of Source Data', 
                    os.path.expanduser('~'))
        print fname
                
if __name__ == "__main__":
    import sys
    a = QtGui.QApplication(sys.argv)
    QtCore.QObject.connect(a,QtCore.SIGNAL("lastWindowClosed()"),
                           a,QtCore.SLOT("quit()"))
    vwr = viewer()
    vwr.show()
    a.exec_()
 
 
