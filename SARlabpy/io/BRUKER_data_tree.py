# -*- coding: utf-8 -*-
"""
BRUKER_data_tree based on the example in pyqtgraph for the simple use of 
DataTreeWidget to display a structure of nested dicts, lists, and arrays
"""

import os
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import numpy as np

import SARlabpy.io.BRUKER_classes as sar 
import SARlabpy

app = QtGui.QApplication([])

pdata  = sar.PDATA_file(os.path.join(sar.dataroot,'readfidTest.ix1',
                                        '3','pdata','1'))
           
dicts = {'filename':pdata.filename, 
         'reco':pdata.reco,
         'd3proc':pdata.d3proc}                           
tree = pg.DataTreeWidget(data=dicts)
tree.show()
tree.setWindowTitle('pyqtgraph: ACQP tree')
tree.resize(600,600)


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()