# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 01:06:33 2013

@author: stefan
"""
from __future__ import print_function

from PyQt4 import QtGui  # (the example applies equally well to PySide)
import pyqtgraph as pg
import numpy
from matplotlib.cbook import flatten
import SARlabpy
import SARlabpy.io.BRUKER_classes as sar
import SARlabpy.io.SARlogger
SARlabpy.io.SARlogger.initiate_logging(sar)

#read all the data
alldata = []
allnames = []
current_data = 0        

NecS3Exp = sar.Experiment('NecS3')
for study in NecS3Exp.studies:
    print('-'*40+'\n'+study.subject.SUBJECT_id)
    for scan in study.scans:
        print(' '+scan.acqp.ACQ_protocol_name)
        if scan.acqp.ACQ_protocol_name == '06_FLASH2D.6sl-DCE':
            data = scan.pdata[0].data
            print(data.shape)
            alldata.append(data)
            allnames.append(scan.acqp.FILE_LOCATION)


class QButtonAdv(QtGui.QWidget):
    def __init__(self, imv=None, lbl=None, parent=None):
        self.imv=imv
        self.lbl=lbl
        QtGui.QWidget.__init__(self, parent)
        self.button = QtGui.QPushButton('advance', self)
        self.button.clicked.connect(self.advance_image)
    def advance_image(self):
        global current_data
        data = alldata[current_data]
        print(current_data, data.shape,allnames[current_data],sep=' - ')
        self.imv.setImage(data)
        self.lbl.setText(allnames[current_data])
        current_data = (current_data + 1) % len(alldata)

class QButtonRef(QtGui.QWidget):
    def __init__(self, imv=None, parent=None):
        self.imv=imv
        QtGui.QWidget.__init__(self, parent)
        self.button = QtGui.QPushButton('reformat', self)
        self.button.clicked.connect(self.reformat)
    def reformat(self):
        data = alldata[current_data]
        NR = data.shape[0]/6
        idx_nested = [numpy.arange(NR)*6+i for i in numpy.arange(6)]
        idx = [i for i in flatten(idx_nested)]
        newdata = data[idx, :, :]
        self.imv.setImage(newdata)


## Always start by initializing Qt (only once per application)
app = QtGui.QApplication([])

## Define a top-level widget to hold everything
w = QtGui.QWidget()

## Create some widgets to be placed inside
imv = pg.ImageView()
imv.setImage(alldata[0])
lbl = QtGui.QLabel()
btn1 = QButtonAdv(imv, lbl)
btn2 = QButtonRef(imv)

## Create a grid layout to manage the widgets size and position
layout = QtGui.QGridLayout()
w.setLayout(layout)

## Add widgets to the layout in their proper positions
layout.addWidget(btn1, 0, 0)   # button goes in upper-left
layout.addWidget(lbl, 1, 0)   # button goes in upper-left
layout.addWidget(btn2, 2, 0)   # button goes in upper-left
layout.addWidget(imv, 0, 1, 3, 1)  # plot goes on right side, spanning 3 rows

## Display the widget as a new window
w.show()

## Start the Qt event loop
app.exec_()

