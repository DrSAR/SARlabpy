# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 18:24:59 2013

@author: fmoosvi
"""

## This is to compare B1 maps from Flash, MSME and VFA.

import sarpy
import sarpy.fmoosvi.analysis
import pylab

########################## Start Example #########################
### This is for a phantom study using LL, BSB1 MSME, BSB1 flash
#
#AP1Exp = sarpy.Experiment('DiLL.iu2')
#BS_scans = AP1Exp.find_scan_by_protocol('07\_bSB1mapFLASH\(modified\)')
#
#zero_BS_minus = BS_scans[0]
#zero_BS_plus = BS_scans[1]
#power_BS_minus = BS_scans[2]
#power_BS_plus = BS_scans[3]
#
#b1map = sarpy.fmoosvi.analysis.h_BS_B1map(zero_BS_minus,\
#                                          zero_BS_plus,\ 
#                                          power_BS_minus,\
#                                          power_BS_plus,\
#                                          pdata_num=0)
#                                          
#pylab.imshow(b1map[:,:,0,0], aspect = 0.5)
########################## End Example #########################