# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 21:51:39 2013

@author: firas
"""

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


######################## Start Example #########################
# This is for one phantom study using LL an BSB1 MSME
APExp = sarpy.Experiment('DiLL.iI1')
BS_scans = APExp.find_scan_by_protocol('09\-2DBSB1map\-MSME_bas\(modified\)')[-4:]


zero_BS_minus = BS_scans[0]
zero_BS_plus = BS_scans[1]
power_BS_minus = BS_scans[2]
power_BS_plus = BS_scans[3]

b1map = sarpy.fmoosvi.analysis.h_BS_B1map(zero_BS_minus,\
                                          zero_BS_plus,\
                                          power_BS_minus,\
                                          power_BS_plus,\
                                          pdata_num=0)
                                          
pylab.imshow(b1map[:,:,0,0], aspect = 0.5)
######################## End Example #########################

from sarpy.fmoosvi import analysis

minus_image = analysis.h_phase_from_fid(power_BS_minus)
plus_image =  analysis.h_phase_from_fid(power_BS_plus)

fig1 = pylab.figure()
pylab.imshow(minus_image[:,:,0])
pylab.colorbar()
fig = pylab.figure()
pylab.imshow(plus_image[:,:,0])
pylab.colorbar()



#APExp = sarpy.Experiment('dBloch')
#
#zero_BS_minus = APExp.studies[0].scans[0]  # -4000 45
#zero_BS_plus = APExp.studies[0].scans[1] #  4000 46
#power_BS_minus = APExp.studies[0].scans[2]     #-4000 35
#power_BS_plus = APExp\.studies[0].scans[3]   # 4000 36
#
#
#b1map = sarpy.fmoosvi.analysis.h_BS_B1map(zero_BS_minus,\
#                                          zero_BS_plus,\
#                                          power_BS_minus,\
#                                          power_BS_plus,\
#                                          pdata_num=1)
#                                          
#pylab.imshow(b1map[:,:,20])