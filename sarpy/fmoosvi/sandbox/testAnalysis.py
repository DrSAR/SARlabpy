# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 00:24:57 2013

@author: firas
"""
# Add sarpy to the di

import sarpy
import sarpy.fmoosvi.wrappers

reload(sarpy)

NecS1Exp = sarpy.Experiment('NecS1')
NecS1dce = NecS1Exp.find_scan_in_experiment('06_FLASH2D.6sl-DCE')


#Experiment
auc = sarpy.fmoosvi.wrappers.calculate_AUC(NecS1Exp)



##Single Study #works!
#auc = sarpy.fmoosvi.wrappers.calculate_AUC(NecS1Exp.studies[0])
#
##Single Scan #works!
#auc = sarpy.fmoosvi.wrappers.calculate_AUC(NecS1Exp.studies[0].scans[7])
#
##List of Studies #works!
#auc = sarpy.fmoosvi.wrappers.calculate_AUC(NecS1Exp.studies[0:5])
#
##List of Scans #works!
#NecS1auc = sarpy.fmoosvi.wrappers.calculate_AUC(NecS1dce)
