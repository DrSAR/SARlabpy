# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 00:49:44 2013

@author: stefan
"""

import sarpy

exp = sarpy.Experiment('*')
print exp.__repr__()
for study in exp.studies:
    for scan in study.scans:
        for pdata in scan.pdata:
            try:
                if pdata.reco.RECO_base_image_uid != pdata.visu_pars.VisuUid:
                    print('not consistent === ',scan)
                else:
                    print(scan)
            except IOError:
                print('not found: ',scan)
