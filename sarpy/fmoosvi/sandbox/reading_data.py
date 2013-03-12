# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 22:51:48 2013

@author: firas
"""

import sarpy

readfidExp = sarpy.Experiment('readfid')
for study in readfidExp.studies:
    print('-'*40+'\n'+study.subject.SUBJECT_id)
    for scan in study.scans:
        print("  "+scan.acqp.ACQ_protocol_name)



