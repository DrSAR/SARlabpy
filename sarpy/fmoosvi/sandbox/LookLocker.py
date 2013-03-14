# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 22:41:46 2013

@author: firas
"""

## Porting over LookLocker code from Annie

import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters
import sarpy.fmoosvi.wrappers


NecS3 = sarpy.Experiment('NecS3')

LLdata = NecS3.find_scan_by_protocol('04_ubcLL2')