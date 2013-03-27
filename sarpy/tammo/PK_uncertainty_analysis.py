# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 10:58:28 2013

@author: tammo
"""

# Setting up for data input and analysis
from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import pywt
import scipy


import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.wrappers

Experiment = sarpy.Experiment('NecS3')

fid = Experiment.studies[1].scans[5].fid

auc = sarpy.fmoosvi.wrappers.calculate_AUC(Experiment.studies[0])