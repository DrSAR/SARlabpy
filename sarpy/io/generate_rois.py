#!/usr/bin/env python 
# -*- coding: utf-8 -*-

"""
Created on Thu May 30 17:49:35 2013

@author: fmoosvi
"""

import sarpy.fmoosvi.analysis
import argparse
import os

#TODO: for some reason the required args still show up as optional

parser = argparse.ArgumentParser()

parser.add_argument('-m', '--masterlist_name', type=str, required=True,
                   help='Name of the masterlist file. e.g. NecS3 \
                   Usage: python generate_rois.py -m HerP2 -d axref -a roi -i export ')

parser.add_argument('-d', '--data_label', type=str, required=True,
                   help='Data label, usually anatomy or IR_A')                   

parser.add_argument('-i', '--iotype', type=str, required=True,
                   choices = ['import','export'], help='IOType, import or export')    

parser.add_argument('-a', '--adata_label', type=str,
                   help='Optional: adata label, defaults to roi') 
                                     
parser.add_argument('-p', '--path', type=str,
                   help='Optional: Path to/from export/import data, default to mdata')
                   
parser.add_argument('-f', '--force', type=str, choices = ['True','False'],
                    help='Optional: Replace data? True or False, defaults to False')   
                   
args = parser.parse_args()
masterlist_name = args.masterlist_name
data_label = args.data_label
ioType = args.iotype

#TODO: Integrate this intothe parser for cleaner code!

try:
    adata_label = args.adata_label
except AttributeError:
    adata_label = 'roi'

try:
    path = os.path.expanduser(os.path.join(args.path))
except AttributeError:
    path = os.path.expanduser(os.path.join('~','mdata',masterlist_name))

try:
    force = args.force
    force = True
except AttributeError:
    force = False

sarpy.fmoosvi.analysis.h_generate_ROI(masterlist_name, data_label, 
                                    adata_label, ioType, path, forceVal = force)