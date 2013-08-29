# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 15:43:31 2013

@author: stefan
"""

import sarpy
import pandas
import parsing

def df_from_masterlist(masterlist):
    '''
    Take in a masterlist (which is a dict after reading from json file)
    and parse patient names into 
        - patient name
        - tumour type
        - tumour location
        - patient number
    '''
    raise NotImplemented

def calculate_dfSeries_from_adata(
        masterlist,
        roots_scanlabels,
        adata_labels,
        roi_labels):
    '''
    Take in previously calculated adata and work our some 
    descriptive summary data.
    
    Returns: pandas.Dataframe
    '''
    
    raise NotImplemented

    
