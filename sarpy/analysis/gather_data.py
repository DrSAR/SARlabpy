# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 15:43:31 2013

@author: stefan
"""

import sarpy
import pandas
from pyparsing import (Word, alphas, nums, Group, SkipTo, StringEnd)


def df_from_masterlist(masterlist):
    '''
    Take in a masterlist (which is a dict after reading from json file)
    and parse patient names into 
        - patient name
        - tumour type
        - tumour location
        - patient number
    '''

    studyname = Word(alphas)+Word(nums)
    tumourtype = Word(alphas, exact=1)
    location = Word(alphas, exact=1)
    patientnumber = Word(nums)
    patientname = (studyname('studyname') + tumourtype('tumourtype') + 
                   location('location') + patientnumber('patientnumber'))
   
   
    patient_df = pandas.DataFrame(columns=('studyname', 'patientnumber'))
    for patname in masterlist.keys():
        print patname
#        print line
        parsed_data = patientname.parseString(patname)
        row = pandas.DataFrame([{'studyname':''.join(parsed_data['studyname']),
                                 'tumourtype':parsed_data['tumourtype'],
                                 'location':parsed_data['location'],
                                 'patientnumber':parsed_data['patientnumber']}],
                                 index=(patname,))
        patient_df = patient_df.append(row)

    return patient_df

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

    
