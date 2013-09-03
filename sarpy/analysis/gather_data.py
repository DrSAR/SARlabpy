 # -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 15:43:31 2013

@author: stefan
"""

import sarpy
import pandas
from pyparsing import (Word, alphas, nums, Group, SkipTo, StringEnd)
import re
import collections
import os
import json
import scipy
import scipy.stats
import numpy


def df_from_masterlist(masterlist_name):
    '''
    Take in a masterlist (which is a dict after reading from json file)
    and parse patient names into 
        - patient name
        - tumour type
        - tumour location
        - patient number
    '''
    root = os.path.expanduser(os.path.join('~/sdata',
                        masterlist_name,
                        masterlist_name))
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)    
                           
                           
    studyname = Word(alphas)+Word(nums)
    tumourtype = Word(alphas, exact=1)
    location = Word(alphas, exact=1)
    patientnumber = Word(nums)
    patientname = (studyname('studyname') + tumourtype('tumourtype') + 
                   location('location') + patientnumber('patientnumber'))
   
   
    patient_df = pandas.DataFrame(columns=('studyname', 'patientnumber'))
    for patname in master_list.keys():
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
        masterlist_name,
        scan_label_list,
        adata_label,
        roi_label
        ):
    '''
    Take in previously calculated adata and work out some 
    descriptive summary data.
    
    Returns: pandas.Dataframe
    '''
   # Open the masterlist and read in all the crap

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

# This initializes the list of data frames that will be generated                           
    df = []                            

## Here there will be a confusing loop over all scan_labels, if scan_label is
# just a string, make it into a list just to make things easier for us later
# the loop doesn't work for strings

    if isinstance(scan_label_list,str):
        scan_label_list = []
        scan_label_list.append(scan_label)

    for scan_label in scan_label_list:

        # Set the column name for each fo the adata_labels in the list
            # get the prefix of the datalabel, which corresponds to, 
            # 0h-, 24h-, +24h- etc...
    
        prefix = re.match(".*-",scan_label) 
        dfSeries_label = prefix.group() + adata_label
    
        # Initialize the dictionaries in which summary data is stored
        avg_data = collections.OrderedDict()
        slice_position = collections.OrderedDict()
        roi_px_count = collections.OrderedDict()
        slices = collections.OrderedDict()
    
        animal = []
        slc = []
        avg = []
        px = []
        pos = []
    
        # Loop through each animal
        for k,v in master_list.iteritems():
            
            try:
                scan = sarpy.Scan(v[scan_label][0])           
                data = scan.adata[adata_label].data
                roi = sarpy.Scan(v[scan_label][0]).adata[roi_label].data
            except(KeyError, IOError):
                print('{0}: Something not found or {1} in patient {2}'.format(adata_label,scan_label,k) )
                continue
           
            roi_data = data * roi
            
            # Get the slices and patient names and repeat patient names for the slice
            # Create the data for the data series
            
            #slice_position.update(dict(zip([k]*len(slices),scan.acqp.ACQ_slice_offset)))
            
            for slice in xrange(roi_data.shape[-1]):
    
                animal.append(k)
                slc.append(slice+1)
                avg.append(scipy.stats.nanmean(roi_data[:,:,slice].flatten()))
                px.append(numpy.nansum(roi[:,:,slice].flatten()))
                pos.append(scan.acqp.ACQ_slice_offset[slice])          
    
        series_of_dicts = {dfSeries_label  + '_slice' : pandas.Series(slc, index = numpy.array(animal)),
                           dfSeries_label + '_avg': pandas.Series(avg, index = numpy.array(animal)),                       
                           dfSeries_label + '_pos': pandas.Series(pos, index = numpy.array(animal)),
                           dfSeries_label + '_px': pandas.Series(px, index = numpy.array(animal))}
                  
        df.append(pandas.DataFrame(series_of_dicts))
        
    return pandas.concat(df)