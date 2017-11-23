#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 13:53:03 2013

@author: analysis
"""

import json
import os
import collections
import csv
import scipy.stats
import sarpy
import time
import getpass # used to get the current username
import numpy
import sarpy.analysis.getters


def write_csv(masterlist_name, 
              data_label, 
              adata_label,
              roi_override = None):

    import time

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    fname_to_open = root+'.json'
    
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str) 

    # So write csv looks at the right adata label to get the roi. primarily for deltaT1 adata
    if roi_override is None:
        roi_override = adata_label+'_roi'        
    
    # Get Max number of slices, this is so that the average column is always the last column
    # and so slice averages get zero-filled

    numSlices = []
    for k,v in list(master_list.items()):

        try:        
            data = sarpy.Scan(v[data_label][0]).adata[adata_label].data
            roi = sarpy.Scan(v[data_label][0]).adata[roi_override].data  
            try:
                numSlices.append(data.shape[2])
            except IndexError:
                numSlices.append(1)
        except(IOError,KeyError):
            continue

    # Create empty list so that it can be populated with rows
    export_data = []

    # Find the largest slice - really not sure why this is done this way
    #if numSlices:
    #    maxSlices = numpy.max(numSlices)
    
    for k,v in list(master_list.items()):

        scn_to_analyze = v[data_label][0]

        avgL = adata_roi_average(scn_to_analyze,adata_label,roi_override)
        
        ######     

        export_data.append([k] + avgL)

    time = time.strftime("%Y-%m-%d %H:%M")       				       
    header =['Animal ID'] + ['Slice '+str(x+1) for x in range(
                                        len(max(export_data,key=len))-3)] + ['Average'] + ['Volume']

    footers = []
    footers.append(['# Data generated on {0} by {1}'.format(time,getpass.getuser())])
    footers.append(["#Legend"])
    footers.append(["#inf:slice wasn't acquired"])
    footers.append(["#nan: no roi in slice"])
    footers.append(["#AnalysisErr: No Data available"])
    
    for f in footers:
        export_data.append(f)    
    

    filename = os.path.expanduser(os.path.join(
                                          '~','sdata',masterlist_name,'export',
                                          data_label+'_'+adata_label+'.csv'))
   
    # open output file
    outfile = open(filename, "w" )
    
    # get a csv writer
    writer = csv.writer(outfile)
    
    # write header
    writer.writerow(header)
    
    # write data
    [ writer.writerow(x) for x in export_data ]
    
    # close file
    outfile.close()

def adata_roi_average(scn_to_analyze,
                      adata_label,
                      roi_override = None):

    try:        
        data = sarpy.Scan(scn_to_analyze).adata[adata_label].data
        roi = sarpy.Scan(scn_to_analyze).adata[roi_override].data

    except(KeyError,IOError):
        print(('adata_roi_average: Not found {0} and {1}'.format(scn_to_analyze,adata_label) ))
        data = numpy.empty([1])*numpy.nan
        weights= numpy.empty([1])*numpy.nan
        avgL = ['AnalysisErr']

    else:
        data_roi = data*roi
        #data[numpy.isnan(data)] = 0
        weights = sarpy.analysis.getters.get_roi_weights(roi)
        weights = list(weights)
        avg = []

        try:
            tumour_volume = sarpy.analysis.getters.get_tumour_volume(scn_to_analyze,roi_override)
        except AssertionError:
            tumour_volume = numpy.nan

        maxSlices = sarpy.analysis.getters.get_num_slices(scn_to_analyze)            

        if maxSlices>1:

            for slice in range(maxSlices): # Zero filling non-existent slices

                if slice < maxSlices:
                    avg.append(numpy.nanmean(data_roi[:,:,slice].flatten()))
                else:
                    avg.append(numpy.inf)
                    weights.append(0)
        else:
            avg.append(numpy.nanmean(data_roi[:,:].flatten()))

        # Removing infs and nans from the avg to get a proper weighted avg
        # Also removing 0 weights                    

        weights = numpy.array(weights)
        weights = weights[weights>0]

        avgW = numpy.array(avg)
        avgW = list(avgW[numpy.isfinite(avgW)])
        
        try:
            avg.append((numpy.average(avgW, weights=weights)))
        except TypeError:
            avg.append(numpy.nan)

        avg.append(tumour_volume)
        avgL=[str(e) for e in avg]

    return avgL

def create_export_csv(exp_abbreviation = 'HPGP4',
                      day_label = 'dce-0h',
                      data_scan_label = 'dce.HPG',
                      roi_scan_label = 'roi',
                      BATscreened_adata_label = 'BAT',
                      roi_adata = 'roi',
                      data_adata = 'auc60'):

    import time
    import getpass
    import os
    import csv

    export_data = []

    # First check to make sure that the number of slices for all patients is the same, 
    # if not, print an error message and record the exception in a dictionary for processing later
    
    allExperiment = sarpy.Experiment.from_masterlist(exp_abbreviation+'.config')
    
    maxSlicesList = []
    
    for scn in sorted(allExperiment.labels[data_scan_label]):
                
        maxSlicesList.append(sarpy.analysis.getters.get_num_slices(scn))

    maxSlices = numpy.max(maxSlicesList)
    
    for pat in sorted(allExperiment.patients.keys()):

        avg_data = []

        data_scn_to_analyze = allExperiment.patients[pat][data_scan_label]
        roi_scn_to_analyze = allExperiment.patients[pat][roi_scan_label]

        #data = sarpy.Scan(data_scn_to_analyze).adata[adata_label].data
        #roi = sarpy.Scan(roi_scn_to_analyze).adata[roi_adata].data   

        # Create empty list so that it can be populated with rows
        avg_data = determine_averages(data_scn_to_analyze,
                                      data_adata,
                                      BATscreened_adata_label,                                        
                                      roi_scn_to_analyze,
                                      roi_adata,
                                      maxSlices = maxSlices)
        
        avg_data.insert(0,data_scn_to_analyze)

        export_data.append(avg_data)

    timer = time.strftime("%Y-%m-%d %H:%M")   
    header =['Animal ID'] + ['Slice '+str(x+1) for x in range(
                                        len(max(export_data,key=len))-3)] + ['Weighted Average'] + ['Volume']

    footers = []
    footers.append(['# Data generated on {0} by {1}'.format(timer,getpass.getuser())])
    footers.append(["#Legend"])
    footers.append(["#inf:slice wasn't acquired"])
    footers.append(["#nan: no roi in slice"])
    footers.append(["#AnalysisErr: No Data available"])

    for f in footers:
        export_data.append(f)    


    filename = os.path.expanduser(os.path.join(
                                          '~','sdata',exp_abbreviation,'export',
                                          day_label+'_'+data_adata+'.csv'))

    # open output file
    outfile = open(filename, "w" )

    # get a csv writer
    writer = csv.writer( outfile )

    # write header
    writer.writerow(header)

    # write data
    [ writer.writerow(x) for x in export_data ]

    # close file
    outfile.close()    
       
def determine_averages(data_scn_to_analyze,
                       adata_label,
                       BAT_adata_label,
                       roi_scn_to_analyze,
                       roi_label,
                       maxSlices):
    
    data = sarpy.Scan(data_scn_to_analyze).adata[adata_label].data
    
    BAT = sarpy.Scan(data_scn_to_analyze).adata[BAT_adata_label].data
    roi = sarpy.Scan(roi_scn_to_analyze).adata[roi_label].data      
    
    data_roi = data*roi*BAT
    weights = sarpy.analysis.getters.get_roi_weights(roi)
    weights = list(weights)
    avg = []

    try:
        tumour_volume = sarpy.analysis.getters.get_tumour_volume(roi_scn_to_analyze,roi_label)
    except AssertionError:
        tumour_volume = numpy.nan

    totalSlices = sarpy.analysis.getters.get_num_slices(data_scn_to_analyze)            

    if totalSlices>1:

        for slc in range(maxSlices): # Zero filling non-existent slices

            if slc < totalSlices:
                avg.append(numpy.nanmean(data_roi[:,:,slc].flatten()))
            else:
                avg.append(numpy.inf)
                weights.append(0)
    else:
        avg.append(numpy.nanmean(data_roi[:,:].flatten()))
        
    # Removing infs and nans from the avg to get a proper weighted avg
    # Also removing 0 weights                    

    weights = numpy.array(weights)
    weights = weights[weights>0]

    avgW = numpy.array(avg)
    avgW = list(avgW[numpy.isfinite(avgW)])

    try:
        avg.append((numpy.average(avgW, weights=weights)))

    except TypeError:
        avg.append(numpy.nan)

    avg.append(tumour_volume)
    avgL=[str(e) for e in avg]

    return avg
           
if __name__ == '__main__':

    import argparse
    import os

    # we have been called as a script
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-m', '--masterlist_name', type=str, required=True,
                       help='Name of the masterlist file. e.g. NecS3 \
                       Usage: python write_csv.py -m HerP2 -d axref -a roi -i export ')
    
    parser.add_argument('-d', '--data_label', type=str, required=True,
                       help='Data label, usually anatomy or IR_A')

    parser.add_argument('-a', '--adata_label', type=str,required=True,
                       help='Data label, usually anatomy or IR_A')                          
                       
    args = parser.parse_args()
    masterlist_name = args.masterlist_name
    data_label = args.data_label
    adata_label = args.adata_label
    
    
    write_csv(masterlist_name,data_label,adata_label)

#######
# 
#######
