#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 13:53:03 2013

@author: fmoosvi
"""
def write_csv(masterlist_name, 
              data_label, 
              adata_label,
              roi_override = None):
    
    import json
    import os
    import collections
    import csv
    import scipy.stats
    import sarpy
    import time
    import getpass # used to get the current username
    import numpy
    import sarpy.fmoosvi.getters
    
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
    
    export_data = []

    # Get Max number of slices, this is so that the average column is always the last column
    # and so slice averages get zero-filled

    numSlices = []
    for k,v in master_list.iteritems():

        try:        
            data = sarpy.Scan(v[data_label][0]).adata[adata_label].data
            roi = sarpy.Scan(v[data_label][0]).adata[roi_override].data  

            numSlices.append(data.shape[-1])
        except(IOError,KeyError):
            continue     

    if numSlices:
        maxSlices = numpy.max(numSlices)
    
    for k,v in master_list.iteritems():
        
        try:        
            data = sarpy.Scan(v[data_label][0]).adata[adata_label].data
            roi = sarpy.Scan(v[data_label][0]).adata[roi_override].data

        except(KeyError,IOError):
            print('write_csv: Not found {0} and {1},{2}'.format(k,data_label,adata_label) )
            data = numpy.empty([1])*numpy.nan
            weights= numpy.empty([1])*numpy.nan
            avgL = ['AnalysisErr']

        else:
            data_roi = data*roi
            #data[numpy.isnan(data)] = 0
            weights = sarpy.fmoosvi.getters.get_roi_weights(roi)
            weights = weights.tolist()
            avg = []

            tumour_volume = sarpy.fmoosvi.getters.get_tumour_volume(v[data_label][0],roi_override)

            for slice in xrange(maxSlices): # Zero filling non-existent slices

                if slice < data_roi.shape[-1]:
                    avg.append(scipy.stats.nanmean(data_roi[:,:,slice].flatten()))
                else:
                    avg.append(numpy.inf)
                    weights.append(0)

            # Removing infs and nans from the avg to get a proper weighted avg
            # Also removing 0 weights                    

            weights = numpy.array(weights)
            weights = weights[weights>0]

            avgW = numpy.array(avg)
            avgW = avgW[numpy.isfinite(avgW)].tolist()
            
            avg.append((numpy.average(avgW, weights=weights)))
            avg.append(tumour_volume)
            avgL=[str(e) for e in avg]

        export_data.append([k] + avgL)

    time = time.strftime("%Y-%m-%d %H:%M")       				       
    header =['Animal ID'] + ['Slice '+str(x+1) for x in xrange(
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
    outfile = open(filename, "wb" )
    
    # get a csv writer
    writer = csv.writer( outfile )
    
    # write header
    writer.writerow(header)
    
    # write data
    [ writer.writerow(x) for x in export_data ]
    
    # close file
    outfile.close()
       
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
    
    
