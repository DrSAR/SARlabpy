#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 13:53:03 2013

@author: fmoosvi
"""
def write_csv(masterlist_name, data_label, adata_label):
    
    import json
    import os
    import collections
    import csv
    import scipy.stats
    import sarpy
    import time
    import getpass # used to get the current username
    import numpy
    
    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    fname_to_open = root+'.json'
    
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str) 
    
    export_data = []
    
    for k,v in master_list.iteritems():
        
        try:        
            data = sarpy.Scan(v[data_label][0]).adata[adata_label].data
            weights = sarpy.Scan(v[data_label][0]).adata[adata_label+'_weights'].data

            data[numpy.isnan(data)] = 0
            weights[numpy.isnan(weights)] = 0
            avg = numpy.average(data, weights=weights)

        except(KeyError,IOError):
            print('write_csv: Not found {0} and {1},{2}'.format(
                                                    k,data_label,adata_label) )
            data = numpy.empty([1])*numpy.nan                                     
            avg = 'AnalysisErr'
            
            pass

        export_data.append([k] + data.tolist() + [avg])

    time = time.strftime("%Y-%m-%d %H:%M")       				       
    header =['Animal ID'] + ['Slice '+str(x+1) for x in xrange(
                                        len(max(export_data,key=len))-2)] + ['Average']

    footer =  ['# Data generated on {0} by {1}'.format(time,getpass.getuser())]
    
    export_data.append(footer)    
    

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
    
    