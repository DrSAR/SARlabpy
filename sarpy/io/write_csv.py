#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 13:53:03 2013

@author: fmoosvi
"""
def write_csv(masterlist, scan_label, adata_label):
    
    import json
    import os
    import collections
    import csv
    import scipy.stats
    import sarpy
    
    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist,
                        masterlist)
    fname_to_open = root+'.json'
    
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str) 
    
    export_data = []
    
    for k,v in master_list.iteritems():
        
        try:        
            data = sarpy.Scan(v[scan_label][0]).adata[adata_label].data
            avg = scipy.stats.nanmean(data)

        except KeyError:
            print('write_csv: Not found {0} and {1},{2}'.format(
                                                    k,scan_label,adata_label) )
            data = numpy.empty([1])*nan                                     
            avg = 'AnalysisErr'
            
            pass

        export_data.append([k] + data.tolist() + [avg])
#        print export_data
        				       
    header = ['Animal ID'] + ['Slice '+str(x+1) for x in xrange(
                                        len(export_data[0])-2)] + ['Average']

    filename = scan_label + '_' + adata_label + '.csv'
   
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
    
write_csv('HerS10', '0h-LL', 'T1avg')
    
# if __name__ == '__main__':    