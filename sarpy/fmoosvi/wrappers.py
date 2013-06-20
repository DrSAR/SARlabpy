# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:11:03 2013

@author: fmoosvi
"""

import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters
import pylab
import numpy
import os
import json
import re
import sarpy.ImageProcessing.resample_onto

def bulk_analyze(masterlist_name, data_label, analysis_label, forceVal = False):
    
    mdata = os.path.expanduser(os.path.join('~','mdata',masterlist_name+'.json'))
    
    with open(mdata,'r') as master_file:
        master_list = json.load(master_file)
        
    if re.match('auc60', analysis_label):

        for k,v in master_list.iteritems():
            
            try:
                scan = sarpy.Scan(v[data_label][0])
                bbox = sarpy.fmoosvi.getters.get_bbox(v, data_label)
                
                if (not analysis_label in [analysis_label]) or forceVal is True:
                
                    curr_auc = sarpy.fmoosvi.analysis.h_calculate_AUC(scan, bbox)
                    scan.store_adata(key=analysis_label, data = curr_auc, force = forceVal)
                
                else:
                    print('adata already exists {0}'.format(scan.shortdirname))
                    pass 
                
            except IOError:
                
                print('Not found: {0} and {1} and {2}'.format(k,data_label,analysis_label) )
                
                pass
        
    elif re.match('T1map_LL', analysis_label):
        
        for k,v in master_list.iteritems():
            
            try: 
                scan = sarpy.Scan(v[data_label][0])           
                bbox = sarpy.fmoosvi.getters.get_bbox(v, data_label)
            
                if (not analysis_label in [analysis_label]) or forceVal is True:
                    T1map_LL, T1map_fitdict = sarpy.fmoosvi.analysis.h_fit_T1_LL(scan,bbox)
                    scan.store_adata(key=analysis_label, data = T1map_LL,force = forceVal)
                    scan.store_adata(key=analysis_label+'_fitdict', data = T1map_fitdict, force = forceVal)
                else:
                    print('adata already exists {0}'.format(scan.shortdirname))
                    pass
                
            except IOError:
                print('Not found: {0} and {1} and {2}'.format(k,data_label,analysis_label) )
                pass

                
    elif re.match('vtc', analysis_label):
        
        for k,v in master_list.iteritems():
            
            try:
                scan = sarpy.Scan(v[data_label][0])           
                bbox = sarpy.fmoosvi.getters.get_bbox(v, data_label)

                if (not analysis_label in [analysis_label]) or forceVal is True:
                    vtc = sarpy.fmoosvi.analysis.h_generate_VTC(scan, bbox)
                    scan.store_adata(key=analysis_label, data = vtc, force = forceVal)
                else:
                    print('adata already exists {0}'.format(scan.shortdirname))
                    pass
                               
            except IOError:
                
                print('Not found: {0} and {1} and {2}'.format(k,data_label,analysis_label) )
                pass

    elif re.match('roi_check',analysis_label):
        
        for k,v in master_list.iteritems():
            
            try:
                
                if (not analysis_label in [analysis_label]) or forceVal is True:
                    scan = sarpy.Scan(v[data_label][0])
    
                    roi = scan.adata['roi'].data
    
                    masked_data = roi * scan.pdata[0].data
                    
                    scan.store_adata(key=analysis_label, data = masked_data, force = forceVal)
                else:
                    print('adata already exists {0}'.format(scan.shortdirname))
                    pass
                                               
            except IOError:
                
                print('Not found: {0} and {1} and {2}'.format(k,data_label,analysis_label) )
                pass            
            
    else:
        print('This type of analysis has not yet been implemented. \
                Do so in the wrappers file.')     

def bulk_transfer_roi(masterlist_name, src_data_label, dest_data_label, 
                 analysis_label, forceVal = False):

    mdata = os.path.expanduser(os.path.join('~','mdata',masterlist_name+'.json'))
    
    with open(mdata,'r') as master_file:
        master_list = json.load(master_file)
          
    for k,v in master_list.iteritems():
        
        try:
            src_scan = sarpy.Scan(v[src_data_label][0])
            dest_scan = sarpy.Scan(v[dest_data_label][0])
            
            curr_adata = dest_scan.adata.keys()
                            
            if (not analysis_label in curr_adata) or forceVal is True:

                pd = dest_scan.pdata[0]
                roi = src_scan.adata['roi']

                resampled_roi = sarpy.ImageProcessing.resample_onto.\
                resample_onto_pdata(roi,pd)
                                                 
                dest_scan.store_adata(key=analysis_label, data = resampled_roi,
                                 force = forceVal)
            else:
                print('adata already exists {0}'.format(dest_scan.shortdirname))
                pass                
        except IOError:
            
            print('Not found: {0} and {1} and {2}'.format(k,data_label,analysis_label) )
            pass

def roi_distribution(data, roi, bins,  display_histogram = True, 
                     save_histogram = False, save_name = 'hist'):
   
    if type(roi) == numpy.ndarray and data == numpy.ndarray :
        masked_data = data * roi
    elif type(roi) == numpy.ndarray :
        masked_data = data.data * roi
    elif type(roi) == sarpy.io.AData_classes.AData:
        roi_mask = sarpy.fmoosvi.analysis.h_image_to_mask(roi.data)
        masked_data = data.data * roi_mask
      
        
    if display_histogram:
        #fig = pylab.figure()
        pylab.hist(masked_data.flatten(), bins, alpha = 0.5)
        pylab.title('Distribution of T1s in the ROI')
        
    if save_histogram:
        
        filename = save_name + '.png'                
        pylab.savefig(filename, bbox_inches=0, dpi=300)
          
def plotVTC(masterlist_name, key, data_label, adata_label = None):

    mdata = os.path.expanduser(os.path.join('~','mdata',masterlist_name + '.json'))
    
    with open(mdata,'r') as master_file:
        master_list = json.load(master_file)
        
    value = master_list[key]
        
    data = sarpy.Scan(value[data_label][0])
    bbox_px = sarpy.fmoosvi.getters.get_bbox(value,data_label)
    
    if adata_label is None:
        nrdata = sarpy.Scan(value[data_label][0]).adata['vtc'].data        
    else:      
        nrdata = sarpy.Scan(value[data_label][0]).adata[adata_label].data
        
    num_slices = nrdata.shape[2]
 
    x_bbox_px = bbox_px[1] - bbox_px[0] +1

    # Handle the special case of a corner bbox at origin
    if bbox_px[0]==0: x_bbox_px - 1
    
#    pylab.figure(figsize = [10,data.method.PVM_FovCm[1]/data.method.PVM_FovCm[0]*10])
#    pylab.imshow(data.pdata[0].data[:,:,3,80])
#   
    for slice in xrange(num_slices):
        G = pylab.matplotlib.gridspec.GridSpec(x_bbox_px, 1)
        fig = pylab.figure(figsize = [10,data.method.PVM_FovCm[1]/data.method.PVM_FovCm[0]*10])       
    
        i = 0
        
        for x in xrange(bbox_px[0],bbox_px[1]):
         
            fig.add_subplot(G[i],frameon=False,xticks=[],yticks=[])
        
            d = nrdata[x,:,slice]
        
            pylab.plot(d,'-',linewidth=0.1)
            pylab.ylim([-0.2,2])
            
            i+=1 # increment counter for the subplot (x)             
        
        
        
        
### Retired Code

def create_summary(data_list, key_list, clims = None, 
                   colour_map = 'jet'):
     
    num_slices = data_list[0].shape[-1]
    data_num = len(data_list)
    
    fig = pylab.figure(figsize = (14,3))
    G = pylab.matplotlib.gridspec.GridSpec(data_num,num_slices)   
    
    for slice in xrange(num_slices):
        
        for n in xrange(data_num):    
            fig.add_subplot(G[n,slice],frameon=False, xticks=[], yticks=[])
            
            try:
                data = data_list[n][0,:,:,slice]  
            except:
                data = data_list[n][:,:,slice]
            
            a = pylab.imshow(data, cmap = colour_map )
            
            if isinstance(clims,list):                
                a.set_clim(clims)
            else:
                clims = sarpy.fmoosvi.getters.get_image_clims(data)
                a.set_clim(clims[0],clims[1])
            if n == 0:
                pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
    

    # Colobar set up
    cax = fig.add_axes([0.9, 0.10, 0.01, 0.7])
    cax.set_title(key_list[1], fontsize = 12)       
    fig.colorbar(a, cax=cax)
    
    # Figure spacing adjustments
    #fig.subplots_adjust(right = 0.85, wspace = 0.0001, hspace=0.0001)
    G.tight_layout(fig, h_pad = 0.1, w_pad = 0.001)
    G.update(right = 0.87)
    
    # Saving Figure    
    filename = key_list[1] + key_list[0] + '.png'                
    pylab.savefig(filename, bbox_inches=0, dpi=300)
    pylab.close('all')        
        
def calculate_BSB1map(Bruker_object, BS_protocol_name = '07_bSB1mapFLASH',
                      POI_protocol_name = '04_ubcLL+'):
    
    # Can take in an experiment, a list of scans (4 - BS + 1 scan with poi),
    # a study
    
    print ("Why are you using this!? It's not fully implemented yet")
    
    try:
        if type(Bruker_object) == sarpy.io.BRUKER_classes.Experiment \
          or type(Bruker_object) == sarpy.io.BRUKER_classes.Study:
              
            scan_list = Bruker_object.find_scan_by_protocol(BS_protocol_name)
            LL_scan_list = Bruker_object.find_scan_by_protocol(POI_protocol_name)
        else:
            scan_list = []                                
            for study in Bruker_object:
                scan_list.append(study.find_scan(BS_protocol_name))
                LL_scan_list = Bruker_object.find_scan_by_protocol(POI_protocol_name)
            print('You put in a list of studies, try to avoid that')
    except:
        raise

    ## Now to calculate the B1 map
    b1map = []
    
    for scan in LL_scan_list:
        zero_BSminus = scan_list[0]
        zero_BSplus = scan_list[1]
        high_BSminus = scan_list[2]
        high_BSplus = scan_list[3]
        scan_with_POI = LL_scan_list
    
    
    b1map.append(sarpy.fmoosvi.analysis.h_BS_B1map(zero_BSminus, zero_BSplus, \
                                                   high_BSminus, high_BSplus, \
                                                   scan_with_POI))

    # Get rid of annoying extra dimension if scan_list contains only one element
    # Also return arrays instead of a list
            
    if numpy.array(b1map).shape[0] == 1:     
        return numpy.array(b1map)[0,:,:,:]
    else:
        return numpy.array(b1map)                                                

def create_plot(data_list, key_list):
        
    num_slices = data_list[0].shape[0]
    data_num = len(data_list)
    
    fig = pylab.figure(figsize = (14,3))
    G = pylab.matplotlib.gridspec.GridSpec(data_num,num_slices)   
    
    for slice in xrange(num_slices):
        
        for n in xrange(data_num):    
        
#            fig.add_subplot(G[n,slice],frameon=False, xticks=[], yticks=[])
            fig.add_subplot(G[n,slice])
        
            try:
                data = data_list[n][slice]  
            except:
                raise
        
        pylab.plot(data)
        
        pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
      
    # Figure spacing adjustments
    #fig.subplots_adjust(right = 0.85, wspace = 0.0001, hspace=0.0001)
    G.tight_layout(fig, h_pad = 0.3, w_pad = 0.3)
    #G.update(right = 0.87)
    
    # Saving Figure    
    filename = key_list[1] + key_list[0] + '.png'                
    pylab.savefig(filename, bbox_inches=0, dpi=300)
    pylab.close('all')
    
            
        
        
        
        
        
        
        
        

    