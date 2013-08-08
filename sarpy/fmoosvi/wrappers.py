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
import collections
import scipy
import scipy.stats

def bulk_analyze(masterlist_name, 
                 data_label, 
                 analysis_label, 
                 forceVal = False):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str) 


    if re.match('auc60', analysis_label):

        for k,v in master_list.iteritems():
            
            try:
                scan = sarpy.Scan(v[data_label][0])
                bbox = sarpy.fmoosvi.getters.get_bbox(v, data_label)

                if (not analysis_label in scan.adata.keys()) or forceVal is True:
                
                    curr_auc = sarpy.fmoosvi.analysis.h_calculate_AUC(scan, bbox)           
                    scan.store_adata(key=analysis_label, data = curr_auc, force = forceVal)
                
                else:
                    print('{0}: adata already exists {1}'.format(
                    analysis_label,scan.shortdirname))
                    pass 
                
            except(IOError):
                
                print('{0}: Not found {1} and {2}'.format(
                analysis_label,k,data_label) )                
                pass
            
            except ZeroDivisionError:
                pass            
            
                            
    elif re.match('T1map_LL', analysis_label):
        
        for k,v in master_list.iteritems():
            
            try: 
                scan = sarpy.Scan(v[data_label][0])           
                bbox = sarpy.fmoosvi.getters.get_bbox(v, data_label)
            
                if (not analysis_label in scan.adata.keys()) or forceVal is True:
                    T1map_LL, T1map_fitdict = sarpy.fmoosvi.analysis.h_fit_T1_LL(scan,bbox)
                    scan.store_adata(key=analysis_label, data = T1map_LL,force = forceVal)
                    scan.store_adata(key=analysis_label+'_fitdict', data = T1map_fitdict, force = forceVal)
                else:
                    print('{0}: adata already exists {1} '.format(
                    analysis_label,scan.shortdirname))
                    pass 
                
            except(IOError):
                print('{0}: Not found {1} and {2}'.format(
                analysis_label,k,data_label) )                
                pass

                
    elif re.match('vtc', analysis_label):
        
        for k,v in master_list.iteritems():
            
            try:
                scan = sarpy.Scan(v[data_label][0])           
                bbox = sarpy.fmoosvi.getters.get_bbox(v, data_label)

                if (not analysis_label in scan.adata.keys()) or forceVal is True:
                    vtc = sarpy.fmoosvi.analysis.h_generate_VTC(scan, bbox)
                    scan.store_adata(key=analysis_label, data = vtc, force = forceVal)
                else:
                    print('{0}: adata already exists {1}'.format(analysis_label,scan.shortdirname))
                    pass
                               
            except(IOError):
                print('{0}: Not found {1} and {2}'.format(
                analysis_label,k,data_label) )                
                pass

    elif re.match('roi_check',analysis_label):
        
        for k,v in master_list.iteritems():

#TODO: Change 'roi' to use roi_src_adata_label. Edit function            
            try:
                scan = sarpy.Scan(v[data_label][0])
                if (not analysis_label in scan.adata.keys()) or forceVal is True:
                    roi = scan.adata['roi'].data
                    masked_data = roi * scan.pdata[0].data   
                    scan.store_adata(key=analysis_label, data = masked_data, force = forceVal)
                else:
                    print('{0}: adata already exists {1}'.format(analysis_label,scan.shortdirname))
                    pass                                       
            except(IOError,KeyError):
                print('{0}: Not found {1} and {2}'.format(
                analysis_label,k,data_label) )                
                pass   
            

            
    else:
        print('This type of analysis has not yet been implemented. \
                Do so in the wrappers file.')
                
def calc_AUGC(masterlist_name,
              data_label,
              adata_scan_label,
              analysis_label='augc60',
              forceVal = False):
                  
    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)                    
       
    for k,v in master_list.iteritems():
        
        try:
            scan = sarpy.Scan(v[data_label][0])
            bbox = sarpy.fmoosvi.getters.get_bbox(v, data_label)
            
        except(IOError):
            print('{0}: Not found {1} and {2}'.format(
                                                analysis_label,k,data_label) )
            pass     
        
        if (not analysis_label in scan.adata.keys()) or forceVal is True:
        
            curr_augc = sarpy.fmoosvi.analysis.h_calculate_AUGC(scan, adata_scan_label, bbox=bbox)
            scan.store_adata(key=analysis_label, data = curr_augc, force = forceVal)
            print('{0}: Success. Saved {1}'.format(analysis_label, 
                              scan.shortdirname))
        
        else:
            print('{0}: adata already exists {1} '.format(
            analysis_label,scan.shortdirname))
            pass 
            
   

def conc_from_signal(masterlist_name, 
                     data_label, 
                     data_label_T1map, 
                     adata_label = 'T1map_LL', 
                     analysis_label='gd_conc', 
                     forceVal = False):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)                
        
    for k,v in master_list.iteritems():
        
        try:
            scan = sarpy.Scan(v[data_label][0])
            scan_T1map = sarpy.Scan(v[data_label_T1map][0])
            bbox = sarpy.fmoosvi.getters.get_bbox(v, data_label)
            
        except(IOError):          
            print('{0}: Not found {1} and {2}'.format(
            analysis_label,scan.shortdirname,data_label) )
            continue        
            
        if (not analysis_label in scan.adata.keys()) or forceVal is True:
            conc = sarpy.fmoosvi.analysis.h_conc_from_signal(scan, scan_T1map, adata_label, bbox)
            scan.store_adata(key=analysis_label, data = conc, force = forceVal)
            
            print('{0}: Success. Saved {1}'.format(analysis_label, 
                              scan.shortdirname))                
            
        else:
            print('{0}: adata already exists {1} '.format(
            analysis_label,scan.shortdirname))
            pass 
        


def store_deltaT1(masterlist_name,
                  T1map_1,
                  T1map_2,
                  adata_label1 = 'T1map_LL',
                  adata_label2 = None,
                  analysis_label='deltaT1',
                  forceVal = False):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

    if adata_label2 is None:
        adata_label2 = adata_label1

    for k,v in master_list.iteritems():

        try:
            scan1 = sarpy.Scan(v[T1map_1][0])
            scan2 = sarpy.Scan(v[T1map_2][0])
            bbox = sarpy.fmoosvi.getters.get_bbox(v, T1map_2)

    	except IOError:
    	    print('{0}: Not found adata {1} or {2} in patient {3}'.format(
    	          analysis_label,adata_label1,adata_label2,k) )
    	    continue

        if (not analysis_label in scan2.adata.keys()) or forceVal is True:
                
            deltaT1 = scan1.adata[adata_label1].data - scan2.adata[adata_label2].data
            scan2.store_adata(key=analysis_label, data = deltaT1, force = forceVal)

            print('{0}: Success. Saved {1}'.format(analysis_label,
                                  scan2.shortdirname))
                                  
        else:
            print('{0}: adata already exists {1} '.format(
            analysis_label,scan2.shortdirname))
            pass

### ROI based code

def roi_average(masterlist_name,
                data_label,
                adata_label,
                roi_label,
                analysis_label=None,
                forceVal = False):
                   
    if analysis_label is None:
        analysis_label = adata_label + '_avg'                   

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

    for k,v in master_list.iteritems():
        
        try:
            scan = sarpy.Scan(v[data_label][0])           
            data = scan.adata[adata_label].data
            roi = sarpy.Scan(v[data_label][0]).adata[roi_label].data
        except(KeyError, IOError):
            print('{0}: Not found adata {1} or {2} in patient {3}'.format(
                                    analysis_label,adata_label,roi_label,k) )
            continue
       
        roi_data = data * roi
        
        # I think this might be necessary to ensure that values don't get carrier
        # over from previous iterations?
        try:
            del avg_T1
        except NameError:
            pass
        avg_T1 = []     
       
        for slice in xrange(roi_data.shape[-1]):
            avg_T1.append(scipy.stats.nanmean(roi_data[:,:,slice].flatten()))
            
        if (not analysis_label in scan.adata.keys()) or forceVal is True:
                
            scan.store_adata(key=analysis_label, 
                             data = numpy.array(avg_T1), 
                             force = forceVal)

            print('{0}: Success. Saved {1}'.format(analysis_label,
                                  scan.shortdirname))
                                  
        else:
            print('{0}: adata already exists {1} '.format(
            analysis_label,scan.shortdirname))
            pass          

def bulk_transfer_roi(masterlist_name, dest_adata_label, roi_src_adata_label = None, tag = None, forceVal = False):
    '''
    Move an ROI from one scan to another. E.g., Moving an roi from an anatomy 
    scan to a LL scan. 
    
    dest_adata_label: Label specifying the destination of the roi transfer
    
    tag: string to prepend before ROI. useful if auc60_roi exists and you want 
            auc60_old_roi
    
    Note: this routine looks for the special 'roi' label. This is the canonical
    label used to specify ROIs and other scans will 'inherit' - by way of resample
    - this roi. 

    Example:
        
    sarpy.fmoosvi.wrappers.bulk_transfer_roi('NecS3',dest_adata_label,
                                             forceVal = False)
    '''
    
    # Get the name of the source roi
    if roi_src_adata_label is None:
        roi_src_adata_label = 'roi'

    # Set the name of the destination scan

    if tag is None:
        dest_label = dest_adata_label + '_roi'
    else:
        dest_label = dest_adata_label + str(tag) + '_roi'

#    # Get all the studies in the experiment

#    exp = sarpy.Experiment(exp_name)
    
    
    
    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)

    fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)

    for k in master_list.iteritems():
        
        for study in sarpy.Patient(k[0]).studies:
            adata_dict = study.find_adata_scans()
        
            try:
                src = sarpy.Scan(adata_dict[roi_src_adata_label][0])
    
            except(IOError):
                print('bulk_transfer_roi: IO Problem with roi src scan for study {0}'.format(study.shortdirname))
                raise
            except KeyError:
                print('bulk_transfer_roi: Key roblem with roi src scan for study {0}, maybe this study is extra'.format(study.shortdirname))


    
            if dest_adata_label not in adata_dict:
                # Skip this Patient
                print('bulk_transfer_roi: {0} not found in study {1}'.format(dest_adata_label,study.shortdirname))                
                continue
                
                
            for d in adata_dict[dest_adata_label]:
                dest = sarpy.Scan(d)
           
                if (not dest_label in dest.adata.keys()) or forceVal is True:
                    
                    dest_pd = dest.pdata[0]
                    roi = src.adata[roi_src_adata_label]
                    
                    resampled_roi = sarpy.ImageProcessing.resample_onto.\
                        resample_onto_pdata(roi, dest_pd, replace_nan=0)
                        
                    places = numpy.where(resampled_roi < .5)
                    other_places = numpy.where(resampled_roi >= .5)
                    resampled_roi[places] = numpy.nan
                    resampled_roi[other_places] = 1
    
                    dest.store_adata(key = dest_label, data = resampled_roi, 
                                     force = forceVal)   
                                     
                    print('bulk_transfer_roi: Success. Saved {0} in {1}'.format(dest_label, 
                                  dest.shortdirname))                                 
                else:
                    print('bulk_transfer_roi: {0} adata already exists {1}'.format(dest_label, 
                          dest.shortdirname))
                    pass      

    #TODO: Fix this function so it uses the masterlist    
    print('Firas: you really need to fix this function so it uses the masterlist')

def calc_enhancement_curve(masterlist_name, 
                           data_label,
                           analysis_label = 'ec', 
                           roi_label = 'auc60_roi',
                           forceVal = False):


    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str) 
        
    for k,v in master_list.iteritems():
        
        try:
            scan = sarpy.Scan(v[data_label][0])
            
            try: # check if the roi_label is good or not
                roi = scan.adata[roi_label].data
                 
                if (not analysis_label in scan.adata.keys()) or forceVal is True:
            
                    conc = sarpy.fmoosvi.analysis.h_enhancement_curve(scan, roi_label)
                    scan.store_adata(key=analysis_label, data = conc, force = forceVal)

                else:
                    print('{0} adata already exists {1}'.format(analysis_label,scan.shortdirname))
                    pass 
                                                 
            except KeyError:
                print('enhancement_curve: {0} adata for scan {1} does not exist'.format(roi_label,scan.shortdirname) )

        except IOError:
            
            print('{0}: Not found {1} and {2}'.format(k,data_label,analysis_label) )
            
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

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)
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
    
            
        
        
        
        
        
        
        
        

    
