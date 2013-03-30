# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:11:03 2013

@author: fmoosvi
"""

import sarpy
import sarpy.fmoosvi.analysis
import sarpy.fmoosvi.getters
import pylab

def calculate_AUC(Bruker_object, protocol_name = '06_FLASH2D+', pdata_num = 0, time = 60):


    if type(Bruker_object) == sarpy.io.BRUKER_classes.Scan:
        
        scan_list = []
        scan_list.append(Bruker_object)
        
    else:
      
        try: # Will work for experiment and study
            scan_list = Bruker_object.find_scan_by_protocol(protocol_name)
        except: # check for list of studies
            try: 
                scan_list = []     
    
                if type(Bruker_object) == list:
                    scan_list = Bruker_object
                    print('You put in a list of scans, try to avoid that')
    
                elif type(Bruker_object[0]) == sarpy.io.BRUKER_classes.Study:
                    for study in Bruker_object:
                        scan_list.append(study.find_scan(protocol_name))
                        print('You put in a list of studies, try to avoid that')              
            except:
                raise

    ## Now to calculate the AUC

    auc = []

    for scan in scan_list:
        try:
            auc.append(sarpy.fmoosvi.analysis.h_calculate_AUC(scan))
        except:
            print('calculate_AUC failed for Scan {0} failed, please fix.'.format(scan.shortdirname))
        
    return auc
    
def calculate_BSB1map(Bruker_object, BS_protocol_name = '07_bSB1mapFLASH', \
                      POI_protocol_name = '04_ubcLL+'):
    
    # Can take in an experiment, a list of scans (4 - BS + 1 scan with poi),
    # a study
    
    print ("Why are you using this!? It's not implemented yet")
    
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
        
    return b1map

    return #b1map
    
def calculate_T1map(Bruker_object, protocol_name = '04_ubcLL2', flip_angle_map = 0):

    try:       
        if type(Bruker_object) == sarpy.io.BRUKER_classes.Scan:
            scan_list = []
            scan_list.append(Bruker_object)
        elif type(Bruker_object) == sarpy.io.BRUKER_classes.Experiment \
          or type(Bruker_object) == sarpy.io.BRUKER_classes.Study:
              
            scan_list = Bruker_object.find_scan_by_protocol(protocol_name)
        else:
            scan_list = Bruker_object
            
    except:
        try: 
            scan_list = []                                
            for study in Bruker_object:
                scan_list.append(study.find_scan_by_protocol(protocol_name))
            print('You put in a list of studies, try to avoid that')
        except:
            raise

    ## Now to calculate the AUC
    T1_map = []

    for scan in scan_list:
        
        print scan
        T1_map.append(sarpy.fmoosvi.analysis.h_fit_T1_LL(Bruker_object,\
                                                      flip_angle_map = flip_angle_map))
    return T1_map


def create_summary(data_list, key_list, clims = None, colour_map = 'jet'):
     
    num_slices = data_list[0].shape[-1]
    data_num = len(data_list)
    
    fig = pylab.figure(figsize = (14,3))
    G = pylab.matplotlib.gridspec.GridSpec(data_num,num_slices)   
    #clims = sarpy.fmoosvi.getters.get_image_clims(data_list[0])
    
    for slice in xrange(num_slices):
        
        for n in xrange(data_num):
    
            fig.add_subplot(G[n,slice],frameon=False, xticks=[], yticks=[])
            
            try:
                data = data_list[n][0,:,:,slice]  
            except:
                data = data_list[n][:,:,slice] 
            
            a = pylab.imshow(data, cmap = colour_map )
            
            if isinstance(clims[0],(int, float, long)):
                
                a.set_clim(clims)
                
            else:
                a.set_clim(600,3500)

            #a.set_clim(clims)
            #pylab.axis('off')

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





























