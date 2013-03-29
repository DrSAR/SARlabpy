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


def create_summary(data_list, key_list):
     
    num_slices = data_list[0].shape[2]
    
    fig = pylab.figure(figsize = (30,30), dpi = 300)
    data_num = len(data_list)
    
    clims = sarpy.fmoosvi.getters.get_image_clims(data_list[0])
    
    for slice in xrange(num_slices):
        
        for n in xrange(data_num):
    
            fig.add_subplot(data_num,num_slices,slice+1 + (n*num_slices) )
            a = pylab.imshow(data_list[n][0,:,:,slice])
            a.set_clim(clims)
            pylab.axis('off')
            
            if n == 0:
                pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
    
    # Figure spacing adjustments
    fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
    
    # Colobar set up
    cax = fig.add_axes([0.89, 0.10, 0.03, 0.7])
    cax.set_title(key_list[len(key_list)-1], fontsize = 12)       
    fig.colorbar(a, cax=cax)
    
    # Saving Figure    
    filename = key_list[len(key_lis)-1] + '-' + data_list[0] + '.png'                
    pylab.savefig(filename, bbox_inches=0, dpi=300)
    pylab.close('all')        

#for k,v in master_sheet.iteritems():
#    try:           
#        data1 = sarpy.Scan(master_sheet[k]['0h-LL']).adata['T1map_LL'].data.get_data()
#        data2 = sarpy.Scan(master_sheet[k]['24h-LL']).adata['T1map_LL'].data.get_data()
#        
#        fig = pylab.figure(figsize=(14, 11), dpi=300)
#        
#        num_slices = sarpy.fmoosvi.getters.get_num_slices(sarpy.Scan(master_sheet[k]['0h-LL']))
#        
#        for slice in xrange(num_slices):
#            
#            fig.add_subplot(2,6,slice+1)
#            a = pylab.imshow(data1[0,:,:,slice])
#            a.set_clim(600,3500)
#            pylab.axis('off')
#            pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
#            fig.show()
#            
#            fig.add_subplot(2,6,slice+num_slices+1)
#            a = pylab.imshow(data2[0,:,:,slice])
#            a.set_clim(600,3500)
#            pylab.axis('off') 
#            #pylab.title('Slice {0}'.format(slice+1))
#            fig.show()
#
#        # Figure spacing adjustments
#        fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
#        
#        # Colobar set up
#        cax = fig.add_axes([0.89, 0.10, 0.03, 0.7])
#        cax.set_title('T$_1$ (ms)', fontsize = 12)       
#        fig.colorbar(a, cax=cax)
#    
#       # Saving Figure    
#        filename = 'T1-' + k + '.png'                
#        pylab.savefig(filename, bbox_inches=0, dpi=300)
#        pylab.close('all')        
# 
#    except:
#        print k, 'did not produce LL summary'
#    
#