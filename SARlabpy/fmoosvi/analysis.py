# Copyright (C) 2012-2013 Stefan A Reinsberg and SARlab members
# full license details see LICENSE.txt
"""Collection of analysis routines for Bruker data


"""

import numpy
import os
import SARlabpy as sar

def calculateAUC(study_path, series_num = 0):

    # Read in data and headers

    procdirname = os.path.join(study_path,series_num)


    if os.path.exists(procdirname):
        
        data_dict = sar.read2dseq(procdirname,typecast=True)
        
        repetition_time = float(data_dict['header']['method']['PVM_RepetitionTime'])*1E-3
        num_slices = int(data_dict['header']['method']['PVM_SPackArrNSlices'])
        reps = int(data_dict['header']['method']['PVM_NRepetitions'])
        
        # Determine point of injection by averaging
        
        
    else:
        
        print "This directory doesn't exist, please enter a valid dir"




def enhancementCurve(data_dict, mask=False):
        
    if mask:
        
        print "I don't know how to deal with masks quite yet"
        
    else:
        
        try:
            img_mean = data_dict['data'][:,:,:,:].mean(1).mean(1)
            
            return img_mean
            
        except IndexError:
            print "You might only have 3D or 3D data, check data source"