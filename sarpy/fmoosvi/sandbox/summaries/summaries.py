# -*- coding: utf-8 -*-

# Testing on T1 maps from NecS3

from sarpy.fmoosvi import getters
import sarpy
import pylab

NecS3_patients = getters.get_patients_from_experiment('NecS3', check=True)

for patient in NecS3_patients:

    LL = patient.find_scan_by_protocol('04')
    
        fig = pylab.figure()

    for session in LL:
        try:
            
            ax1 = fig.add_subplt(1,2,1)
            ax1.imshow(LL[0].adata['T1map_LL'].data.get_data()[:,:,0])
            
            ax2 = fig.add_subplt(1,2,2)
            ax1.imshow(LL[1].adata['T1map_LL'].data.get_data()[:,:,1])
        except:
           # print('You sir, are unhelpful')