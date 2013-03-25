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

            data1 = LL[0].adata['T1map_LL'].data.get_data()
            data2 = LL[1].adata['T1map_LL'].data.get_data()
            
            ax1 = fig.add_subplot(1,2,1)
            a = pylab.imshow(data1[0,:,:,0])
            a.set_clim(600,3200)
            pylab.colorbar()
            fig.show()
            
            ax2 = fig.add_subplot(1,2,2)
            b = pylab.imshow(data2[0,:,:,1])
            b.set_clim(600,3200)
            pylab.colorbar()
            fig.show()

        except:
            print('You sir, are unhelpful')
            
#            
#pylab.figure()
#img = pylab.imshow(proc_data[0][:,:,0], aspect = 0.5)
#img.set_clim(0.0,1200)
#pylab.colorbar()
#pylab.show()