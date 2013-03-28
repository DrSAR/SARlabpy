# -*- coding: utf-8 -*-

# Testing on T1 maps from NecS3

from sarpy.fmoosvi import getters
import sarpy
import pylab
import numpy

NecS3_patients = getters.get_patients_from_experiment('NecS3', check=True)

exceptions = ['NecS3Hs01a', 'NecS3Hs03']


## Saving IR images
for patient in NecS3_patients:
    
    IR = [study.find_scan_by_protocol('05')[0] for study in patient.studies if len(study.find_scan_by_protocol('06'))>0]
    
    
    if patient.patient_id in exceptions: continue
        
    try:     
        fig = pylab.figure(figsize=(8.5, 12), dpi=300)
    
        data1 = IR[0].pdata[0].data
        data2 = IR[1].pdata[0].data       
        
        num_slices = sarpy.fmoosvi.getters.get_num_slices(IR[0])
        
        for slice in xrange(num_slices):
            
            fig.add_subplot(2,6,slice+1)
            a = pylab.imshow(data1[:,:,slice])
            #a.set_clim(600,3500)
            pylab.axis('off')
            pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
            fig.show()
            
            fig.add_subplot(2,6,slice+num_slices+1)
            a = pylab.imshow(data2[:,:,slice])
            #a.set_clim(600,3500)
            pylab.axis('off') 
            fig.show()

        # Figure spacing adjustments
        fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
        
       # Saving Figure    
        filename = 'IR-' + patient.patient_id + '.png'                
        pylab.savefig(filename, bbox_inches=0, dpi=300)
        pylab.close('all')
                     
    except:
        print('You sir, {0}, are unhelpful'.format(patient.patient_id))
        

### Saving T1 maps
#for patient in NecS3_patients:
#    LL = patient.find_scan_by_protocol('04')
# 
#    try:     
#        fig = pylab.figure(figsize=(8.5, 12), dpi=300)
#    
#        data1 = LL[0].adata['T1map_LL'].data.get_data()
#        data2 = LL[1].adata['T1map_LL'].data.get_data()
#        
#        num_slices = sarpy.fmoosvi.getters.get_num_slices(LL[0])
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
#        filename = 'T1-' + patient.patient_id + '.png'                
#        pylab.savefig(filename, bbox_inches=0, dpi=300)
#        pylab.close('all')
#                     
#    except:
#        print('You sir, {0}, are unhelpful'.format(patient.patient_id))
#        
### Saving AUC60 figures
#
#for patient in NecS3_patients:
#    AUC60 = patient.find_scan_by_protocol('06')
# 
#    try:     
#        fig = pylab.figure(figsize=(8.5, 12), dpi=300)
#    
#        data1 = AUC60[0].adata['AUC60'].data.get_data()
#        #data2 = AUC60[2].adata['AUC60'].data.get_data()
#        
#        num_slices = sarpy.fmoosvi.getters.get_num_slices(AUC60[0])
#        
#        for slice in xrange(num_slices):
#            
#            fig.add_subplot(1,6,slice+1)
#            a = pylab.imshow(numpy.transpose(data1[0,:,:,slice], [1,0]))
#            a.set_clim(0,30)
#            pylab.axis('off')
#            pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
#            fig.show()
#            
##            fig.add_subplot(2,6,slice+num_slices+1)
##            a = pylab.imshow(numpy.transpose(data2[0,:,:,slice], [1,0,2]))
##            a.set_clim(0,60)
##            pylab.axis('off') 
##            #pylab.title('Slice {0}'.format(slice+1))
##            fig.show()
#
#        # Figure spacing adjustments
#        fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
#        
#        # Colobar set up
#        cax = fig.add_axes([0.89, 0.10, 0.03, 0.4])
#        cax.set_title('AUC60', fontsize = 12)       
#        fig.colorbar(a, cax=cax)
#    
#       # Saving Figure    
#        filename = 'AUC60' + patient.patient_id + '.png'                
#        pylab.savefig(filename, bbox_inches=0, dpi=300)
#        pylab.close('all')
#                     
#    except Exception, e: 
#        print('You sir, {0}, are unhelpful'.format(patient.patient_id))
#        print str(e)
#        pass