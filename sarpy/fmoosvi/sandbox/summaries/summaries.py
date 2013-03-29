# -*- coding: utf-8 -*-

# Testing on T1 maps from NecS3

import sarpy
import sarpy.fmoosvi.getters
import sarpy.fmoosvi.wrappers
import pylab
import numpy
import json

with open('/Volumes/Data/Dropboxes/PhD./Dropbox/Studies/NecS3/NecS3.json','r') as infile:
    master_sheet = json.load(infile)
    
###################### LL T1 Maps ##############################################  



#data_list = []
#for k,v in master_sheet.iteritems():
#
#    try:
#        key_list = []
#        key_list.append(k)
#        key_list.append('0h-LL')
#        key_list.append('24h-LL')
#        key_list.append('T1map_LL')
#        key_list.append('LL T1 Map')
#        
#        data1 = sarpy.Scan(master_sheet[k][key_list[1]]).adata[key_list[3]].data.get_data()
#        data2 = sarpy.Scan(master_sheet[k][key_list[2]]).adata[key_list[3]].data.get_data()
#        
#        data_list.append(data1)
#        data_list.append(data2)
#        
#        sarpy.fmoosvi.wrappers.create_summary(data_list, key_list)
#        
#    except:
#        print k
#        raise
    
    
###################### LL T1 Maps ##############################################     
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
###################### AUC60 Maps ##############################################       
#for k,v in master_sheet.iteritems():
#    
#    try: ## AUC60 maps taken twice
#        try:
#            data1 = sarpy.Scan(master_sheet[k]['0h-DCE1'][0]).adata['AUC60'].data.get_data()
#            data2 = sarpy.Scan(master_sheet[k]['24h-DCE2'][0]).adata['AUC60'].data.get_data()
#    
#            fig = pylab.figure(figsize=(14, 11), dpi=300)
#            
#            num_slices = sarpy.fmoosvi.getters.get_num_slices(sarpy.Scan(master_sheet[k]['0h-DCE1']))
#            
#            for slice in xrange(num_slices):
#                
#                fig.add_subplot(2,6,slice+1)
#                a = pylab.imshow(data1[0,:,:,slice])
#                a.set_clim(10,50)
#                pylab.axis('off')
#                pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
#                fig.show()
#                
#                fig.add_subplot(2,6,slice+num_slices+1)
#                a = pylab.imshow(data2[0,:,:,slice])
#                a.set_clim(10,50)
#                pylab.axis('off') 
#                #pylab.title('Slice {0}'.format(slice+1))
#                fig.show()
#    
#            # Figure spacing adjustments
#            fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
#            
#            # Colobar set up
#            cax = fig.add_axes([0.89, 0.10, 0.03, 0.7])
#            cax.set_title('AUC60', fontsize = 12)       
#            fig.colorbar(a, cax=cax)
#        
#           # Saving Figure    
#            filename = 'AUC60-' + k + '.png'
#            pylab.savefig(filename, bbox_inches=0, dpi=300)
#            pylab.close('all')
#            
#        except:
#            pass
#    
#        try: ## AUC60 maps taken once
#            data1 = sarpy.Scan(master_sheet[k]['24h-DCE1'][0]).adata['AUC60'].data.get_data()
#    
#            fig = pylab.figure(figsize=(14, 11), dpi=300)
#            
#            num_slices = sarpy.fmoosvi.getters.get_num_slices(sarpy.Scan(master_sheet[k]['24h-DCE1'][0]))
#            
#            for slice in xrange(num_slices):
#                
#                fig.add_subplot(2,6,slice+num_slices+1)
#                a = pylab.imshow(data1[0,:,:,slice])
#                a.set_clim(10,50)
#                pylab.axis('off') 
#                #pylab.title('Slice {0}'.format(slice+1))
#                fig.show()
#    
#            # Figure spacing adjustments
#            fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
#            
#            # Colobar set up
#            cax = fig.add_axes([0.89, 0.10, 0.03, 0.7])
#            cax.set_title('AUC60', fontsize = 12)       
#            fig.colorbar(a, cax=cax)
#        
#           # Saving Figure    
#            filename = 'AUC60-' + k + '.png'                
#            pylab.savefig(filename, bbox_inches=0, dpi=300)
#            pylab.close('all')     
#        except:
#            pass
#        
#    except:
#        print k, 'did not produce AUC60 summary'
#
##################### IR-RARE ##############################################       
for k,v in master_sheet.iteritems():
    
    try: ## Two IR RARE on Day1
        try:
            data1 = sarpy.Scan(master_sheet[k]['0h-IR_A'][0]).pdata[0].data
            data2 = sarpy.Scan(master_sheet[k]['0h-IR_B'][0]).pdata[0].data
            data12 = numpy.vfig = pylab.figurestack([data1,data2])
    
            data3 = sarpy.Scan(master_sheet[k]['24h-IR_A'][0]).pdata[0].data
            data4 = sarpy.Scan(master_sheet[k]['24h-IR_B'][0]).pdata[0].data
            data34 = numpy.vstack([data3,data4])
 
            fig = pylab.figure(figsize=(14, 11), dpi=300)
            
            num_slices = sarpy.fmoosvi.getters.get_num_slices(sarpy.Scan(master_sheet[k]['0h-IR_A'][0]))
            
            for slice in xrange(num_slices):
                
                fig.add_subplot(2,6,slice+1)
                a = pylab.imshow(data12[:,:,slice], cmap='gray')
                #a.set_clim(10,50)
                pylab.axis('off')
                pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
                fig.show()
                
                fig.add_subplot(2,6,slice+num_slices+1)
                a = pylab.imshow(data34[:,:,slice], cmap='gray')
                #a.set_clim(10,50)
                pylab.axis('off') 
                #pylab.title('Slice {0}'.format(slice+1))
                fig.show()
    
            # Figure spacing adjustments
            fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
                    
           # Saving Figure    
            filename = 'IR-' + k + '.png'
            pylab.savefig(filename, bbox_inches=0, dpi=300)
            pylab.close('all')
            
        except:
            pass
    
        try:## One IR RARE on Day1
            data1 = sarpy.Scan(master_sheet[k]['0h-IR_A'][0]).pdata[0].data
            data2 = numpy.zeros(shape = data1.shape)
            data12 = numpy.vstack([data1,data2])
    
            data3 = sarpy.Scan(master_sheet[k]['24h-IR_A'][0]).pdata[0].data
            data4 = sarpy.Scan(master_sheet[k]['24h-IR_B'][0]).pdata[0].data
            data34 = numpy.vstack([data3,data4])
            
            fig = pylab.figure(figsize=(14, 11), dpi=300)
            
            num_slices = sarpy.fmoosvi.getters.get_num_slices(sarpy.Scan(master_sheet[k]['0h-IR_A'][0]))
            
            for slice in xrange(num_slices):
                
                fig.add_subplot(2,6,slice+1)
                a = pylab.imshow(data12[:,:,slice], cmap='gray')
                #a.set_clim(10,50)
                pylab.axis('off')
                pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
                fig.show()
                
                fig.add_subplot(2,6,slice+num_slices+1)
                a = pylab.imshow(data34[:,:,slice], cmap='gray')
                #a.set_clim(10,50)
                pylab.axis('off') 
                #pylab.title('Slice {0}'.format(slice+1))
                fig.show()
    
            # Figure spacing adjustments
            fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
            
            # Saving Figure    
            filename = 'IR-' + k + '.png'
            pylab.savefig(filename, bbox_inches=0, dpi=300)
            pylab.close('all')
            
        except:
            pass
    except:
        print k, 'did not produce IR-RARE summary'


### Saving IR images
#for patient in NecS3_patients:
#    
#    IR = [study.find_scan_by_protocol('05')[0] for study in patient.studies if len(study.find_scan_by_protocol('05'))>0]
#    
#    
#    if patient.patient_id in exceptions: continue
#        
#    try:     
#        fig = pylab.figure(figsize=(8.5, 12), dpi=300)
#    
#        data1 = IR[0].pdata[0].data
#        data2 = IR[1].pdata[0].data       
#        
#        num_slices = sarpy.fmoosvi.getters.get_num_slices(IR[0])
#        
#        for slice in xrange(num_slices):
#            
#            fig.add_subplot(2,6,slice+1)
#            a = pylab.imshow(data1[:,:,slice])
#            #a.set_clim(600,3500)
#            pylab.axis('off')
#            pylab.title('Slice {0}'.format(slice+1), fontsize = 14)
#            fig.show()
#            
#            fig.add_subplot(2,6,slice+num_slices+1)
#            a = pylab.imshow(data2[:,:,slice])
#            #a.set_clim(600,3500)
#            pylab.axis('off') 
#            fig.show()
#
#        # Figure spacing adjustments
#        fig.subplots_adjust(right = 0.85, wspace = 0, hspace=0)
#        
#       # Saving Figure    
#        filename = 'IR-' + patient.patient_id + '.png'                
#        pylab.savefig(filename, bbox_inches=0, dpi=300)
#        pylab.close('all')
#                     
#    except:
#        print('You sir, {0}, are unhelpful'.format(patient.patient_id))
#        

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