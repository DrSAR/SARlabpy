# -*- coding: utf-8 -*-
'''
Loads a slew of nii that have been generated from BRUKER exported DICOM
(BRUKER ->(datamanager)-> DICOM ->(dcm2nii)-> nii.gz )

and compares the headers in those nii files with the visu_pars information
'''

import os
import re
import glob
import numpy
import nibabel
import sarpy

'''
A case of 3D 
##$VisuUid=( 65 ) <2.16.756.5.5.100.1384712661.16188.1363828008.19>
##$VisuCoreFrameCount=1
##$VisuCoreDim=3
##$VisuCoreSize=( 3 ) 150 150 75
##$VisuCoreDimDesc=( 3 ) spatial spatial spatial
##$VisuCoreExtent=( 3 ) 80 80 40
##$VisuCoreFrameThickness=( 1 ) 40
##$VisuCoreUnits=( 3, 65 ) <mm> <mm> <mm>
##$VisuCoreOrientation=( 1, 9 ) 1 0 0 0 1 0 0 0 1
##$VisuCorePosition=( 1, 3 ) -40 -40 -16
##$VisuCoreDataMin=( 1 ) 0
##$VisuCoreDataMax=( 1 ) 32767
##$VisuCoreDataOffs=( 1 ) 0
##$VisuCoreDataSlope=( 1 ) 664.446562232508
##$VisuCoreFrameType=MAGNITUDE_IMAGE
##$VisuCoreWordType=_16BIT_SGN_INT
##$VisuCoreByteOrder=littleEndian
##$VisuCoreDiskSliceOrder=disk_normal_slice_order

A case of multi-slice 2D
-------------------------
##$VisuUid=( 65 )<2.16.756.5.5.100.1384712661.16188.1363826979.13>
##$VisuCoreFrameCount=5
##$VisuCoreDim=2
##$VisuCoreSize=( 2 ) 256 256
##$VisuCoreDimDesc=( 2 ) spatial spatial
##$VisuCoreExtent=( 2 ) 80 80
##$VisuCoreFrameThickness=( 1 ) 2
##$VisuCoreUnits=( 2, 65 ) <mm> <mm>
##$VisuCoreOrientation=( 5, 9 ) 1 0 0 0 1 0 0 0 1 1 0 0 0 1 0 0 0 1 1 0 0 0 1 0 0 0 1 1 0 0 0 1 0 0 0 1 1 0 0  0 1 0 0 0 1
##$VisuCorePosition=( 5, 3 ) -40 -40 -17.8 -40 -40 -7.8 -40 -40 2.2 -40 -40 12.2 -40 -40 22.2
##$VisuCoreDataMin=( 5 ) 0 0 0 0 0
##$VisuCoreDataMax=( 5 ) 19844 32767 28689 18439 8917
##$VisuCoreDataOffs=( 5 ) 0 0 0 0 0
##$VisuCoreDataSlope=( 5 ) 187.593280693836 187.593280693836 187.593280693836 187.593280693836  187.593280693836
##$VisuCoreFrameType=MAGNITUDE_IMAGE
##$VisuCoreWordType=_16BIT_SGN_INT

'''

def load_all_nii(dirname=os.path.expanduser('~'),
                 visu_dicts=None):
    nii_list = glob.glob(os.path.join(dirname, '*.nii.gz'))
    
    for filename in nii_list:
        a = nibabel.Nifti1Image.from_filename(filename)
        print('\n'+filename)
        print a.get_affine()[0:3, 0:3].reshape(9).round(3)
        uid = re.sub('.*\/','',filename)
        uid = re.sub('.nii.gz$','',uid)
        v = visu_dicts[uid]
        print numpy.array(v.VisuCoreOrientation[0:9]).astype('float').round(3)
        pixdims = numpy.array(v.VisuCoreExtent).astype('float')/ \
                  numpy.array(v.VisuCoreSize)
        # for 2D we still need to figure out the 3rd dimension
        if v.VisuCoreDim == 2:
            # check distance of neigbouring Frames
            if v.VisuCoreFrameCount == 1:
                d = v.VisuCoreFrameThickness
            else:
                p1 = numpy.array(v.VisuCorePosition[0:3])                
                p2 = numpy.array(v.VisuCorePosition[3:6])
                d = numpy.sqrt(((p1 - p2)**2).sum())
            pixdims = numpy.hstack([pixdims, d])

        if len(pixdims) != 3:
            raise ValueError('unexpected value for VisuCoreDim')
        
        print pixdims
        

if __name__ == '__main__':

    visu_dicts = {}
    Exp = sarpy.Experiment('PhantomOrientation')    
    for study in Exp.studies:
        for scan in study.scans:
            for pdata in scan.pdata:
                visu_dicts[pdata.visu_pars.VisuUid]=pdata.visu_pars
    
    load_all_nii(dirname = '/tmp/nii', visu_dicts=visu_dicts)