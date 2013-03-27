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

def nii_geom(filename=os.path.expanduser('~'),
                 visu_dict=None):
    v = visu_dict
    nii_img = nibabel.Nifti1Image.from_filename(filename)
    R_nii = nii_img.get_affine()[0:3, 0:3].reshape(9)
    pixdims = nii_img.get_header()['pixdim']
    print('affine:\n{0}'.format(R_nii.round(3)))
    
    M = numpy.matrix(v.VisuCoreOrientation[0]).reshape(3,3)
    M_inv = M.I
    LPS_2_RAS = numpy.array([-1,1,-1, -1,1,-1, 1,-1,1])
    M_inv_RAS = LPS_2_RAS * numpy.array(M_inv.reshape(9))

    pixdims = numpy.array(v.VisuCoreExtent).astype('float')/ \
              numpy.array(v.VisuCoreSize)
    # for 2D we still need to figure out the 3rd dimension
    if v.VisuCoreDim == 2:
        # check distance of neigbouring Frames
        if v.VisuCoreFrameCount == 1:
            d = v.VisuCoreFrameThickness
        else:
            p1 = numpy.array(v.VisuCorePosition[0])                
            p2 = numpy.array(v.VisuCorePosition[1])
            d = numpy.sqrt(((p1 - p2)**2).sum())
        pixdims = numpy.hstack([pixdims, d])
        k_size = v.VisuCoreFrameCount
    else:
        k_size = v.VisuCoreSize[2] 
        
    if len(pixdims) != 3:
        raise ValueError('unexpected value for VisuCoreDim')
        
    fktr = numpy.tile(pixdims,3)
    R_visupars = M_inv_RAS.reshape(9) * fktr
    print('visupars:\n{0}'.format(R_visupars.round(3)))
    print('-'*80)
    
    if not numpy.allclose(R_nii, R_visupars, atol=1e-4, rtol=1e-4):
        raise AssertionError('Rotation Matrices differ')    
    
    shift_nii = nii_img.get_affine()[0:3,3]
    print('affine shift     : {0}'.format(shift_nii.round(3)))
    # I think this is screwed up. Where, I wonder, is (0,0,0)?    
    # -> bruker and nifti appears to be assuming different corners into 
    # which to put the origin
    
    # are we dealing with ax, sag, cor? find out by looking at the through-
    # plane vector and checking where it predominantly points to. Note that
    # we have to get rid of the pixel scaling.
    through_plane_vctr = R_visupars[2::3] / pixdims[2]
    # ori = ['sag', 'cor', 'ax']
    ori_num = numpy.where(abs(through_plane_vctr) == 
                          max(abs(through_plane_vctr)))[0][0]                          
    # depending on orientation we have to adjust the origin of our images.    
    if ori_num == 1:
        addtl_offset = (R_visupars[1::3]*(v.VisuCoreSize[1] - 1) + 
                        R_visupars[2::3]*(k_size - 1))
    elif (ori_num == 2) or (ori_num ==0): # axial or sagittal
        addtl_offset = R_visupars[1::3]*(v.VisuCoreSize[1] - 1)
    else:
        raise ValueError('unexpected scan orientation detected')

    shift_visu =  v.VisuCorePosition[0,:] * [-1, -1, 1] - addtl_offset
    print('VisuCorePosition : {0}'.format(shift_visu))
    if not numpy.allclose(shift_nii, shift_visu, atol=1e-4, rtol=1e-4):
        raise AssertionError('offset vectors differ')




if __name__ == '__main__':

    dirname = os.path.join(sarpy.dataroot,'nii')
    nii_list = glob.glob(os.path.join(dirname, '*.nii.gz'))
    
    visu_dicts = {}
    Exp = sarpy.Experiment('PhantomOrientation')    
    for study in Exp.studies:
        for scan in study.scans:
            for pdata in scan.pdata:
                nii_file = os.path.join(dirname,
                                          pdata.uid()+'.nii.gz')
                if os.path.isfile(nii_file):
                    print('\n'+nii_file)
                    print('{0}\ ({1})'.format(scan.acqp.ACQ_protocol_name,                                    ' '.join(scan.method.PVM_SPackArrSliceOrient)))
                    nii_geom(nii_file, pdata.visu_pars)
                                
#    for filename in nii_list:
#        uid = re.sub('.*\/','',filename)
#        uid = re.sub('.nii.gz$','',uid)
#        load_all_nii(filename, visu_dict=visu_dicts[uid])
                