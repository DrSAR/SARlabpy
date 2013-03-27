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

    if len(pixdims) != 3:
        raise ValueError('unexpected value for VisuCoreDim')
        
    fktr = numpy.tile(pixdims,3)
    R_visupars = M_inv_RAS.reshape(9) * fktr
    print('visupars:\n{0}'.format(R_visupars.round(3)))
    print('-'*80)
    
    if not numpy.allclose(R_nii, R_visupars, atol=1e-4, rtol=1e-4):
        raise AssertionError('Rotation Matrices differ')

    # quaternion approach
#    b = nii_img.get_header()['quatern_b']
#    c = nii_img.get_header()['quatern_c']
#    d = nii_img.get_header()['quatern_d']
#    a=numpy.sqrt(1-b**2-c**2-d**2)
#    quatern_mat = nibabel.quaternions.quat2mat([a, b,c,d])
#    M_inv_RAS = LPS_2_RAS * numpy.array(M_inv.reshape(9))
#    print quatern_mat.round(5).reshape(9)
#    print M_inv_RAS.round(5)[0]
#    if numpy.allclose(quatern_mat.reshape(9),  M_inv_RAS[0],
#                      rtol=1e-3, atol=1e-3):
#        print ('close!')
#    print('nifti-quat = [{0} {1} {2} {3}]'.format(a,b,c,d))
#    print(nibabel.quaternions.mat2quat(M_inv_RAS))
    
    
    shift_nii = nii_img.get_affine()[0:3,3]
    print('affine shift     : {0}'.format(shift_nii.round(3)))
    shift_visu = v.VisuCorePosition[0,:] * [-1, 1, 1]
    print('VisuCorePosition : {0}'.format(shift_visu))
    # I think the DICOM exporter screwed this one up. Where, I wonder,
    # is (0,0,0)?    
    
    
    # all corners
    TL = numpy.matrix([0,0,0]).reshape(3,1)
    TR = numpy.matrix([v.VisuCoreSize[0],0,0]).reshape(3,1)
    BL = numpy.matrix([0,v.VisuCoreSize[1],0]).reshape(3,1)
    BR = numpy.matrix([v.VisuCoreSize[0],v.VisuCoreSize[1],0]).reshape(3,1)
#    print((numpy.matrix(R_visupars.reshape(3,3))*TL).reshape(3))
#    print((numpy.matrix(R_visupars.reshape(3,3))*TR).reshape(3))
#    print((numpy.matrix(R_visupars.reshape(3,3))*BL).reshape(3))
    print((0.5*numpy.matrix(R_visupars.reshape(3,3))*BR).reshape(3))

if __name__ == '__main__':

    dirname = os.path.join(sarpy.dataroot,'nii')
    nii_list = glob.glob(os.path.join(dirname, '*.nii.gz'))
    
    visu_dicts = {}
    Exp = sarpy.Experiment('PhantomOrientation.iY1')    
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
                