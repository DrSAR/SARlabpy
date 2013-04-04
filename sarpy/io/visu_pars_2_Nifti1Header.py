# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 00:36:03 2013

@author: stefan
"""

import numpy
import nibabel

import logging
logger=logging.getLogger('sarpy.io.visu_pars_2_Nifti1Header')

def visu_pars_2_Nifti1Header(visu_pars):
    '''
    Take visu_pars header and extract all useful information for  Nifti1Header.
    This attempts to get all the geometry information out without requiring 
    any special flip etc in the data beyond the initial column-major-to-row- 
    major flip when reading in the 2dseq file.
    
    :param dict visu_pars:
        visu_pars structure as previded by readJCAMP
    
    Note: this results in  a headr that's different from the one obtained by 
    exporting BRUKER data to dicom. Not clear why that is...
    
    `Original Header information <http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/`_
    But `this site is by far the best resource for info on the Nifti header:
    <http://brainder.org/2012/09/23/the-nifti-file-format/>`

    There are 3 different methods by which continuous coordinates can be 
    attached to voxels. We are using method 2:
    METHOD 2 (used when qform_code > 0, which should be the "normal" case):
    `Ref: <http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html#ref3>`_
    ---------------------------------------------------------------------
    The (x,y,z) coordinates are given by the pixdim[] scales, a rotation
    matrix, and a shift.  This method is intended to represent
    "scanner-anatomical" coordinates, which are often embedded in the
    image header (e.g., DICOM fields (0020,0032), (0020,0037), (0028,0030),
    and (0018,0050)), and represent the nominal orientation and location of
    the data.  This method can also be used to represent "aligned"
    coordinates, which would typically result from some post-acquisition
    alignment of the volume to a standard orientation (e.g., the same
    subject on another day, or a rigid rotation to true anatomical
    orientation from the tilted position of the subject in the scanner).
    The formula for (x,y,z) in terms of header parameters and (i,j,k) is:


    [ x ]   [ R11 R12 R13 ] [        pixdim[1] * i ]   [ qoffset_x ]
    [ y ] = [ R21 R22 R23 ] [        pixdim[2]  j ] + [ qoffset_y ]
    [ z ]   [ R31 R32 R33 ] [ qfac  pixdim[3] * k ]   [ qoffset_z ]

    In methods 2 and 3, the (x,y,z) axes refer to a subject-based coordinate system,
   with
   +x = Right  +y = Anterior  +z = Superior.
   This is a right-handed coordinate system.  However, the exact direction
   these axes point with respect to the subject depends on qform_code

    Examples:
        >>> import tempfile
        >>> scn = Scan(os.path.join('readfidTest.ix1','3'))
        >>> fname = os.path.join(tempfile.gettempdir(),'readfid.nii')
        >>> print('writing tempfile %s' % fname) #doctest:+ELLIPSIS
        writing tempfile ...readfid.nii
        >>> scn.pdata[0].write2nii(fname)

    '''
    header = nibabel.nifti1.Nifti1Header()
    # setting units: it's mm = 2, s = 8, ms = 16
    header.set_xyzt_units(xyz=2,t=16) 
    # Still trying to figure out what the easiest way to do this is.
    header.set_dim_info(freq=0, phase=1, slice=2) 
    
    # Potentially non-existent settings, hence we wrap in try blocks
    try:
        # check whether all slopes are the same 
        # (i.e. set of list has length 1) -> exception
        #TODO: Deal with case when slope/offset are not identical throughout            
        if len(set(visu_pars.VisuCoreDataSlope)) > 1:
            raise ValueError("Don't know how to deal with VisuCoreDataSlope"+
                            "that vary from frame to frame")
        if len(set(visu_pars.VisuCoreDataOffs)) > 1:
            raise ValueError("Don't know how to deal with VisuCoreDataOffs"+
                            "that vary from frame to frame")
                
        slope = visu_pars.VisuCoreDataSlope[0] 
        inter = visu_pars.VisuCoreDataOffs[0]
    except AttributeError:
        logger.warn('Could not set Data slope, or Intercept\n'+
                'assuming identity.')
        slope = 1
        inter = 0

    header.set_slope_inter(slope = slope, inter = inter)

    # Let's figure out the rotation matrix. You ready? Here we go ...
    try:
        rot_mat = visu_pars.VisuCoreOrientation.flatten()
        for i in xrange(9):
            if len(set(rot_mat[i::9])) > 1:
                raise ValueError("Different slicepacks with different " +
                            "orientations!")
                            
        M = numpy.matrix(visu_pars.VisuCoreOrientation[0]).reshape(3,3)
        M_inv = M.I
        # success of the following line also hinges on the ritually
        # correct sacrifice of a chicken over the keyboard (axially).
        LPS_2_RAS = numpy.array([-1,1,-1, -1,1,-1, 1,-1,1])
        M_inv_RAS = LPS_2_RAS * numpy.array(M_inv.reshape(9))

        pixdims = numpy.array(visu_pars.VisuCoreExtent).astype('float')/ \
                  numpy.array(visu_pars.VisuCoreSize)
        # swap in-plane coordinates
        pixdims[0], pixdims[1] = pixdims[1], pixdims[0]
        # for 2D we still need to figure out the 3rd dimension
        if visu_pars.VisuCoreDim == 2:
            # check distance of neigbouring Frames
            if visu_pars.VisuCoreFrameCount == 1:
                d = visu_pars.VisuCoreFrameThickness
            else:
                p1 = numpy.array(visu_pars.VisuCorePosition[0])                
                p2 = numpy.array(visu_pars.VisuCorePosition[1])
                d = numpy.sqrt(((p1 - p2)**2).sum())
            pixdims = numpy.hstack([pixdims, d])
            
        if len(pixdims) != 3:
            raise ValueError('unexpected value for VisuCoreDim')
            
        R_visupars = M_inv_RAS.reshape(9) * numpy.tile(pixdims,3)

    except AttributeError:
        logger.warn('Could not set rotation \nassuming Identity.')
        R_visupars = numpy.eye(3).reshape(9)

    # Now we have to do the rotation business correctly. So far we have a
    # rotation matrix under the assumption of the somewhat screwy sign
    # conventions that exist for LPS systems (not RAS as in nifti)

    # Are we dealing with ax, sag, cor? Find out by looking at the 
    # through-plane vector and check where it predominantly points to. 
    # Note that we have to get rid of the pixel scaling.
    through_plane_vctr = R_visupars[2::3] / pixdims[2]
    ori_num = numpy.where(abs(through_plane_vctr) == 
                          max(abs(through_plane_vctr)))[0][0]                          
    # ori = ['sag', 'cor', 'ax'][ori_num]
    # depending on orientation we have to adjust which axis the VisuCorePosition
    # refers to. In-plane coordinates are x/y for axial, y/z for sag and x/z for
    # cor.
    if ori_num == 0: 
        # sagittal
        M = numpy.matrix([visu_pars.VisuCoreOrientation[0][3:6], 
                          -visu_pars.VisuCoreOrientation[0][0:3],
                          -visu_pars.VisuCoreOrientation[0][6:9]]).reshape(3,3)
    elif ori_num == 1:
        # coronal
        M = numpy.matrix([visu_pars.VisuCoreOrientation[0][3:6], 
                          -visu_pars.VisuCoreOrientation[0][0:3],
                          visu_pars.VisuCoreOrientation[0][6:9]]).reshape(3,3)
    elif ori_num == 2:
        # axial
        M = numpy.matrix([-visu_pars.VisuCoreOrientation[0][3:6], 
                          -visu_pars.VisuCoreOrientation[0][0:3],
                          visu_pars.VisuCoreOrientation[0][6:9]]).reshape(3,3)

    R_visupars = numpy.array(M.I).reshape(9) * numpy.tile(pixdims,3)

    try:
        qoffset =  visu_pars.VisuCorePosition[0,:] * [-1, -1, 1]
    except AttributeError:
        logger.warn('Could not find positional offset\nassuming zero.')
        qoffset = [0,0,0]
        
    aff = numpy.empty((4,4))
    aff[0:3,0:3] = R_visupars.reshape(3,3)
    aff[3,:] = [0, 0, 0, 1]
    aff[0:3,3] = qoffset
    header.set_qform(aff, code='scanner')
    header.set_sform(aff, code='scanner')

    return aff, header
    
    
if __name__ == '__main__':
    print('no doctests available yet')
