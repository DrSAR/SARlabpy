# -*- coding: utf-8 -*-
"""
Copyright: SARlab members, UBC, Vancouver
"""
import numpy
import SimpleITK as sitk
from sarpy import (visu_pars_2_Nifti1Header,
                    visu_pars_2_matrix_size)

def atleast_4d(arr):
    '''
    Return at least a 4d array and fill the missing axis with 1

    :param numpy.ndarray arr:
        numpy array

    Example:
        >>> atleast_4d(numpy.arange(720).reshape(2,3,4,5,6)).shape
        (2, 3, 4, 5, 6)
        >>> atleast_4d(numpy.arange(120).reshape(2,3,4,5)).shape
        (2, 3, 4, 5)
        >>> atleast_4d(numpy.arange(120).reshape(2,3,20)).shape
        (2, 3, 20, 1)
        >>> atleast_4d(numpy.arange(120).reshape(2,60)).shape
        (2, 60, 1, 1)
        >>> atleast_4d(numpy.arange(120).reshape(120)).shape
        (120, 1, 1, 1)
    '''
    stshape = arr.shape
    while len(stshape)<4: stshape+=(1,)
    return arr.reshape(stshape)

def resample_onto(source_fname, target_fname):
    '''
    Resamples the data in source_fname onto the grid as defined in target_fname

    :param string source_fname:
        This should be pointing to a Nifti1 file. The data and header
        will be read and used to resample data on a new grid.
    :param string source_fname:
        This should also be pointing to a Nifti1 file. Only the header with
        the enclosed geometry information will be considered to determine the
        new grid to resample source_fname.data onto.
    :return:
        resampled image as numpy.array

    Example:
        >>> import os, nibabel
        >>> import sarpy
        >>> scan_a = sarpy.Scan("PhantomOrientation.iY1/7")
        >>> scan_ref = sarpy.Scan("PhantomOrientation.iY1/4")
        >>> fname = '3D.nii.gz'
        >>> fname_ref = 'coronal.nii.gz'

        #>>> scan_a.pdata[0].export2nii(fname)
        #>>> scan_ref.pdata[0].export2nii(fname_ref)
        #>>> nibabel.load(fname).shape
        #(150, 150, 75)
        #>>> nibabel.load(fname_ref).shape
        #(256, 256, 5)
        #>>> new = resample_onto(fname, fname_ref)
        #>>> new.shape
        #(256, 256, 5)
        #>>> scan_a = sarpy.Experiment('NecS3Hs04').studies[-1].scans[7]
        #>>> scan_ref = sarpy.Experiment('NecS3Hs04').studies[-1].scans[8]
        >>> fname = '3D.nii.gz'
        >>> fname_ref = 'coronal.nii.gz'
        >>> scan_a.pdata[0].export2nii(fname)
        >>> scan_ref.pdata[0].export2nii(fname_ref)

        #>>> nibabel.load(fname).shape
        #(64, 64, 64, 2)
        #>>> nibabel.load(fname_ref).shape
        #(64, 128, 6, 25)
        #>>> new = resample_onto(fname, fname_ref)
        #>>> new.shape

    SimpleITK.Image is either 2D or 3D not more...
    check out possibilty of creating SimpleITK:image orientations from visu_para

    '''
    img_input = sitk.ReadImage(source_fname)

    ref_input = sitk.ReadImage(target_fname)
    resample_filter = sitk.ResampleImageFilter()
    resample_filter.SetDefaultPixelValue(float('nan'))
    resample_filter.SetReferenceImage(ref_input)

    # return the resampled image as an SimpleITK Image object
    output=resample_filter.Execute(img_input)
    # make a numpy array from that image
    # NB: this appears to return the thing in reverse order of dimensions
    dims = range(output.GetDimension())
    dims.reverse()
    return (sitk.GetArrayFromImage(output).transpose(dims), output)

def resample_onto_pdata(source_pdata, target_pdata):
    '''
    Resamples the data 
    
    :param string source_pdata:
        This should be pointing to pdata set. The data and header (visu_pars)
        will be read and used to resample data on a new grid.
    :param string target_pdata:
        This should also be pointing to a pdata set. Only the header (visu_pars)
        with the enclosed geometry information will be considered to determine 
        the new grid to resample source_pdata.data onto. This means that the 
        dimensions beyond the first 3 do not matter. It also means that image 
        stacks with more than just one orientation will be ignored (e.g Tri-
        Pilot)G
    :return:
        resampled image as numpy.array, or maybe ITK image? or maybe nifti1?

    Example:
        >>> import sarpy
        >>> source_pdata = sarpy.Scan("readfidTest.ix1/7").pdata[0]
        >>> source_pdata.data.shape
        (105, 133, 5, 25)
        >>> target_pdata = sarpy.Scan("readfidTest.ix1/3").pdata[0]
        >>> target_pdata.data.shape
        (105, 133, 25)
        >>> resample_onto_pdata(source_pdata, target_pdata).shape
        (105, 133, 25, 25)

    SimpleITK.Image is either 2D or 3D not more...
    check out possibilty of creating SimpleITK:image orientations from visu_para

    '''
    # the affine transform information has the pixel dimensions (sizes) rolled
    # into the. We need information from the header to extricate that.
    aff_discard, header = visu_pars_2_Nifti1Header(source_pdata.visu_pars)
    aff = header.get_qform()
    assert sum(abs(aff_discard-aff).flatten()) < 1e-5, 'affine not the same as qform'
    aff_discard, header_ref = visu_pars_2_Nifti1Header(target_pdata.visu_pars)
    aff_ref = header_ref.get_qform()
    assert sum(abs(aff_discard-aff_ref).flatten()) < 1e-5, 'affine not the same as qform'
    
    matrix_size, d1, d2 = visu_pars_2_matrix_size(source_pdata.visu_pars)
    matrix_size_ref, d1, d2 = visu_pars_2_matrix_size(target_pdata.visu_pars)

    pixdim = numpy.tile(header['pixdim'][1:4], (3,1))
    pixdim_ref = numpy.tile(header_ref['pixdim'][1:4], (3,1))
    # convert the affine transformation to a SimpleITK.Image direction, offset
    # and Origin.
    sign_matrix = numpy.array([[-1,-1,-1],[-1,-1,-1],[1,1,1]])
    direction = aff[0:3,0:3] / pixdim * sign_matrix
    direction_ref = aff_ref[0:3,0:3] / pixdim_ref * sign_matrix
    
    sign_vector = numpy.array([-1,-1,1])
    ori = aff[0:3,3] * sign_vector
    ori_ref = aff_ref[0:3,3] * sign_vector
    # setup resample filter
    resample_filter = sitk.ResampleImageFilter()
    resample_filter.SetOutputDirection(direction_ref.flatten())
    resample_filter.SetOutputOrigin(ori_ref)
    resample_filter.SetSize(matrix_size_ref.tolist())
    resample_filter.SetOutputSpacing(header_ref['pixdim'][1:4].tolist())
    resample_filter.SetDefaultPixelValue(float('nan'))

    # loop over all non-3D dimensions:
    # We have to flatten all higher dimensions and then go through them.
    flat_input_img = atleast_4d(atleast_4d(source_pdata.data)[:,:,:,...])
    dims = source_pdata.data.shape
    extra_dims = int(numpy.prod(dims[3:]))
    # output dimensions depend in a weird way on input and output dimensions:
    # the first three (spatial) dimensions are taken from the reference.
    # There are as many of these 3D volumes as there were non-spatial frames in the
    # input.
    mixed_dims = list(matrix_size_ref[0:3])
    mixed_dims.reverse() # ITK reverse spatial dimensions
    mixed_dims.append(extra_dims)
    output = numpy.empty(mixed_dims)
    for i in xrange(extra_dims):
        tobetrfx = flat_input_img[:,:,:,i]
        src = sitk.GetImageFromArray(tobetrfx.transpose([2,1,0])) # this is a 3D array
        src.SetDirection(direction.flatten())
        src.SetSpacing(header['pixdim'][1:4].tolist())
        src.SetOrigin(ori)

        res_ITK = resample_filter.Execute(src)
        res = sitk.GetArrayFromImage(res_ITK)
        output[:,:,:,i] = res
        
    # Now: reverse the spatial dimensions again from the itk to our xyz order
    dims = [2,1,0,3]
    output = output.transpose(dims)
    result_dim = list(matrix_size_ref[0:3])+list(matrix_size[3:])
    return output.reshape(result_dim).squeeze()

if __name__ == "__main__":
#    import doctest
#    doctest.testmod()
    import sarpy
    necs3 = sarpy.Experiment('NecS3').studies[7]
    dce = necs3.find_scan_by_protocol('06')[0].pdata[0]
    rare = necs3.find_scan_by_protocol('05')[0].pdata[0]
    rare.export2nii('/tmp/rare.nii')
    dce.export2nii('/tmp/dce.nii')
    # the following line only works for 3D input data (use resample_onto_pdata)
#    res_old = resample_onto('/tmp/rare.nii','/tmp/dce.nii')
    res_self = resample_onto_pdata(rare, dce)
    