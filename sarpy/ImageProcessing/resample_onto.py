# -*- coding: utf-8 -*-
"""
Copyright: SARlab members, UBC, Vancouver
"""
import numpy
import SimpleITK as sitk
from visu_pars_2_Nifti1Header import visu_pars_2_Nifti1Header


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
        >>> scan_a.pdata[0].export2nii(fname)
        >>> scan_ref.pdata[0].export2nii(fname_ref)
        >>> nibabel.load(fname).shape
        (150, 150, 75)
        >>> nibabel.load(fname_ref).shape
        (256, 256, 5)
        >>> new = resample_onto(fname, fname_ref)
        >>> new.shape
        (256, 256, 5)
        >>> scan_a = sarpy.Experiment('NecS3Hs04').studies[-1].scans[7]
        >>> scan_ref = sarpy.Experiment('NecS3Hs04').studies[-1].scans[8]
        >>> fname = '3D.nii.gz'
        >>> fname_ref = 'coronal.nii.gz'
        >>> scan_a.pdata[0].export2nii(fname)
        >>> scan_ref.pdata[0].export2nii(fname_ref)
        >>> nibabel.load(fname).shape
        (64, 64, 64, 2)
        >>> nibabel.load(fname_ref).shape
        (64, 128, 6, 25)
        >>> new = resample_onto(fname, fname_ref)
        >>> new.shape

    SimpleITK.Image is either 2D or 3D not more...
    check out possibilty of creating SimpleITK:image orientations from visu_para

    '''
    img_input = sitk.ReadImage(source_fname)
    ref_input = sitk.ReadImage(target_fname)
    resample_filter = sitk.ResampleImageFilter()
    resample_filter.SetReferenceImage(ref_input)
    # return the resampled image as an SimpleITK Image object
    output=resample_filter.Execute(img_input)
    # make a numpy array from that image
    # NB: this appears to return the thing in reverse order of dimensions
    dims = range(output.GetDimension())
    dims = dims.reverse()
    return sitk.GetArrayFromImage(output).transpose(dims)

def resample_onto2(source_pdata, target_pdata):
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
        Pilot)
    :return:
        resampled image as numpy.array, or maybe ITK image? or maybe nifti1?

    Example:
        >>> import sarpy
        >>> scan_a = sarpy.Scan("PhantomOrientation.iY1/7")
        >>> scan_ref = sarpy.Scan("PhantomOrientation.iY1/4")
        >>> source_pdata = scan_a.pdata[0]
        >>> target_pdata = scan_ref.pdata[0]
        >>> scan_a.pdata[0].data.shape
        (150, 150, 75)
        >>> scan_ref.pdata[0].data.shape
        (256, 256, 5)
        >>> new = resample_onto2(scan_a, scan_ref)
        >>> new.shape
        (256, 256, 5)

    SimpleITK.Image is either 2D or 3D not more...
    check out possibilty of creating SimpleITK:image orientations from visu_para

    '''
    # the affine transform information has the pixel dimensions (sizes) rolled
    # into the. We need information from the header to extricate that.
    aff, header = visu_pars_2_Nifti1Header(source_pdata.visu_pars)
    aff_ref, header_ref = visu_pars_2_Nifti1Header(source_pdata.visu_pars)

    pixdim = numpy.tile(header['pixdim'][1:4], (3,1))
    pixdim_ref = numpy.tile(header_ref['pixdim'][1:4], (3,1))
    # convert the affine transformation to a SimpleITK.Image direction, offset
    # and Origin.
    direction = aff[0:3,0:3] / pixdim
    direction_ref = aff_ref[0:3,0:3] / pixdim_ref
    # setup resample filter
    resample_filter = sitk.ResampleImageFilter()
    resample_filter = resample_filter.SetOutputDirection(direction_ref.flatten())
    resample_filter = resample_filter.SetOutputOrigin(aff_ref[0:3,3])
    resample_filter = resample_filter.SetSize()

    # return the resampled image as an SimpleITK Image object
    output=resample_filter.Execute(img_input)
    # make a numpy array from that image
    # loop over all non-3D dimensions:
    # We have to flatten all higher dimensions and then go through them.
    flat_img = source_pdata.data[:,:,:,-1]

    dims = source_pdata.data.shape
    extra_dims = numpy.prod(dims[3:])
    for i in xrange(extra_dims):
        src = SimpleITK.Image(flat_img[:,:,:,i]) # this is a 3D array
        src.set_direction(direction)
        src.set_size(pixdim)
        src.set_offset(aff[3,0:3].flatten())

        output = 


    img_input = sitk.ReadImage(source_fname)
    ref_input = sitk.ReadImage(target_fname)
    resample_filter = sitk.ResampleImageFilter()
    resample_filter.SetReferenceImage(ref_input)
    # return the resampled image as an SimpleITK Image object
    output=resample_filter.Execute(img_input)
    # make a numpy array from that image
    # NB: this appears to return the thing in reverse order of dimensions
    dims = range(output.GetDimension())
    dims = dims.reverse()
    return sitk.GetArrayFromImage(output).transpose(dims)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

