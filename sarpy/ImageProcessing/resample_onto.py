# -*- coding: utf-8 -*-
"""
Copyright: SARlab members, UBC, Vancouver
"""
import SimpleITK as sitk

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

if __name__ == "__main__":
    import doctest
    doctest.testmod()

