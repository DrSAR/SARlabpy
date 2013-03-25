# -*- coding: utf-8 -*-

# def write_nii


"""
Original Header information:
http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/

This site is by far the best resource for info on the Nifti header:
http://brainder.org/2012/09/23/the-nifti-file-format/


     dim[0] = number of dimensions;
              - if dim[0] is outside range 1..7, then the header information
                needs to be byte swapped appropriately
              - ANALYZE supports dim[0] up to 7, but NIFTI-1 reserves
                dimensions 1,2,3 for space (x,y,z), 4 for time (t), and
                5,6,7 for anything else needed.
                
                
Useful fields [not in the correct order!]:

-Type-	-Name-		-Size-	-Description-
------------------------------------------------------------------------
int	sizeof_hdr	4B		Size of the header. Must be 348 (bytes).
char	dim_info	1B		Encoding directions (phase, frequency, slice) -x
short	dim[8]		16B		Data array dimensions, entries: num,x,y,z,t
short	slice_start	2B		First slice index
float	pixdim[8]	32B		Grid spacings, units per dimension 
float	scl_slope	4B		Data scaling, slope -x
float 	scl_offset	4B		Data scaling, offset -x
char	xyzt_units	1B		Units of pixdim[1..4] -x

--- deals with orientation ---
short	qform_code	2B		Use of the quaternion fields -x
float 	quatern_b	4B		Quaternion params
float 	quatern_c	4B		Quaternion params
float 	quatern_d	4B		Quaternion params
float 	qoffset_x	4B		Quaternion params
float 	qoffset_y	4B		Quaternion params
float 	qoffset_z	4B		Quaternion params

Orientation information

The most visible improvement of the nifti format over the previous analyze format is the ability to unambiguously store information orientation. The file standard assumes that the voxel coordinates refer to the center of each voxel, rather than at any of its corners. The world coordinate system is assumed to be RAS: +x is Right, +y is Anterior and +z is Superior, which is precisely different than the coordinate system used in analyze, which is LAS.

   There are 3 different methods by which continuous coordinates can
   attached to voxels.  The discussion below emphasizes 3D volumes, and
   the continuous coordinates are referred to as (x,y,z).  The voxel
   index coordinates (i.e., the array indexes) are referred to as (i,j,k),
   with valid ranges:
     i = 0 .. dim[1]-1
     j = 0 .. dim[2]-1  (if dim[0] >= 2)
     k = 0 .. dim[3]-1  (if dim[0] >= 3)
   The (x,y,z) coordinates refer to the CENTER of a voxel.  In methods
   2 and 3, the (x,y,z) axes refer to a subject-based coordinate system,
   with
     +x = Right  +y = Anterior  +z = Superior.
   This is a right-handed coordinate system.  However, the exact direction
   these axes point with respect to the subject depends on qform_code

WHY 3 METHODS?
   --------------
   Method 1 is provided only for backwards compatibility. The intention
   is that Method 2 (qform_code > 0) represents the nominal voxel locations
   as reported by the scanner, or as rotated to some fiducial orientation and
   location.  Method 3, if present (sform_code > 0), is to be used to give
   the location of the voxels in some standard space.  The sform_code
   indicates which standard space is present.  Both methods 2 and 3 can be
   present, and be useful in different contexts (method 2 for displaying the
   data on its original grid; method 3 for displaying it on a standard grid).

<---Note from Firas: seems we need to use Method 2--->

   In Method 2, the origin of coordinates would generally be whatever
   the scanner origin is; for example, in MRI, (0,0,0) is the center
   of the gradient coil.
   
   METHOD 2 (used when qform_code > 0, which should be the "normal" case):
   Ref: http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html#ref3
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


   The qoffset_* shifts are in the NIFTI-1 header.  Note that the center
   of the (i,j,k)=(0,0,0) voxel (first value in the dataset array) is
   just (x,y,z)=(qoffset_x,qoffset_y,qoffset_z).


   The rotation matrix R is calculated from the quatern_* parameters.
   This calculation is described below.


   The scaling factor qfac is either 1 or -1.  The rotation matrix R
   defined by the quaternion parameters is "proper" (has determinant 1).
   This may not fit the needs of the data; for example, if the image
   grid is
     i increases from Left-to-Right
     j increases from Anterior-to-Posterior
     k increases from Inferior-to-Superior
   Then (i,j,k) is a left-handed triple.  In this example, if qfac=1,
   the R matrix would have to be

     [  1   0   0 ]
     [  0  -1   0 ]  which is "improper" (determinant = -1).
     [  0   0   1 ]


   If we set qfac=-1, then the R matrix would be

     [  1   0   0 ]
     [  0  -1   0 ]  which is proper.
     [  0   0  -1 ]


   This R matrix is represented by quaternion [a,b,c,d] = [0,1,0,0]
   (which encodes a 180 degree rotation about the x-axis).


"""
import sarpy
import nibabel


def write_nii(scan_object, filename = 'test.nii'):
    
    header = nibabel.nifti1.Nifti1Header()
    
    header.set_data_shape([128, 64, 16, 1])
    header.set_dim_info(freq=0, phase=1, slice=2) # TODO: Hmm, why does nifti want to know freq/phase/slice
    header.set_xyzt_units(xyz=2,t=16) #TODO: find the units of these number. #mm = 2, s = 8, ms = 16
    header.set_slope_inter(slope = 5, inter = 0)
    
     header.set_qform(affine, code=None, strip_shears=True) #TODO: Okay this is gonna be the hardest

    
    
    
    
    
    
    
    
    
    
    
    