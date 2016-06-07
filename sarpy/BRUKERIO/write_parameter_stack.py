import os
import numpy
import nibabel
import sarpy

def write_parameter_stack(pat_name, 
                          params,
                          masterlist,
                          slice=None,
                          rootdir='/tmp',
                          loglevel=10):
    '''For one patient, collects all images and adatas that are
    listed in paramdict. The format of paramdict uses scan labels
    and adata labels as a(n ordered) dictionary. Typically, parameter
    stacks are written to /tmp separately for all slices unless the
    keyword parameter slice is set to a slice number

    # identify adata (referred to by its label) for a given scan (referred to by its scanlabel)
    params = [('roi-0h', ['','roi']),
              ('dce-0ha_HPG', ['auc60', 'auc300', 'auc600']),
              ('dce-0ha_Gd', ['auc60']),
              ('LL-0h+2', ['deltaR1_HPG2']),
              ('LL-0h+3', ['deltaR1_Gd'])
             ]
    expname = 'HerS11'
    masterlist = json.load(open(os.path.expanduser(os.path.join(
                           '~', 'sdata', expname, expname+'.json')),'r'))
    write_parameter_stack('HerS11Bs05', params, masterlist, slice=4)

    or, to get all slices for all patients:
    for patname in masterlist:
        for sl in arange(6):
        write_parameter_stack(patname, params, masterlist, slice=sl)
    '''
    expt = sarpy.Experiment(pat_name)
    
    # construct the array to hold the parameter stack
    first_scan = masterlist[pat_name][params[0][0]][0]
    first_maplbl = params[0][1][0]
    if first_maplbl:
        stack_shape = sarpy.Scan(first_scan).adata[first_maplbl].data.shape
    else:
        stack_shape = sarpy.Scan(first_scan).pdata[0].data.shape
    nr_layers = len([i for a in params for i in a[1]])
    img_stack = numpy.zeros(stack_shape[0:-1]+(nr_layers,))
    
    # go through all requested parameters and assemble
    if slice is None:
        slice_range = numpy.arange(stack_shape[-1])
    else:
        slice_range = [slice]
    for slice in slice_range:
        idx = 0
        for scanlbl,maplbls in params:
            for maplbl in maplbls:
                try:
                    scn = sarpy.Scan(masterlist[pat_name][scanlbl][0])
                except IOError:
                    if loglevel>=10:
                        print('Scan "{0}" not found for "{1}'.format(scanlbl, pat_name))
                else:
                    if not maplbl: # string is empty
                        if loglevel>=20:
                            print('image of shape {0}'.format(scn.pdata[0].data.shape))
                        img_stack[:,:,idx] = scn.pdata[0].data[:,:,slice]
                    else:
                        try:
                            if loglevel>=20:
                                print('adata "{0}" for scan "{1}", shape {2}'.format(maplbl, scanlbl,
                                                                             scn.adata[maplbl].data.shape))
                            img_stack[:,:,idx] = scn.adata[maplbl].data[:,:,slice]
                        except KeyError:
                            if loglevel>=10:
                                print('MISSING adata map "{0}" for scan "{1}" (pt: {2})'.
                                                  format(maplbl, scanlbl, pat_name))
                idx += 1
        nii_img = nibabel.Nifti1Image(sarpy.helpers.replace_nan_with(img_stack), None)
        fname = os.path.join(rootdir, pat_name+'-params_sl{0:=02}.nii'.format(slice+1))
        if loglevel>=20:
            print('saving to '+fname)
        nii_img.to_filename(fname)
