import re
import numpy
import sarpy
def get_histogram(masterlist,
                  pat_regex='', # the empty regex fits all
                  scan_regex='',
                  adata_regex='',
                  slice_range=None,
                  roi_label=None, # tuple containing scan and adata label
                  verbose=False):
    '''Returns histogram for a parameter map inside a given ROI
    
    HerP2_dR1_HPG=get_histogram(masterlist, pat_regex='',
                                scan_regex='',
                                adata_regex='multiday_deltaR1_HPG',
                                slice_range=None,
                                roi_label=('roi-24h', 'roi'),
                                verbose=True)
    '''
    
    vals = {}
    for pat in [x for x in masterlist if re.search(pat_regex, x)]:
        if verbose: print(('looking at {0}'.format(pat)))
        for scn_lbl in [x for x in masterlist[pat] if re.search(scan_regex, x)]:
            if verbose: 
                print(('\t - scan {0}'.format(scn_lbl)))
            try:
                scn = sarpy.Scan(masterlist[pat][scn_lbl][0])
            except IOError:
                if verbose: print(('\t !! could not find data for {0}'.format(scn_lbl)))
            else:
                for adata_lbl in [x for x in list(scn.adata.keys()) if re.search(adata_regex, x)]:
                    target_data = scn.adata[adata_lbl].data
                    if verbose: print(('\t - adata {0} {1}'.format(adata_lbl,target_data.shape)))
                    if roi_label is not None:
                        # mask the target_data with a ROI mask
                        
                        # Not sure why this check wasn't here before, 
                        # some scans just don't exist in the masterlist because they weren't acquired
                        print(scn.shortdirname)

                        try: 
                            scn_roi = sarpy.Scan(masterlist[pat][roi_label[0]][0])

                            roi = scn_roi.adata[roi_label[1]].data

                            target_data = target_data * roi
                                                    
                        except:
                            print(('\t !! could not find data for {0} and {1}'.format(pat,scn_lbl)))
                        
                    if slice_range is None:
                        slice_range = numpy.arange(scn.adata[adata_lbl].data.shape[-1])
                    xx = target_data[:,:,slice_range]
                    vals[pat]=xx[numpy.where(numpy.isfinite(xx))]
    return vals

