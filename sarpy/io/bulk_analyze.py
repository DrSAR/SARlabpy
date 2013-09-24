# bulk analysis assuming a certain set of function.
# typical call:
# bulk_analyze(master_list,     # filename
#              scan_criterion,  # function handle, returns T for analysable scan
#              scan_params,     # the function returns a dictionary for: 
#                                    - parameters that might change from scan to scan
#                                    - parameters that are the same for all scans
#              analysis_func)   # function that takes as parameter the:
#                               #      - scan (string)
#                               #      - a dictionary of parameters

import os
import collections
import json

def bulk_analyze(master_list_fname,
                 scan_criterion,
                 scan_params,
                 analysis_func):

    # open masterlist and iterate over all patients and scans

    with open(os.path.join(os.path.expanduser('~/sdata'),master_list_fname),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(object_pairs_hook=collections.OrderedDict
                                      ).decode(json_str)    

    for pat_lbl, pat in master_list.iterlist():
        for scn_lbl, scn_fname in pat.iteritems():
            if scan_criterion(scn_fname):
                # this is a scan we should analyze!
                kwargs = scan_params(scn)
                analysis_func(scn, **kwargs)


import re
import sarpy

# below are example function that achieve the minimum for a run of analysis
def this_is_dce(fname):
    scn = sarpy.Scan(fname)
    if re.search('.*DCE',fname.acqp.protocol_name):
        return True
    else:
        return False
