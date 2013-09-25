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


class bulk_analyzer(object):
    def __init__(self,
                 master_list_fname):

    # open masterlist and iterate over all patients and scans
        with open(os.path.join(os.path.expanduser('~/sdata'),
                               master_list_fname),'r') as master_file:
            json_str = master_file.read()
            self.master_list = json.JSONDecoder(
                            object_pairs_hook=collections.OrderedDict).decode(json_str)
    # this is where we keep a record of the results
        self.results={}

    def scan_criterion(self, scan_fname):
        return true

    def process_params(self, scan_fname):
        return {}
   
    def analysis_func(self, scan_fname):
        return 0

    def process(self):
        for pat_lbl, pat in self.master_list.iterlist():
            for scn_lbl, scn_fname in pat.iteritems():
                if self.scan_criterion(scn_fname):
                    # this is a scan we should analyze!
                    kwargs = scan_params(scn)
                    self.results[scn_lbl] = self.analysis_func(scn, **kwargs)
        return self.results

import re
import sarpy

# below are example function that achieve the minimum for a run of analysis
class dce_NR_counter(bulk_analyzer): # we inherit from the generic analyzer
                                     # object
    def scan_criterioni(self, scan_fname):
        scn = sarpy.Scan(scan_fname)
        if re.search('.*DCE',scn.acqp.protocol_name):
            return True
        else:
            return False
    def analysis_func(self, scan_fname):
        return scn.acqp.NR

res = dce_NR_counter('NecS3').process()

