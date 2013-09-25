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

# This is the most basic analyzer class that will be used to create
# more elaborate classes. 
class BulkAnalyzer(object):
    def __init__(self,
                 master_list_fname):

    # open masterlist and iterate over all patients and scans
        with open(os.path.join(os.path.expanduser('~/sdata'),
                               master_list_fname, master_list_fname+'.json'),'r') as master_file:
            json_str = master_file.read()
            self.master_list = json.JSONDecoder(
                            object_pairs_hook=collections.OrderedDict).decode(json_str)
    # this is where we keep a record of the results
        self.results={}

    def scan_criterion(self, pat_lbl, scn_lbl):
        scn_fname = self.master_list[pat_lbl][scn_lbl][0]
        if scn_fname:
            return scn_fname
        else:
            return None

    def process_params(self, scn):
        ''' Placeholder method to calculate all the parameters required
        for processing'''
        return {}
   
    def analysis_func(self, scn, **kwargs):
        ''' Placeholder method to perform analysis on a scan '''
        return 0

    def process(self):
        for pat_lbl, pat in self.master_list.iteritems():
            for scn_lbl, scn_details in pat.iteritems():
                scn_2_analyse = self.scan_criterion(pat_lbl, scn_lbl)
                if scn_2_analyse is not None:
                    # this is a scan we should analyze!
                    kwargs = self.process_params(scn_2_analyse)
                    self.results[scn_lbl] = self.analysis_func(scn_2_analyse, **kwargs)
        return self.results

import re
import sarpy

class AnalyzerByScanLabel(BulkAnalyzer):
    '''Analyzer class that will only analyze scans whose masterlist
    scan label matches (re.search) the scan_label given upon initialization.
    To make use of this class, one should inherit from it and should only
    require to change method analysis_func. Then, simply instantiate and
    process()'''
    def __init__(self, master_list_fname, scan_label=None):
        super(AnalyzerByScanLabel, self).__init__(master_list_fname)
        if scan_label is None:
            raise ValueError('scan_label requires a string value')
        self.scan_label = scan_label

    def scan_criterion(self, pat_lbl, scn_lbl):
        print('testing: %s (%s) ' % (pat_lbl, scn_lbl))
        scn_fname = super(AnalyzerByScanLabel, self).scan_criterion(pat_lbl, scn_lbl)
        if re.search(self.scan_label, scn_lbl):
            return scn_fname
        return None

# below are example function that achieve the minimum for a run of analysis
class DCE_NR_counter(BulkAnalyzer): # we inherit from the generic analyzer
                                     # object
    def scan_criterion(self, pat_lbl, scn_lbl): 
        print('testing: %s (%s) ' % (pat_lbl, scn_lbl))
        scan_fname = self.master_list[pat_lbl][scn_lbl][0]
        if scan_fname:
            scn = sarpy.Scan(scan_fname)
            if re.search('.*DCE',scn.acqp.ACQ_protocol_name):
                return scn
        return None

    def analysis_func(self, scn, **kwargs):
        print('analyzing: %s' % scn)
        return scn.acqp.NR

res = DCE_NR_counter('NecS3').process()

