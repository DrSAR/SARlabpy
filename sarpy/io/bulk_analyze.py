'''
Definition of classes for bulk analysis.

The general idea is to instantiate a class that most suits your needs
for a given type of analysis and overload the methods that need tweaking.
'''
# bulk_analyze(masterlist,     # filename
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

class BulkAnalyzer(object):
    '''Most basic analyzer class that will be used to create 
    more elaborate classes - unlikely to be instantiated directly '''
    def __init__(self,
                 masterlist_fname):
        ''' init of the bare minimum:
            - find masterlist
            - think about where to store reslults ?!'''
    # open masterlist and iterate over all patients and scans
        fname = os.path.join(os.path.expanduser('~/sdata'),
                             masterlist_fname, masterlist_fname+'_updated.json')
        try:
            master_file = open(fname,'r')
        except IOError:
            fname = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_fname, masterlist_fname+'.json')
            try:
                master_file = open(fname,'r')
            except IOError:
                print('Could not open masterlist file %s' % fname)
                raise
            else:
                print('Warning: could not find updated masterlist')
                print('but I loaded the default instead')
      
        with master_file:
            json_str = master_file.read()
            self.masterlist = json.JSONDecoder(
                            object_pairs_hook=collections.OrderedDict).decode(json_str)

    # this is where we keep a record of the results
        self.results={}

    def scan_criterion(self, pat_lbl, scn_lbl):
        ''' method to return either
        (i) a scan_fname if pat_lbl and scn_lbl exist
        (ii) None otherwise'''
        scn_fname = self.masterlist[pat_lbl][scn_lbl][0]
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

    def store_result(self, result, 
                     scn=None):
        #scn.adata[self.adata_lbl]=result
        return 0

    def process(self):
        ''' This method will get called to run the processing.
        Currently there is no parallelization involved'''
        for pat_lbl, pat in self.masterlist.iteritems():
            for scn_lbl, scn_details in pat.iteritems():
                scn_2_analyse = self.scan_criterion(pat_lbl, scn_lbl)
                if scn_2_analyse is not None:
                    # this is a scan we should analyze!
                    kwargs = self.process_params(scn_2_analyse)
                    self.store_result(self.analysis_func(scn_2_analyse,
                                                          **kwargs),
                                      scn_2_analyse)

import re
import sarpy

class AnalyzerByScanLabel(BulkAnalyzer):
    '''Analyzer class that will only analyze scans whose masterlist
    scan label matches (re.search) the scan_label given upon initialization.
    To make use of this class, one should inherit from it and should only
    require to change method analysis_func. Then, simply instantiate and
    process()'''
    def __init__(self, masterlist_fname, scan_label=None):
        super(AnalyzerByScanLabel, self).__init__(masterlist_fname)
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
class DCE_NR_counter(BulkAnalyzer):
    '''example of a class that finds all protocols with DCE in them 
    and prints the NR value from the acqp file. Simple run like so:
    >>>> DCE_NR_counter('NecS3').process()
    '''
    def scan_criterion(self, pat_lbl, scn_lbl): 
        scan_fname = self.masterlist[pat_lbl][scn_lbl][0]
        if scan_fname:
            scn = sarpy.Scan(scan_fname)
            if re.search('.*DCE',scn.acqp.ACQ_protocol_name):
                return scn
        return None

    def analysis_func(self, scn, **kwargs):
        print('NR (%s) = %s' % (scn.shortdirname,scn.acqp.NR))
        return scn.acqp.NR

#res = DCE_NR_counter('NecS3').process()

