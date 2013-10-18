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
import time
import collections
import json
import re
import numpy
import IPython.parallel
import sarpy


class BulkAnalyzer(object):
    '''Most basic analyzer class that will be used to create 
    more elaborate classes - unlikely to be instantiated directly 
    
    only scans that fit the scan_label regexp criterion, i.e. whose masterlist
    scan label matches (re.search) the scan_label given upon initialization
    will be analysed.  
    '''
    def __init__(self,
                 masterlist_fname=None,
                 adata_lbl='testing',
                 force_overwrite=False,
                 scan_label=None):
        ''' init of the bare minimum:
            - masterlist (*_updated.json is preferrred over *.json)
            - adata_lbl ... label under which to store the analysed result
            - force_overwrite ... handed to keyword force in call to store_adate
            - scan_label (default ".*"
        '''
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
        self.adata_lbl = adata_lbl
        self.force_overwrite = force_overwrite
        
        if scan_label is None:
            self.scan_label = '.*'
        else:
            self.scan_label = scan_label


    def scan_criterion(self, pat_lbl, scn_lbl):
        ''' method to return either
        (i) a scan_fname if pat_lbl and scn_lbl exist
        (ii) None otherwise'''
        scn_fname = self.masterlist[pat_lbl][scn_lbl][0]
        if scn_fname and re.search(self.scan_label, scn_lbl):
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
        scn.store_adata(key=self.adata_lbl, 
                        data=result,
                        force=self.force_overwrite)

    def process(self):
        ''' This method will get called to run the processing.
        Currently there is no parallelization involved'''

        start1 = time.time()
        for pat_lbl, pat in self.masterlist.iteritems():
            for scn_lbl, scn_details in pat.iteritems():
                scn_2_analyse = self.scan_criterion(pat_lbl, scn_lbl)
                if scn_2_analyse is not None:
                    # this is a scan we should analyze!
                    kwargs = self.process_params(scn_2_analyse)
                    result = self.analysis_func(scn_2_analyse,
                                                          **kwargs)
                    # try:
                    #     self.store_result(result, scn_2_analyse)
                    # except AttributeError as e:
                    #     print(scn_2_analyse)
                    #     print(e)                        
        end1 = time.time()
        print 'Without parallelization : {0} s'.format(end1 - start1)   


class ParallelBulkAnalyzer(BulkAnalyzer):
    def __init__(self, 
                 masterlist_fname=None,
                 adata_lbl='testing',
                 force_overwrite=False,
                 ipython_profile='sarlab'):
        super(ParallelBulkAnalyzer, self).__init__(masterlist_fname=masterlist_fname,
                                                   adata_lbl='testing',
                                                   force_overwrite=force_overwrite)
        self.clients = IPython.parallel.Client(profile=ipython_profile)
        self.dview = self.clients[:]
        print(self.clients.ids)

        self.dview.execute('import sys,os')
        
        self.dview.execute("sys.path.append(os.path.join("+
                           "os.path.expanduser('~'),'sarpy'))")        
                
        # ensuring all engines have the same version of the imports
        with self.dview.sync_imports():
            import os
            import sys
            import sarpy
#            import sarpy.fmoosvi.parallel
            import sarpy.fmoosvi.analysis
            
        self.lv = self.clients.load_balanced_view()   # this object represents the engines (workers)
        self.lv.block = True


    def parallel_analysis_func(self, **kwargs):
        ''' Takes an iterable list and returns a lambda function for parallelization '''

        print('Shame on you for not implementing this. RTFM!')
        return None


    def process(self):
        ''' This method will get called to run the processing.
        '''
        start1 = time.time()

        list_of_scans = []
        list_of_scan_names = []
        for pat_lbl, pat in self.masterlist.iteritems():
            for scn_lbl, scn_details in pat.iteritems():
                scn_2_analyse = self.scan_criterion(pat_lbl, scn_lbl)
                if scn_2_analyse is not None:
                    list_of_scans.append(sarpy.Scan(scn_2_analyse))
                    list_of_scan_names.append(scn_2_analyse)

        func = self.parallel_analysis_func()

        results = self.lv.map(func, list_of_scan_names)
        end1 = time.time()
        print 'With parallelization : {0} s'.format(end1 - start1)    
        print results

	for res, scn in zip(results, list_of_scans):
            try:
                scn.store_adata(key=self.adata_lbl, data=res,
                                force=self.force_overwrite)
            except AttributeError as e:
                print(scn)
                print(e)                        

##################################################33333
#%%
class ParallelBulkAnalyzerFactory(BulkAnalyzer):
    def __init__(self, 
                 ipython_profile='sarlab',
                 **kwargs):
        super(ParallelBulkAnalyzerFactory, self).__init__(**kwargs)
        self.clients = IPython.parallel.Client(profile=ipython_profile)
        self.dview = self.clients[:]
        print(self.clients.ids)

        self.dview.execute('import sys,os')
        
        self.dview.execute("sys.path.append(os.path.join("+
                           "os.path.expanduser('~'),'sarpy'))")        
                
        # ensuring all engines have the same version of the imports
        with self.dview.sync_imports():
            import os
            import sys
            import sarpy
#            import sarpy.fmoosvi.parallel
            import sarpy.fmoosvi.analysis
            
        self.lv = self.clients.load_balanced_view()   # this object represents the engines (workers)
        self.lv.block = True

    @classmethod
    def from_function(cls, func, *args, **kwargs):
        print func
        print args
        print kwargs
        instance = cls(*args, **kwargs) 
        instance.parallel_analysis_func = func
        return instance

    def process(self, idx=None):
        ''' This method will get called to run the processing.
        '''
        start1 = time.time()

        list_of_scans = []
        list_of_scan_names = []
        for pat_lbl, pat in self.masterlist.iteritems():
            for scn_lbl, scn_details in pat.iteritems():
                scn_2_analyse = self.scan_criterion(pat_lbl, scn_lbl)
                if scn_2_analyse is not None:
                    list_of_scans.append(sarpy.Scan(scn_2_analyse))
                    list_of_scan_names.append(scn_2_analyse)


        if idx is not None:  
            list_of_scans = list(numpy.array(list_of_scans)[[idx]])
            list_of_scan_names = list(numpy.array(list_of_scan_names)[[idx]])

        results = self.lv.map(self.parallel_analysis_func, list_of_scan_names)
        end1 = time.time()
        print 'With parallelization : {0} s'.format(end1 - start1)    
        print results
        
        for res, scn in zip(results, list_of_scans):
            try:
                scn.store_adata(key=self.adata_lbl, data=res,
                                force=self.force_overwrite)
            except AttributeError as e:
                print(scn)
                print(e)                        

# below are example function that achieve the minimum for a run of analysis
class T1map_from_LL(BulkAnalyzer):
    '''example of a class that finds all protocols with DCE in them 
    and prints the NR value from the acqp file. Simple run like so:
    >>>> DCE_NR_counter('NecS3').process()
    '''
    def scan_criterion(self, pat_lbl, scn_lbl): 
        scan_fname = self.masterlist[pat_lbl][scn_lbl][0]
        if scan_fname:
            scn = sarpy.Scan(scan_fname)
            if re.search('.*LL',scn.acqp.ACQ_protocol_name):
                return scan_fname
        return None

    def analysis_func(self, scn, **kwargs):
        ''' Placeholder method to perform analysis on a scan '''

        print('analysing (%s): %s' % (scn.shortdirname,
                                       scn.acqp.ACQ_protocol_name))
        
        T1map, fit_dict = sarpy.fmoosvi.analysis.h_fit_T1_LL_FAind(scn.shortdirname)
        return T1map, fit_dict


class T1map_from_LLP(ParallelBulkAnalyzer):
    '''a class that finds all protocols with DCE in them 
    and prints the NR value from the acqp file. Simple run like so:
    >>>> DCE_NR_counter('NecS3').process()
    '''
    def scan_criterion(self, pat_lbl, scn_lbl): 
        scan_fname = self.masterlist[pat_lbl][scn_lbl][0]
        if scan_fname:
            scn = sarpy.Scan(scan_fname)
            if re.search('.*LL',scn.acqp.ACQ_protocol_name):
                return scan_fname
        return None

    def parallel_analysis_func(self):
        ''' Takes an iterable list and returns a lambda function for parallelization '''

        func = IPython.parallel.interactive(lambda sname:
                    sarpy.fmoosvi.analysis.h_fit_T1_LL_FAind(sname))
        return func






















