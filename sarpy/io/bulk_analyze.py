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
                 adata_save_label='testing',
                 force_overwrite=False,
                 scan_label=None,
                 **kwargs):
        ''' init of the bare minimum:
            - masterlist (*_updated.json is preferrred over *.json)
            - adata_save_label ... label under which to store the analysed result
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
        self.adata_save_label = adata_save_label
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

    def process_params(self, scn_name):
        ''' Placeholder method to calculate all the parameters required
        for processing'''
        for pat,scans in self.masterlist.iteritems():
            for lbl,dbl_list in scans.iteritems():
                if dbl_list[0] == scn_name:
                    bbox = numpy.array(dbl_list[1])
        return {'bbox':bbox}
   
    def analysis_func(self, scn, **kwargs):
        ''' Placeholder method to perform analysis on a scan '''
        return 0

    def store_result(self, result, 
                     scn=None):
        scn.store_adata(key=self.adata_save_label, 
                        data=result,
                        force=self.force_overwrite)

    def process(self):
        ''' This method will get called to run the processing.
        Currently there is no parallelization involved'''

        start1 = time.time()
        for pat_lbl, pat in self.masterlist.iteritems():
            for scn_lbl, scn_details in pat.iteritems():
                scn_to_analyse = self.scan_criterion(pat_lbl, scn_lbl)
                if scn_to_analyse is not None:
                    # this is a scan we should analyze!
                    kwargs = self.process_params(scn_to_analyse)
                    result = self.analysis_func(scn_to_analyse,
                                                          **kwargs)
                    # try:
                    #     self.store_result(result, scn_to_analyse)
                    # except AttributeError as e:
                    #     print(scn_to_analyse)
                    #     print(e)                        
        end1 = time.time()
        print 'Without parallelization : {0} s'.format(end1 - start1)   


class ParallelBulkAnalyzer(BulkAnalyzer):
    def __init__(self, 
                 masterlist_fname=None,
                 adata_save_label='testing',
                 force_overwrite=False,
                 ipython_profile='sarlab'):
        super(ParallelBulkAnalyzer, self).__init__(masterlist_fname=masterlist_fname,
                                                   adata_save_label='testing',
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
            
        self.balanced = self.clients.load_balanced_view()   # this object represents the engines (workers)
        self.balanced.block = False


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
                scn_to_analyse = self.scan_criterion(pat_lbl, scn_lbl)
                if scn_to_analyse is not None:
                    list_of_scans.append(sarpy.Scan(scn_to_analyse))
                    list_of_scan_names.append(scn_to_analyse)

        func = self.parallel_analysis_func()

        results = self.balanced.map_async(func, list_of_scan_names, ordered=False)
        end1 = time.time()
        print 'With parallelization : {0} s'.format(end1 - start1)    

        for res, scn in zip(results, list_of_scans):
            try:
                scn.store_adata(key=self.adata_save_label, data=res,
                                force=self.force_overwrite)
            except AttributeError as e:
                print(scn)
                print(e)                        

##################################################33333
#%%
class ParallelBulkAnalyzerFactory(BulkAnalyzer):
    def __init__(self, 
                 ipython_profile='sarlab',
                 other_scan_label=None,
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
            
        self.balanced = self.clients.load_balanced_view()   # this object represents the engines (workers)
        self.balanced.block = False
        
        self.scan_independent_pparams = kwargs

        # if this is empty process_params will know what to do.
        self.other_scan_label = other_scan_label


    @classmethod
    def from_function(cls, func, *args, **kwargs):
        print func
        print args
        print kwargs
        instance = cls(*args, **kwargs) 
        instance.parallel_analysis_func = func
        return instance

    def process_params(self, scn_name):
        other_params = super(ParallelBulkAnalyzerFactory, self).process_params(scn_name)
        #now calculate add'l params
        if self.other_scan_label is None:
             return other_params   
        else:
            for pat,scans in self.masterlist.iteritems():
                for lbl,dbl_list in scans.iteritems():
                    if dbl_list[0] == scn_name:
                        # bingo we have found the scn_name
                        # check whether the other_scan_label exists
                        try:
                            other_scan_name = self.masterlist[pat][self.other_scan_label][0]
                        except KeyError:
                            other_scan_name = ''
            if other_scan_name:
                return dict(other_params.items() + {'other_scan_name':other_scan_name}.items())
            else:
                return None

    def process(self):
        ''' This method will get called to run the processing.
        '''
        start1 = time.time()

        list_of_scans = []
        list_of_dict_of_params = []
        for pat_lbl, pat in self.masterlist.iteritems():
            for scn_lbl, scn_details in pat.iteritems():
                scn_to_analyse = self.scan_criterion(pat_lbl, scn_lbl)
                if scn_to_analyse is not None:
                    process_params = self.process_params(scn_to_analyse)
                    if process_params is not None:
                        list_of_scans.append(sarpy.Scan(scn_to_analyse))
                        list_of_dict_of_params.append(dict(
                            {'scn_to_analyse':scn_to_analyse}.items() + 
                             process_params.items()+
                             self.scan_independent_pparams.items()))

        # the approch for asynchronous processing taking below is discussed
        # on http://stackoverflow.com/questions/19509059/processing-results-from-asyncmap-as-they-come-in
        # got a nod of approval from Min Ragan-Kelley (of ipython fame)
        if list_of_dict_of_params:
            asyncmap = self.balanced.map_async(self.parallel_analysis_func, 
                                   list_of_dict_of_params,
                                   ordered=False)

            #create original mapping of msg_ids to parameters
            msg_ids_to_parameters = dict(zip(asyncmap.msg_ids, 
                                             list_of_dict_of_params))
            msg_ids_to_scans = dict(zip(asyncmap.msg_ids, 
                                             list_of_scans))

            
            pending = set(asyncmap.msg_ids) # all queued jobs are pending
            while pending:
                try:
                    self.clients.wait(pending, 1e-3)
                except IPython.parallel.TimeoutError:
                    pass # ignore timeouterrors,  it means at least one isn't done

                # finished is the set of msg_ids that are complete
                finished = pending.difference(self.clients.outstanding)
                # update pending to exclude those that just finished
                pending = pending.difference(finished)
                for msg_id in finished:
                    # we know these are done, so don't worry about blocking
                    ar = self.clients.get_result(msg_id)
                    try:
                        ar.get()
                    except Exception as e:
                        print('%s for %i ' % (e, 
                                              msg_ids_to_parameters[msg_id]))                                  
                    else:
                        print("job id %s finished on engine %i " % 
                                (msg_id, ar.engine_id))
                        print("and results for parameter %s :" %
                                msg_ids_to_parameters[msg_id])
                        # note that each job in a map always returns a list of 
                        # length chunksize even if chunksize == 1
                        # since we need the one-to-one map of scan to result,
                        # we don't want chunks bigger than 1
                        for resdict in ar.result:
                            for k,v in resdict.iteritems():                               
                                lbl = self.adata_save_label+k
                                msg_ids_to_scans[msg_id].store_adata(
                                    key=lbl, 
                                    data=v,
                                    force=self.force_overwrite)                            
                                              
        else:
            print('no scans fit criterion and, hence, nothing to process')
        end1 = time.time()
        print 'With parallelization : {0} s'.format(end1 - start1)    

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






















