# -*- coding: utf-8 -*-
"""
Class Definitions for BRUKER data
"""
import os
import re
import glob
import BRUKERIO
import AData_classes

import logging
logger=logging.getLogger('sarpy.io.BRUKER_classes')

dataroot = os.path.expanduser('~/data')

def natural_sort(l): 
    '''
    Sort a list by a natural sort order (number ascending) and even if
    bracketed by blocks of alpha-characters.
    
    This is based on code a blog post by `Jeff Atwood <http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html>`_.
    It is also discussed on `stackoverflow <http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort>`_.
    '''    
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def strip_all_but_classname(obj, class_str):
    '''
    Nicely format the class_str that is used in some of the __str__ calls
    '''
    m = re.search("'[^']+'", str(obj.__class__))
    if m:
        class_name = m.group(0).strip("'")
    else:
        class_name = class_str
    return class_name
                           
class lazy_property(object):
    '''
    Meant to be used for lazy evaluation of an object attribute.
    property should represent non-mutable data, as it replaces itself.
    
    An example illustrates the use. First define a class that contains
    a lazy_property. That call it and observe how the calculation is 
    only performed once.
    
    >>> class Test(object):
    ...     @lazy_property
    ...     def results(self):
    ...         calcs = 42 # do a lot of calculation here
    ...         print('phew, this took a lot of effort - glad I do this so rarely...')
    ...         return calcs
    >>> A=Test()
    >>> A.__dict__
    {}
    >>> A.results
    phew, this took a lot of effort - glad I do this so rarely...
    42
    >>> A.results  # on a second call the calculation is not needed !
    42
    >>> A.__dict__  # see how results is no a proper attribute
    {'results': 42}

    '''

    def __init__(self,fget):
        self.fget = fget
        self.func_name = fget.__name__

    def __get__(self,obj,cls):
        ''' this is called on attribute access '''
        if obj is None:
            return None
        value = self.fget(obj)
        # the following overwrites the dict in the calling object
        # thereby erasing any access to this function (which)
        # calls the expensive calculation once.
        setattr(obj,self.func_name,value)
        return value

# ===========================================================
                           
class JCAMP_file(object):
    '''
    Represents a JCAMP encoded parameter file.
    
    Parameters become attributes in this class
    '''
    def __init__(self, filename):
        if not os.path.isfile(filename):
            raise IOError('File "%s" not found' % filename)
        self.filename = filename
        acqp = BRUKERIO.readJCAMP(self.filename)
        for k,v in acqp.iteritems():
            self.__dict__[k] = v
                                                  
class PDATA_file(object):
    '''
    Initialize a processed data set that usually sits in `*/pdata/[1-9]`.

    constructor expects a directory name (filename). It will look
    for the 2dseq file, the d3proc and th reco file.
    '''
    def __init__(self, filename):
        self.filename = filename 
        
    @lazy_property
    def reco(self):
        return JCAMP_file(os.path.join(self.filename,'reco'))

    @lazy_property
    def visu_pars(self):
        return JCAMP_file(os.path.join(self.filename,'visu_pars'))
                
    @lazy_property
    def d3proc(self):
        return JCAMP_file(os.path.join(self.filename,'d3proc'))

    @lazy_property
    def data(self):
        dta = BRUKERIO.read2dseq(os.path.join(self.filename),
                                 reco=self.reco.__dict__,
                                 d3proc=self.d3proc.__dict__,
                                 visu_pars=self.visu_pars.__dict__)
        return dta['data']

    def uid(self):
        return self.visu_pars.VisuUid
    
    def store_adata(self, *args, **kwargs):
        '''
        Store some secondary data for this PData scan.

        See call signature of AData_classes.AData.fromdata
        
        Would be nice to trigger an update in the adata list of the 
        parent object (Scan). This might be hard. Sounds like Traits to me.        
        '''
        return AData_classes.AData.fromdata(self, *args, **kwargs)



class Scan(object):
    '''
    Object to represent a BRUKER scan consisting, typically, of FID_file 
    and PDATA_file(s).

    The __init__() method expects the filename to be a directory. It will
    try around a bit (filenames inside the directory) before throwing an 
    exception. It will also give up if it can't find any of the
    acqp, method, fid or one 2dseq file.
    
    Simple Example:
        
        >>> import os
        >>> exp3 = Scan('readfidTest.ix1/3')
        >>> dir(exp3)    # doctest: +NORMALIZE_WHITESPACE
        ['__class__', '__delattr__', '__dict__', '__doc__', '__format__', 
         '__getattribute__', '__hash__', '__init__', '__module__', '__new__', 
         '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', 
         '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'acqp', 
         'adata', 'dirname', 'fid', 'method', 'pdata', 'pdata_uids', 
         'shortdirname', 'store_adata']

        >>> exp3.acqp.ACQ_protocol_name
        'FLASH_bas(modified)'
        >>> exp3.acqp.ORIGIN
        'Bruker BioSpin MRI GmbH'
        >>> exp3.acqp.DATE
        'Thu Feb 21 19:01:26 2013 PST (UT-8h)'
        
    Also, Scan has a pretty way of falling over when no fid file is detected:
        
        >>> Scan('readfidTest.ix1/9') # doctest:+ELLIPSIS
        Traceback (most recent call last):
            ...
        IOError:...
    '''
    @lazy_property
    def acqp(self):
        return JCAMP_file(os.path.join(self.dirname,'acqp'))

    @lazy_property
    def method(self):
        return JCAMP_file(os.path.join(self.dirname,'method'))

    @lazy_property
    def fid(self):
        kspace = BRUKERIO.readfid(os.path.join(self.dirname,'fid'),
                                  acqp=self.acqp.__dict__,
                                  method=self.method.__dict__)['data']
        return kspace
        
    @lazy_property
    def adata(self):
        adata_dict = {}
        logger.info('delayed loading of adata (list) now forced ...')
        for pdata_uid in self.pdata_uids:
            adata_potential = os.path.join(AData_classes.adataroot, pdata_uid)
            if os.path.isdir(adata_potential):
                for adata_sets in os.listdir(adata_potential):
                    adata_candidate = AData_classes.AData.fromfile(
                            os.path.join(adata_potential, adata_sets))
                    adata_dict[adata_candidate.key]=adata_candidate
        if len(adata_dict) == 0:
            logger.info('No analysis data in directory "{0}"'.
                         format(self.dirname))
        return adata_dict
        
    @lazy_property
    def pdata(self):
        pdata_list = []
        for f in natural_sort(glob.glob(os.path.join(self.dirname,'pdata','*'))):
            try:
                pdata = PDATA_file(f)
                pdata_list.append(pdata)
            except IOError:
                pass
        return pdata_list

    @lazy_property
    def pdata_uids(self):
        return [pdata.uid() for pdata in self.pdata]
    
    def __init__(self, root, absolute_root=False):
        '''
        Is the filename a direcotry and can we at least find
        one of acqp, fid, or 2dseq
        
        :param string filename: directory name for the BRUKER data set
        '''
        if absolute_root:
            filename = root
        else:
            filename = os.path.join(dataroot, root)
        if os.path.isdir(filename):
            # good, that could be it
            self.dirname = filename
        elif os.path.isdir(os.path.dirname(filename)):
            logger.warning(('Scan initialized with a filename (%s) not a' +
                            ' directory.') % filename)
            # did the punter point the filename to some file inside
            # the directory?
            self.dirname = os.path.dirname(filename)
        else:
            # no good
            raise IOError(
                'Filename "%s" is neither a directory nor a file inside one' %
                filename)
        self.shortdirname = re.sub(dataroot+os.path.sep, '', self.dirname)
        # see whether we can find an fid file
        # in all likelihood this means that an acqp and method file
        # is also present - this was true for ca 9000 scans we tested thi in
        if not(os.path.isfile(os.path.join(self.dirname,'fid'))):
            raise IOError(
                ('Directory "{0}" did not contain a BRUKER fid '+
                'file').format(self.dirname))            

        # if there are no 2dseq files (processed data) this is potentially
        # solvable of the data is retro-actively reconstructed. In the
        # meantime -> warning
        if len(glob.glob(os.path.join(self.dirname,'pdata','*','2dseq'))) == 0:
            logger.warning(('Directory "{0}" did not contain processed data '+
                         '(2dseq)').format(self.dirname))
               

    def __str__(self):
        '''
        Simple print representation of the Scan object:
            
            >>> a=Scan('readfidTest.ix1/5')
            >>> print(a)
            __main__.Scan("readfidTest.ix1/5")
        '''
        return '{0}("{1}")'.format(
                    strip_all_but_classname(self, 'Scan'),
                    self.shortdirname)
        
    def __repr__(self):
        '''
        More elaborate string representation of the Scan object:
            
            >>> a=Scan('readfidTest.ix1/5')
            >>> a
            Scan object: Scan("readfidTest.ix1/5")
               protocol: MSME_bas
               TE = [14]
               TR = [50]
               FA = 180
               FoV = [3, 4, 2.5]
               matrix = [266, 105, 25]
        '''
        try:
            
            try: 
                return ('Scan object: Scan("{0}")\n'+ 
                                            '   protocol: {1}\n'+
                                            '   TE = {2}\n'+
                                            '   TR = {3}\n'+
                                            '   FA = {4}\n'+
                                            '   FoV = {5}\n'+
                                            '   matrix = {6}\n'
                                            '   -----BS-----\n'+
                                            '   BSFreq = {7}\n'+
                                            '   BSPower = {8}').format(self.shortdirname, 
                                                    self.acqp.ACQ_protocol_name,
                                                    self.acqp.ACQ_echo_time,
                                                    self.acqp.ACQ_repetition_time,
                                                    self.acqp.ACQ_flip_angle,
                                                    self.acqp.ACQ_fov,
                                                    self.acqp.ACQ_size,
                                                    self.method.BSFreqOffset,
                                                    self.method.BSPulse[3])
            except:
                return ('Scan object: Scan("{0}")\n'+ 
                        '   protocol: {1}\n'+
                        '   TE = {2}\n'+
                        '   TR = {3}\n'+
                        '   FA = {4}\n'+
                        '   FoV = {5}\n'+
                        '   matrix = {6}').format(self.shortdirname, 
                                self.acqp.ACQ_protocol_name,
                                self.acqp.ACQ_echo_time,
                                self.acqp.ACQ_repetition_time,
                                self.acqp.ACQ_flip_angle,
                                self.acqp.ACQ_fov,
                                self.acqp.ACQ_size)
        except AttributeError:
            return self.__str__()

    def store_adata(self, pdata_idx=0, **kwargs):
        '''
        Store an AData set for one of the processed children of this scan.
        Typically the first one (pdata_idx=0)
        '''
        self.adata[kwargs['key']] = AData_classes.AData.fromdata(
                    self.pdata[pdata_idx], **kwargs)
                    
class Study(object):
    '''
    A study in BRUKER parlance is a collection of scans performed on
    on subject (typically within a day, without interruption). A study
    directory is expected to contain a JCAMP file 'subject' and have one
    or more subdirectories which are scans.
    '''
    def __init__(self, root, absolute_root=False):
        '''
        Initialize the Study through loading metadata from the subject file.
        
        :param string filename: directory name for the BRUKER study
        '''
        if absolute_root:
            filename = root
        else:
            filename = os.path.join(dataroot, root)

        if os.path.isdir(filename):
            # good, that could be it
            self.dirname = filename
        elif os.path.isdir(os.path.dirname(filename)):
            logger.warning(('Scan initialized with a filename (%s) not a' +
                            ' directory.') % filename)
            # did the punter point the filename to some file inside
            # the directory?
            self.dirname = os.path.dirname(filename)
        else:
            # no good
            raise IOError(
                'Filename "%s" is neither a directory nor a file inside one' %
                filename)
          
        self.shortdirname = re.sub(dataroot+os.path.sep, '', self.dirname)
        self.subject = JCAMP_file(os.path.join(self.dirname,'subject'))
         
    @lazy_property
    def scans(self):
        scans = []
        eligible_dirs = natural_sort(os.listdir(self.dirname))
        for fname in eligible_dirs:
            filename = os.path.join(self.dirname, fname)
            if (os.path.isdir(filename) and
                re.match('[0-9]+', fname)):
                try:
                    scans.append(Scan(filename))
                except IOError:
                    pass
        return scans
        
    def __str__(self):
        '''
        Simple print representation of the Study object
        
            >>> a=Study('readfidTest.ix1')
            >>> print(a)
            __main__.Study("readfidTest.ix1")
        '''
        return '{0}("{1}")'.format(
                    strip_all_but_classname(self, 'Study'),
                    self.shortdirname)
        
    def __repr__(self):
        '''
        More elaborate string representation of the Study object:

        >>> a=Study('readfidTest.ix1')
        >>> a
        Study object: Study("readfidTest.ix1")
           subject: Moosvi, readfidTest
           scans: 1-TriPilot-multi
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  MSME_bas
                  MSME_bas
                  MSME_turbo_bas
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  MSME_bas
                  EPI_bas
                  EPI_nav_bas
                  DtiStandard_bas
                  UTE2D
                  UTE3D
                  ZTE
                  FLASH_bas(modified)
        '''
        try:
            scan_list_repr = [scan.acqp.ACQ_protocol_name for 
                                scan in self.scans]
            return ('Study object: Study("{0}")\n'+ 
                    '   subject: {1}\n'+
                    '   scans: {2}').format(
                                self.shortdirname,
                                self.subject.SUBJECT_name_string,
                                '\n          '.join(scan_list_repr))
        except AttributeError:
            return self.__str__()
            
    def find_scan_by_protocol(self, protocol_name):
        found_scans = []

        for s in self.scans: 
            try:
                if re.match(protocol_name, s.acqp.ACQ_protocol_name):
                    found_scans.append(s)
            except AttributeError:
                print('Warning: Scan in dir %s has no acqp attribute' %str(s.shortdirname))
        return(found_scans)
                                
            
        
class StudyCollection(object):
    '''
    A StudyCollection can have multiple studies (in BRUKER speak) which can in 
    turn have multiple scans. It can be initialised by pointing it 
    to a scan or study. This should be a superclass too (?), e.g. the Patient and 
    the Experiment class.
    '''
    
    def __init__(self):
        '''
        Initialize without loading any data.
        '''
        self.study_instance_uids = []
        self.studies = []
    
    def __str__(self):
        '''
        Simple print representation of the Study object
        
            >>> a=StudyCollection()
            >>> print(a)
            __main__.StudyCollection()
        '''
        return '{0}()'.format(
                    strip_all_but_classname(self, 'StudyCollection'))
                    
    def __repr__(self):
        '''
        More elaborate string representation of the StudyCollection object:

            >>> a=StudyCollection()
            >>> a
            __main__.StudyCollection()
        '''
        lines_to_print = [self.__str__()]
        lines_to_print.extend([str(s) for s in self.studies])
        return  '\n'.join(lines_to_print)
        
    def add_study(self, study=None, by_name=None, by_uid=None):
        '''search through database to find anything that matches root
        
        :param Study study: 
            add the indicated study
        :param string by_name:
            search through database with the indicated string
        :param string by_uid:
            search through database by UID. This is potentially very slow
        '''

        done = False
        if study:       
            assert isinstance(study, Study), (
                'Object added to StudyCollection is not a Study')
            done = True

        if by_name:
            if done:
                raise ValueError(
                    'Object added to StudyCollection identified by more than '+
                    'one criterion')                    
            raise NotImplementedError
            done = True
        
        if by_uid:
            if done:
                raise ValueError(
                    'Object added to StudyCollection identified by more than '+
                    'one criterion')                    
            raise NotImplementedError
            done = True
            
        if done:
            if (study.subject.SUBJECT_patient_instance_uid in 
                self.study_instance_uids):
                logger.warning('study previously added')
            else:
                self.studies.append(study)
                try:
                    self.study_instance_uids.append(
                            study.subject.SUBJECT_patient_instance_uid)
                except AttributeError:
                    logger.warning(
                        'SUBJECT_patient_instance_uid not found in study: %s' %
                        study.dirname)
                try:
                    logger.info('study "%s" added to StudyCollection' % 
                                study.subject.SUBJECT_study_name)
                except AttributeError:
                    logger.warning(
                        'SUBJECT_study_name not found in study: %s' %
                        study.dirname)
        else:
            logger.warning('no study added to StudyCollection')

    def get_SUBJECT_id(self):
        return [x.subject.SUBJECT_id for x in self.studies]
        
    def find_scan_by_protocol(self, protocol_name):
        found_scans = []
        for study in self.studies:
            found_scans.extend(study.find_scan_by_protocol(protocol_name))
        return(found_scans)



      
class Patient(StudyCollection):
    '''
    A Patient is a special Collection of studies in that the subject_id
    has to agree from one study to the next
    '''
    def __init__(self, patient_name):
        super(Patient, self).__init__()
        self.patient_id = None
        
        searchdir = os.path.join(dataroot, patient_name) + '*'
        directories = natural_sort(glob.glob(searchdir))
        for dirname in directories:
            study = Study(dirname)
            if not self.patient_id:
                self.patient_id = study.subject.SUBJECT_id
            if self.patient_id != study.subject.SUBJECT_id:
                logger.warning(
                    'loading studies of different SUBJECT_id(s) not permitted')
            else:
                self.add_study(study)
                
    def __str__(self):
        '''
        Simple print representation of the Patient object
        
            >>> a=Patient('readfidTest')
            >>> print(a)
            __main__.Patient("readfidTest")
        '''        
        return '{0}("{1}")'.format(
                strip_all_but_classname(self, 'Patient'),
                self.patient_id)
                
    def __repr__(self):
        '''
        More elaborate string representation of the Patient object:

        >>> a=Patient('readfidTest')
        >>> a
        Patient object: Patient("readfidTest")
           studies: Study object: Study("readfidTest.ix1")
           subject: Moosvi, readfidTest
           scans: 1-TriPilot-multi
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  MSME_bas
                  MSME_bas
                  MSME_turbo_bas
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  MSME_bas
                  EPI_bas
                  EPI_nav_bas
                  DtiStandard_bas
                  UTE2D
                  UTE3D
                  ZTE
                  FLASH_bas(modified)
        '''
        try:
            study_list_repr = [study.__repr__() for study in self.studies]
            return ('Patient object: Patient("{0}")\n'+ 
                    '   studies: {1}').format(
                                self.patient_id,
                                '\n          '.join(study_list_repr))
        except AttributeError:
            return self.__str__()
  
class Experiment(StudyCollection):
    '''
    An Experiment is not BRUKER terminology. It is a grouping of several 
    studies (that have scans themselves) by various criteria.
    Typically, the patient name contains a root pointing to a common goal
    of the performed scans. It is a bit like a patient with the relaxation
    of the requirement of having identical SUBJECT_ids.
    '''
    def __init__(self, root=None, absolute_root=False):
        super(Experiment, self).__init__()
        if root:
            self.find_studies(root=root, absolute_root=absolute_root)
            
    def __str__(self):
        '''
        Simple print representation of Experiment object
        
        >>> d=Experiment('NecS3Hs10')
        >>> print(d)
        __main__.Experiment("NecS3Hs10")
        '''
        return '{0}("{1}")'.format(
                strip_all_but_classname(self, 'Experiment'),
                self.root)

    def __repr__(self):
        '''
        More elaborate string representation of the Experiment object:

        >>> a=Experiment('NecS3')
        >>> a
        Experiment object: Experiment("NecS3")
           studies: --Total (30)--
                    NecS3Hs01a.iJ1
                    ...
        '''
        try:
            
            study_list = [study.shortdirname for study in self.studies]
            
            return ('Experiment object: Experiment("{0}")\n'+ 
                    '   studies: --Total ({1})--\n'+
                    '            {2}').format(self.root, len(study_list),
                                '\n            '.join(study_list))
        except AttributeError:
            return self.__str__()
  


    def find_studies(self, root=None, absolute_root=False):
        if absolute_root:
            searchdir = root
        else:
            searchdir = os.path.join(dataroot, root) + '*'

        self.root = root
        directories = natural_sort(glob.glob(searchdir))
        for dirname in directories:
            study = Study(dirname)
            self.add_study(study)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
