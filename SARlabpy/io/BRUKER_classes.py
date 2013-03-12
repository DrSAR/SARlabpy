# -*- coding: utf-8 -*-
"""
Class Definitions for BRUKER data
"""
import os
import re
import glob
import BRUKERIO

import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

dataroot = os.path.expanduser('~/data')

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
                            
class JCAMP_file(object):
    '''
    Represents a JCAMP encoded parameter file.
    
    Parameters become attributes in this class
    '''
    def __init__(self, filename, lazy=True):
        if not os.path.isfile(filename):
            raise IOError('File "%s" not found' % filename)
        self.filename = filename
        self.__yet_loaded = False
        if not lazy:
            self.__forced_load()
            
    def __getattr__(self, attr):
        '''
        intercepts only undefined attributes. this will be used to 
        trigger the loading of the data
        '''
        if not self.__yet_loaded:
            self.__forced_load() 
        return object.__getattribute__(self, attr)
        
    def __forced_load(self):
        logger.info('loading %s now forced' % self.filename)        
        acqp = BRUKERIO.readJCAMP(self.filename)
        for k,v in acqp.iteritems():
            self.__dict__[k] = v
        self.__yet_loaded = True
        
class ACQP_file(JCAMP_file):
    '''
    Thin specialization of JCAMP_file. Mostly just a place holder and
    a place to store the fact that these files are called "acqp".
    '''
    def __init__(self, dirname, lazy=True):
        super(ACQP_file, self).__init__(os.path.join(dirname,'acqp'), 
                                        lazy=lazy)

class METHOD_file(JCAMP_file):
    '''
    Thin specialization of JCAMP_file. Mostly just a place holder and
    a place to store the fact that these files are called "method".
    '''
    def __init__(self, dirname, lazy=True):
        super(METHOD_file, self).__init__(os.path.join(dirname,'method'),
                                          lazy=lazy)
                                          
class dict2obj(object):
    '''
    Helper class to turn dictionary keys into class attributes.
    It's so much nicer to refer to them that way.
    '''
    def __init__(self,dictionary):
        for k,v in dictionary.iteritems():
            self.__dict__[k] = v

class PDATA_file(object):
    '''
    Initialize a processed data set that usually sits in */pdata/[1-9].

    __init__(filename) expects a directory name in pdata. It will look
    for the 2dseq file, the d3proc and th reco file.
    '''
    def __init__(self, filename, lazy=True):
        self.__yet_loaded = False
        self.filename = filename 
        if not lazy:
            self.__forced_load()
    def  __forced_load(self):
        logger.info('loading 2dseq now forced %s:' % self.filename)
        pdata = BRUKERIO.read2dseq(self.filename)
        self.__yet_loaded = True
        self.data = pdata['data']
        for k,v in pdata['header'].iteritems():
            self.__dict__[k] = dict2obj(v)

    def __getattr__(self, attr):
        '''
        intercepts only undefined attributes. this will be used to 
        trigger the loading of the data
        '''
        if not self.__yet_loaded:
            self.__forced_load() 
        return object.__getattribute__(self, attr)                                          
        
class FID_file(object):
    '''
    Initialize an fid object whih should sit in the scan root directory.
    Really, this is a decorator class.
    
    __init__() does not expect a directory name. It will look
    for the fid file only when the attribue is being accessed. (see __get__)
    The expectation is that the class implementing this decorator will 
    provide an attribute dirname to load.
    '''
    def __init__(self):
        self.__yet_loaded = False
                
    def __set__(*args): raise AttributeError('read only!')
        
    def __get__(self, instance, value):
        if not self.__yet_loaded:
            self.__forced_load(instance.dirname) 
        return self.fid
        
    # TODO: currently, header files acqp and method are loaded twice.
    #       need to remove the redundancy. Somehow provide the headers
    #       if possible.
    def __forced_load(self, filename):
        logger.info('delayed loading of "%s" now forced ...' % filename)
        fid = BRUKERIO.readfid(os.path.join(filename,'fid'))
        self.fid = fid['data']
        self.__yet_loaded = True
        



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
        >>> dir(exp3)
        ['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'acqp', 'dirname', 'fid', 'method', 'pdata']
        >>> exp3.acqp.ACQ_protocol_name
        'FLASH_bas(modified)'
        >>> exp3.acqp.ORIGIN
        'Bruker BioSpin MRI GmbH'
        >>> exp3.acqp.DATE
        'Thu Feb 21 19:01:26 2013 PST (UT-8h)  '

    '''
    fid = FID_file()
    def __init__(self, root, absolute_root=False, lazy=True):
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

        # see whether we can find a useful BRUKER file
        try:
            self.acqp = ACQP_file(self.dirname, lazy=lazy)
        except IOError:
            self.acqp = None
        try:
            self.method = METHOD_file(self.dirname, lazy=lazy)
        except IOError:
            self.method = None
        
        # let's see whether any data has been processed
        self.pdata=[]
        for f in natural_sort(glob.glob(os.path.join(self.dirname,'pdata','*'))):
            try:
                self.pdata.append(PDATA_file(f))
            except IOError:
                pass

        # for success we need at least an image file.
        # or some acqp/method files as a miminum
        if (all([not self.acqp, not self.method, not self.pdata])):
            raise IOError(
                ('Directory "{0}" did not contain any of the typical BRUKER '+
                'files').format(self.dirname))
                
        # this feels kludgy but this will cause AttributeErrors when the 
        # user tries to access objects that don't exist (rather than getting
        # a None)
        for attr in ['acqp', 'method','pdata']:
            if not self.__dict__[attr]:
                self.__delattr__(attr)
                
    def __str__(self):
        return 'Scan object: {0}'.format(self.dirname)
        
    def __repr__(self):
        return ('Scan object: Scan("{0}")\n'+ 
                '   protocol: {1}\n'+
                '   TE = {2}\n'+
                '   TR = {3}\n'+
                '   FA = {4}\n'+
                '   FoV = {5}\n'+
                '   matrix = {6}').format(self.dirname, 
                        self.acqp.ACQ_protocol_name,
                        self.acqp.ACQ_echo_time,
                        self.acqp.ACQ_repetition_time,
                        self.acqp.ACQ_flip_angle,
                        self.acqp.ACQ_fov,
                        self.acqp.ACQ_size)
                
        
class Study(object):
    '''
    A study in BRUKER parlance is a collection of scans performed on
    on subject (typically within a day, without interruption). A study
    directory is expected to contain a JCAMP file 'subject' and have one
    or more subdirectories which are scans.
    '''
    def __init__(self, root, absolute_root=False, lazy=True):
        '''
        Initialize the Study through loading metadata from the subject file.
        
        :param string filename: directory name for the BRUKER study
        :param boolean lazy: 
            do not investigate the number and validity of scans contained 
            with the study directory
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
                
        self.subject = JCAMP_file(os.path.join(self.dirname,'subject'), 
                                      lazy=lazy)
        self.__yet_loaded = False
         
    def __getattr__(self, attr):
        '''
        intercepts only undefined attributes. this will be used to 
        trigger the loading of the data
        '''
        if not self.__yet_loaded:
            self.__forced_load() 
        return object.__getattribute__(self, attr)
        
    def __forced_load(self):
        logger.info('loading %s now forced' % self.dirname)        
        self.scans = []
        eligible_dirs = natural_sort(os.listdir(self.dirname))
        for fname in eligible_dirs:
            filename = os.path.join(self.dirname, fname)
            if (os.path.isdir(filename) and
                re.match('[0-9]+', fname)):
                self.scans.append(Scan(filename, lazy=True))
        self.__yet_loaded = True
            
        
class StudyCollection(object):
    '''
    A StudyCollection can have multiple studies (in BRUKER speak) which can in 
    turn have multiple scans. It can be initialised by pointing it 
    to a scan or study. This should be a superclass to, e.g. the Patient and 
    the Experiment class.
    '''
    
    def __init__(self, lazy=True):
        '''
        Initialize without loading any data.
        
        :param boolean lazy: 
            do not look for other studies, the number and validity of 
            scans contained with the study directory, this parameter is
            handed onwards when adding studies
        '''
        self.study_instance_uids = []
        self.studies = []
        self.lazy = lazy
        
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
                logger.warn('study previously added')
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

      
class Patient(StudyCollection):
    '''
    A Patient is a special Collection of studies in that the subject_id
    has to agree from one study to the next
    '''
    def __init__(self, patient_name, lazy=True):
        super(Patient, self).__init__(lazy=lazy)
        self.patient_id = None
        
        searchdir = os.path.join(dataroot, patient_name) + '*'
        directories = natural_sort(glob.glob(searchdir))
        for dirname in directories:
            study = Study(dirname, lazy=self.lazy)
            if not self.patient_id:
                self.patient_id = study.subject.SUBJECT_id
            if self.patient_id != study.subject.SUBJECT_id:
                logger.warning(
                    'trying to load studies of different SUBJECT_id(s)')
            else:
                self.add_study(study)
        
  
class Experiment(StudyCollection):
    '''
    An Experiment is not BRUKER terminology. It is a grouping of several 
    studies (that have scans themselves) by various criteria.
    Typically, the patient name contains a root pointing to a common goal
    of the performed scans. It is a bit like a patient with the relaxation
    of the requirement of having identical SUBJECT_ids.
    '''
    def __init__(self, root=None, absolute_root=False, lazy=True):
        super(Experiment, self).__init__(root)
        if root:
            self.find_studies(root=root, absolute_root=absolute_root)

    def find_studies(self, root=None, absolute_root=False):
        if absolute_root:
            searchdir = root
        else:
            searchdir = os.path.join(dataroot, root) + '*'

        directories = natural_sort(glob.glob(searchdir))
        for dirname in directories:
            study = Study(dirname, lazy=self.lazy)
            self.add_study(study)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    import numpy
    fname = os.path.expanduser('~/data/readfidTest.ix1/3')
    x = Scan(fname)
    print('mean = {0}'.format(numpy.mean(x.fid)))
