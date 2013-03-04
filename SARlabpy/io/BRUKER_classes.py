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
        
class BRUKER_fid(object):
    '''
    Initialize an fid object whih should sit in the scan root directory.
    Really, this is a decorator class.
    
    __init__(filename) expects a directory name. It will look
    for the fid file, and the acqp and method parameter files.
    '''
    def __init__(self):
        self.__yet_loaded = False
            
#    def __getattr__(self, attr):
#        '''
#        intercepts only undefined attributes. this will be used to 
#        trigger the loading of the fid data
#        '''
#        if not self._yet_loaded:
#            self._forced_load() 
#        return object.__getattribute__(self, attr)
    
    def __set__(*args): raise AttributeError
        
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
        if not lazy:
            self.__forced_load()
    def  __forced_load(self):
        logger.info('loading 2dseq now forced %s:' % self.filename)
        pdata = BRUKERIO.read2dseq(self.filename)
        self.__yet_loaded = True
        self.__dict__['data'] = pdata['data']
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

class Scan(object):
    '''
    Object to represent a BRUKER scan

    The __init__() method expects the filename to be a directory. It will
    try around a bit (filenames inside the directory) before throwing an 
    exception. It will also give up if it can't find any of the
    acqp, method, fid or one 2dseq file.
    
    Simple Example:
        
        >>> import os
        >>> fname = os.path.expanduser('~/data/readfidTest.ix1/3')
        >>> exp3 = Scan(fname)
        >>> dir(exp3)
        ['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'acqp', 'dirname', 'fid', 'method', 'pdata']
        >>> exp3.acqp.ACQ_protocol_name
        'FLASH_bas(modified)'
        >>> exp3.acqp.ORIGIN
        'Bruker BioSpin MRI GmbH'
        >>> exp3.acqp.DATE
        'Thu Feb 21 19:01:26 2013 PST (UT-8h)  '

    '''
    fid = BRUKER_fid()
    def __init__(self, filename, lazy=True):
        '''
        Is the filename a direcotry and can we at least find
        one of acqp, fid, or 2dseq
        
        :param string filename: directory name for the BRUKER data set
        '''
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
        for f in glob.glob(os.path.join(self.dirname,'pdata','*')):
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
        
#TODO: Class for a Study
class Study(object):
    def __init__(self, filename, lazy=True):
        '''
        A study in BRUKER parlance is a collection of scans performed on
        on subject (typically within a day, without interruption). A study
        directory is expected to contain a JCAMP file 'subject' and have one
        or more subdirectories which are scans/
        
        :param string filename: directory name for the BRUKER study
        :param boolean lazy: 
            do not investigate the number and validity of scans contained 
            with the study directory
        '''
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
        for fname in os.listdir(self.dirname):
            filename = os.path.join(self.dirname, fname)
            if (os.path.isdir(filename) and
                re.match('[0-9]+', fname)):
                self.scans.append(Scan(filename, lazy=True))
        self.__yet_loaded = True
            
        
#TODO: Class for a Patient
class Patient(object):
    def __init__(self):
        raise NotImplementedError


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    import numpy
    fname = os.path.expanduser('~/data/readfidTest.ix1/3')
    x = Scan(fname)
    print('mean = {0}'.format(numpy.mean(x.fid)))
    print(x.kasimir)
