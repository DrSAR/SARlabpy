# -*- coding: utf-8 -*-
"""
Class Definitions for BRUKER data
"""
import os
import glob
import BRUKERIO

import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#TODO: Class for a Patient
class Patient(object):
    def __init__(self):
        raise NotImplementedError
        
#TODO: Class for a Study
class Study(object):
    def __init__(self):
        raise NotImplementedError

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
        ['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'acqp', 'fid', 'method', 'pdata']
        >>> exp3.acqp.ACQ_protocol_name
        'FLASH_bas(modified)'
        >>> exp3.acqp.ORIGIN
        'Bruker BioSpin MRI GmbH'
        >>> exp3.acqp.DATE
        'Thu Feb 21 19:01:26 2013 PST (UT-8h)  '

    '''
    def __init__(self, filename):
        '''
        Is the filename a direcotry and can we at least find
        one of acqp, fid, or 2dseq
        
        :param string filename: directory name for the BRUKER data set
        '''
        if os.path.isdir(filename):
            # good, that could be it
            dirname = filename
        elif os.path.isdir(os.path.dirname(filename)):
            logger.warning(('Scan initialized with a filename (%s) not a' +
                            ' directory.') % filename)
            # did the punter point the filename to some file inside
            # the directory?
            dirname = os.path.dirname(filename)
        else:
            # no good
            raise IOError(
                'Filename "%s" is neither a directory nor a file inside one' %
                filename)

        # see whether we can find a useful BRUKER file
        try:
            self.acqp = ACQP_file(dirname)
        except IOError:
            self.acqp = None
        try:
            self.fid = BRUKER_fid(dirname)
        except IOError:
            self.fid = None
        try:
            self.method = METHOD_file(dirname)
        except IOError:
            self.method = None
        
        # let's see whether any data has been processed
        self.pdata=[]
        for f in glob.glob(os.path.join(dirname,'pdata','*')):
            try:
                self.pdata.append(PDATA_file(f))
            except IOError:
                pass

        if not any([self.acqp, self.fid, self.method]):
            #we have not found anything
            raise IOError(
                ('Directory "{0}" did not contain any of the typical BRUKER '+
                'files').format(dirname))
                            
class JCAMP_file(object):
    '''
    Represents a JCAMP encoded parameter file.
    
    Parameters become attributes in this class
    '''
    def __init__(self, filename, lazy=True):
        if not os.path.isfile(filename):
            raise IOError('File "%s" not found' % filename)
        self.filename = filename
        self._yet_loaded = False
        if not lazy:
            self._forced_load()
            
    def __getattr__(self, attr):
        '''
        intercepts only undefined attributes. this will be used to 
        trigger the loading of the data
        '''
        if not self._yet_loaded:
            self._forced_load() 
        return object.__getattribute__(self, attr)
        
    def _forced_load(self):
        logger.info('loading %s' % self.filename)
        acqp = BRUKERIO.readJCAMP(self.filename)
        for k,v in acqp.iteritems():
            self.__dict__[k] = v
        self._yet_loaded = True
        
class ACQP_file(JCAMP_file):
    '''
    Thin specialization of JCAMP_file. Mostly just a place holder and
    a place to store the fact that these files are called "acqp".
    '''
    def __init__(self, dirname):
        super(ACQP_file, self).__init__(os.path.join(dirname,'acqp'))

class METHOD_file(JCAMP_file):
    '''
    Thin specialization of JCAMP_file. Mostly just a place holder and
    a place to store the fact that these files are called "method".
    '''
    def __init__(self, dirname):
        super(METHOD_file, self).__init__(os.path.join(dirname,'method'))
        
class BRUKER_fid(object):
    '''
    Initialize an fid object whih should sit in the scan root directory.
    
    __init__(filename) expects a directory name. It will look
    for the fid file, and the acqp and method parameter files.
    '''
    def __init__(self, dirname, lazy=True):
        self.dirname = dirname
        self._yet_loaded = False
        if not lazy:
            self._forced_load()
            
    def __getattr__(self, attr):
        '''
        intercepts only undefined attributes. this will be used to 
        trigger the loading of the fid data
        '''
        if not self._yet_loaded:
            self._forced_load() 
        return object.__getattribute__(self, attr)
        
    # TODO: currently, header files acqp and method are loaded twice.
    #       need to remove the redundancy. Somehow provide the headers
    #       if possible.
    def _forced_load(self):
        logger.info('loading %s' % self.dirname)
        fid = BRUKERIO.readfid(os.path.join(self.dirname,'fid'))
        self.__dict__['fid'] = fid['data']
        self._yet_loaded = True
        

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
    def __init__(self, filename):
        pdata = BRUKERIO.read2dseq(filename)
        self.__dict__['data'] = pdata['data']
        for k,v in pdata['header'].iteritems():
            self.__dict__[k] = dict2obj(v)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
#    import numpy
#    fname = os.path.expanduser('~/data/readfidTest.ix1/3')
#    x = BRUKER_fid(fname)
#    print('mean = {0}'.format(numpy.mean(x.fid)))
#    print(x.kasimir)