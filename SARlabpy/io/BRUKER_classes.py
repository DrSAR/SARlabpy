# -*- coding: utf-8 -*-
"""
Class Definitions for BRUKER data
"""
import os
import glob
import SARlabpy.io.BRUKERIO

import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

class Scan(object):
    '''
    Object to represent a BRUKER scan

    The __init__() method expects the filename to be a directory
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
        for f in glob.glob(os.path.join(dirname,'pdata','*','2dseq')):
            try:
                self.pdata.append(PDATA_2dseq_file(f))
            except IOError:
                pass

        if not any([self.acqp, self.fid, self.method]):
            #we have not found anything
            raise IOError(
                ('Directory "{0}" did not contain any of the typical BRUKER '+
                'files').format(dirname))
                            
class JCAMP_file(object):
    def __init__(self, filename):
        if not os.path.isfile(filename):
            raise IOError('File "%s" not found' % filename)
        # no go and read the file
        acqp = SARlabpy.io.BRUKERIO.readJCAMP(filename)
        for k,v in acqp.iteritems():
            self.__dict__[k]=v

class ACQP_file(JCAMP_file):
    def __init__(self, dirname):
        super(ACQP_file, self).__init__(os.path.join(dirname,'acqp'))

class METHOD_file(JCAMP_file):
    def __init__(self, dirname):
        super(METHOD_file, self).__init__(os.path.join(dirname,'method'))
        
class BRUKER_fid(object):
    def __init__(self, dirname):
        raise IOError('not yet implemented')

class PDATA_2dseq_file(object):
    def __init__(self, filename):
        raise IOError('not yet implemented')