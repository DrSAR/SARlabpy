# -*- coding: utf-8 -*-
"""
Test Routines for BRUKERIO
"""
import unittest
import doctest
import os

#import matplotlib.pyplot as plt
from SARlabpy.io import BRUKERIO
import logging
from pylab import imshow
logger = logging.getLogger('root')
logger.setLevel(logging.DEBUG)


class Test(unittest.TestCase):
    '''
    Unit tests for readfid.
    '''
    
    def test_doctests(self):
        '''Run the doctests'''
        logger.info('Running the doctests')
        doctest.testmod(BRUKERIO)
        
    def test_readfid(self):
        '''BRUKERIO.io.readfid()'''
        logger.info('test_readfid')
        fname = os.path.expanduser('~/datalink/readfidTest.ix1/3/fid')   
        data = BRUKERIO.readfid(fname)
        filesize = os.stat(fname).st_size
        kspace = data['data']
        logger.info('file size is {0}, and shape {1} '.format(filesize, kspace.shape))
        for i in range(kspace.shape[2]):
            arr = abs(kspace[:,:,i,0])
            imshow('a', arr / arr.max()*10)
            
        self.assertEqual(kspace.shape, (128,64,11,1))

def stresstest_readfid(skip=None):
    '''
    test.test_BRUKERIO.stresstest_readfid(skip=[3,5,8,9,12,13,15,16,17,18])
    '''
    skip = skip or [] # default if sentry is None
    logger.setLevel(logging.INFO)
    fnameroot = os.path.expanduser('~/datalink/readfidTest.ix1/')
    scan_labels = {1:'Tripilot multi', 
                   2:'FLASH 2D', 
                   3:'FLASH 3D', 
                   4:'MSME 2D', 
                   5:'MSME 3D', 
                   6:'MSME 2D-TURBO', 
                   7:'FLASH 2D (NR=25)', 
                   8:'FLASH 2D (NR < 25) partial acq. NR auto reset to 5', 
                   9:'(Scan 9) FLASH 2D (NR < 25) transferred prematurely', 
                   99:'(Scan 99) FLASH 2D (NR < 25) transferred prematurely', 
                   10:'FLASH 2D (MATRIX 32 X 32)', 
                   11:'FLASH 3D (MATRIX 32 X 32)', 
                   12:'EPI "1 segment"', 
                   13:'EPI "16 segments"', 
                   14:'DTI STANDARD', 
                   15:'DTI SPIRAL (did not show)', 
                   16:'UTE 2D',
                   17:'UTE 3D',
                   18:'ZTE 3D'}

    for k in [keys for keys in scan_labels.keys() if keys not in skip]:
        logger.info('opening scan {0} which is "{1}"'.
                     format(k,scan_labels[k]))
        fname = fnameroot+str(k)+'/fid'

        data = BRUKERIO.readfid(fname)
        filesize = os.stat(fname).st_size
        kspace = data['data']
        acqp = data['header']['acqp']
        logger.debug('filesize={0}\ndata shape={1}\nACQ_size = {2}'.
                     format(filesize, kspace.shape,acqp['ACQ_size']))
#        for i in range(kspace.shape[2]):
#            arr = abs(kspace[:,:,i,0])
#            imshow(arr / arr.max()*10)
#            x=input('press key')
    
if __name__ == "__main__":
    unittest.main()
