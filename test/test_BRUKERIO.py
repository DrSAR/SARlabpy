# -*- coding: utf-8 -*-
"""
Test Routines for BRUKERIO
"""
import unittest
import doctest
import os

import cv2
#import matplotlib.pyplot as plt
from SARlabpy.io import BRUKERIO
import logging
logger = logging.getLogger('root')
logger.setLevel(logging.DEBUG)


class Test(unittest.TestCase):
    '''
    Unit tests for readfid.
    '''
    
    def test_doctests(self):
        '''Run the doctests'''
        print('Running the doctests')
        doctest.testmod(BRUKERIO)
        
    def test_readfid(self):
        '''BRUKERIO.io.readfid()'''
        print('test_readfid')
        fname = os.path.expanduser('~/datalink/NecS1Hs08.hk1/3/fid')   
        data = BRUKERIO.readfid(fname)
        filesize = os.stat(fname).st_size
        kspace = data['data']
        print('file size is {0}, and shape {1} '.format(filesize, kspace.shape))
        print (kspace.shape)
        for i in range(kspace.shape[2]):
            arr = abs(kspace[:,:,i,0])
            cv2.imshow('a', arr / arr.max()*10)
            cv2.waitKey(0)
        cv2.destroyAllWindows
            
        self.assertEqual(kspace.shape, (128,64,11,1))

def test_readfid():
    fname = os.path.expanduser('~/datalink/readfidTest.ix1/2/fid')
    data = BRUKERIO.readfid(fname)
    filesize = os.stat(fname).st_size
    kspace = data['data']
    acqp = data['header']['acqp']
    logger.debug('filesize {0}\ndata shape: {1}'.format(filesize, kspace.shape))
    print('{0}'.format(acqp['ACQ_size']))
    for i in range(kspace.shape[2]):
        arr = abs(kspace[:,:,i,0])
        cv2.imshow('a', arr / arr.max()*10)
        cv2.waitKey(0)
        cv2.destroyAllWindows
    
if __name__ == "__main__":
    unittest.main()
