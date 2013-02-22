# -*- coding: utf-8 -*-
"""
Test Routines for BRUKERIO
"""
import unittest
import doctest
import os
import time

import cv2
import matplotlib.pyplot as plt
import SARlabpy
from SARlabpy.io import BRUKERIO
from SARlabpy.io.BRUKERIO import readfid


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
        fname = os.path.expanduser('~/datalink/NecS1Hs08.hk1/3')
        data = readfid(fname+'/fid')
        kspace = data['data']
        print (kspace.shape)
        for i in range(kspace.shape[2]):
            arr = abs(kspace[:,:,i,0])
            cv2.imshow('a', arr / arr.max()*10)
            cv2.waitKey(0)
        cv2.destroyAllWindows
            
        self.assertEqual(kspace.shape, (128,128,3,1))
        
if __name__ == "__main__":
    unittest.main()
