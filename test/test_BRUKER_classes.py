# -*- coding: utf-8 -*-
"""
Test Routines for BRUKERIO
"""
import unittest
import doctest
import sys
import os
import glob

#make the SARlogger features available
from SARlabpy import SARlogger
#test logging for BRUKERIO
from SARlabpy.io import BRUKER_classes
SARlogger.initiate_logging(BRUKER_classes, handler_level=SARlogger.INFO)

print('testing: %s' % __name__)

fnameroot = os.path.expanduser('~/data/readfidTest.ix1/')
scan_labels = {1:'Tripilot multi', 
               2:'FLASH 2D', 
               3:'FLASH 3D', 
               4:'MSME 2D', 
               5:'MSME 3D', 
               6:'MSME 2D-TURBO', 
               7:'FLASH 2D (NR=25)', 
               8:'FLASH 2D (NR < 25) partial acq. NR auto reset to 5', 
               9:'(no fid)! FLASH 2D (NR < 25) transferred prematurely', 
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

ACQ_size_dict = {1: [256, 128],
                10: [64, 32],
                11: [64, 32],
                12: [8192, 1],
                13: [10094, 16],
                14: [266, 105],
                15: [1024],
                16: [128, 402],
                17: [128, 51360, 1],
                18: [1024, 51896, 1],
                2: [266, 105],
                3: [266, 105, 25],
                4: [266, 105],
                5: [266, 105, 25],
                6: [512, 256],
                7: [266, 105],
                8: [266, 105],
                9: [266, 105],
                99: [266, 105]}

class Test(unittest.TestCase):
    '''
    Unit tests for readfid.
    '''
    
    def doctests(self):
        '''Run the doctests'''
        doctest.testmod(BRUKER_classes)
        
    def test_bigslurp(self):
        files = glob.glob(os.path.join(os.path.expanduser('~/data'),'*/[0-9]*'))

        data=[]        
        for datadir in files:
            data.append(BRUKER_classes.Scan(datadir))
            print('loaded {0} (datasize = {1})'.
                  format(datadir, sys.getsizeof(data)))
            
        print('{0} scans loaded.'.format(len(files)))
        
if __name__ == "__main__":
    unittest.main()
#    stresstest_readfid(skip=[6,9,13,14,15])
