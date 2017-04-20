# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 22:26:40 2016

Test routines of lowlevel BRUKERIO functions
"""
import pytest
import numpy
from . import lowlevel

#for the purposes of testing, point dataroot to testmodule
import os
curdir = os.path.dirname(__file__)
# was: dataroot = os.path.expanduser(os.path.join('~','bdata'))
lowlevel.dataroot = os.path.join(curdir,'testdata')
README_file_inside_submodule = os.path.join(lowlevel.dataroot,'README.md')
assert os.path.exists(README_file_inside_submodule),'''
For testing, the submodule gitlab.com:DrSAR/sarpy-testdata.git must be 
present and pulled down into the sarpy/BRUKERIO/testdata directory
On initial cloning of the repo you can achieve this by recursive cloning:
git clone --recursive git@pfeifer.phas.ubc.ca:SARlabpy.git
Later, you can still issue inside the repo:
git submodule update --init --recursive
'''

fnameroot = os.path.join(lowlevel.dataroot,'stefan','nmr','readfidTest.ix1')

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

PDATA_size_dict = {1: [256, 256, 15]  ,
                2: [105, 133, 5],
                3: [105, 133, 25],
                4: [105, 133, 5],
                5: [105, 133, 25],
                6: [256, 256, 3, 5],
                7: [105, 133, 5, 25],
                8: [105, 133, 5, 5],
                10: [32, 32, 5],
                11: [32, 32, 5],
                12: [64, 64, 5],
                13: [256, 256, 3],
                14: [105, 256, 5, 2],
                16: [128, 128, 5],
                17: [128, 128, 128],
                18: [128, 128, 128]}
                
class Test_minors:
    def test_pairwise(self):
        it = [1, 2, 3]
        pw = lowlevel.pairwise(it)
        assert  next(pw) == (1,2)
        assert  next(pw) == (2,3)
        with pytest.raises(StopIteration):
            next(pw)

    def test_convert_int_float_string(self):
        assert lowlevel.convert_int_float_string('13')==13
        assert lowlevel.convert_int_float_string('13.3')==13.3
        assert lowlevel.convert_int_float_string('x1')=='x1'
        
def test_doctests():
    '''Run the doctests'''
    import doctest
    failure_count, test_count = doctest.testmod(m=lowlevel)
    assert failure_count == 0
    assert test_count == 60

def test_read2dseq():
    skip = [9, 15, 17, 18, 99]
    for k in [keys for keys in list(scan_labels.keys()) if keys not in skip]:
        fname_pdata = os.path.join(fnameroot, str(k), 'pdata','1')
        d = lowlevel.read2dseq(fname_pdata)
        print(d['data'].shape, PDATA_size_dict[k])
        assert len(d['dimdesc']) == len(d['dimcomment'])
        assert PDATA_size_dict[k] == list(d['data'].shape)

