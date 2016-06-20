# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 22:26:40 2016

Test routines of lowlevel BRUKERIO functions
"""
import pytest
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
