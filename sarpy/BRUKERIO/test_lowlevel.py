# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 22:26:40 2016

Test routines of lowlevel BRUKERIO functions
"""
import pytest
from . import lowlevel

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
    doctest.testmod(m=lowlevel)
    failure_count, test_count = doctest.testmod(m=lowlevel)
    assert failure_count == 0
    assert test_count == 60
