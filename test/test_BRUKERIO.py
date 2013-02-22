# -*- coding: utf-8 -*-
"""
Test Routines for BRUKERIO
"""
import unittest
import doctest
import os

from pylab import imshow
import SARlabpy
from SARlabpy.io.BRUKERIO import readfid


class Test(unittest.TestCase):
    '''
    Unit tests for readfid.
    '''
    
    def def_doctests(self):
        '''Run the doctests'''
        doctests.testmod(BRUKERIO)
    def test_readfid(self):
        '''BRUKERIO.io.readfid()'''
        fname = os.path.expanduser('~/datalink/Moosvi.ii2/3')
        data = readfid(fname)
        kspace = data['data']