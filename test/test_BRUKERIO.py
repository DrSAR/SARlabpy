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
import pylab
import numpy
logger = logging.getLogger('root')
logger.setLevel(logging.DEBUG)

print('test_BRUKERIO run: %s' % __name__)

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
            pylab.imshow('a', arr / arr.max()*10)
            
        self.assertEqual(kspace.shape, (128,64,11,1))

fnameroot = os.path.expanduser('~/datalink/readfidTest.ix1/')
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

def stresstest_readfid(skip=None):
    '''
    test.test_BRUKERIO.stresstest_readfid(skip=[9,13,15])
    '''
    skip = skip or [] # default if sentry is None
    logger.setLevel(logging.INFO)
 
    for k in [keys for keys in scan_labels.keys() if keys not in skip]:
        logger.info('opening scan {0} which is "{1}"'.
                     format(k,scan_labels[k]))
        fname = fnameroot+str(k)+'/fid'

        data = BRUKERIO.readfid(fname)
        filesize = os.stat(fname).st_size
        kspace = data['data']
        imgspace = BRUKERIO.fftbruker(kspace, 
                                      encoding=data['header']['encoding'])
        acqp = data['header']['acqp']
        logger.debug('filesize={0}\ndata shape={1}\nACQ_size = {2}'.
                     format(filesize, kspace.shape,acqp['ACQ_size']))
                     
        # reshape array by stacking them together
        nslices = kspace.shape[2]
        kspace_stack = numpy.hstack(
                        kspace[:,:,i,0] for i in range(0, nslices))
        img_stack = numpy.hstack(
                        imgspace[:,:,i,0] for i in range(0, nslices))
        panel = numpy.vstack((kspace_stack, img_stack))
        pylab.imshow(numpy.log(abs(panel)+abs(panel).max()/1000))
        pylab.show()

def stresstest_readJCAMP(skip=None):
    '''
    Test readJCAMP for exceptions by reading variety of acqp 
    and method files
    '''
    skip = skip or [] # default if sentry is None
    logger.setLevel(logging.INFO)
 
    for k in [keys for keys in scan_labels.keys() if keys not in skip]:
        logger.info('opening scan {0} which is "{1}"'.
                     format(k,scan_labels[k]))
        fname_acqp = fnameroot+str(k)+'/acqp'
        fname_method = fnameroot+str(k)+'/method'

        hdr = BRUKERIO.readJCAMP(fname_acqp)
        try:
#            print('ACQ_obj_order = {0}'.format(hdr['ACQ_obj_order']))
#            print('ACQ_phase_factor = {0}'.format(hdr['ACQ_phase_factor']))
#            print('ACQ_phase_encoding_mode = {0}'.
#                                format(hdr['ACQ_phase_encoding_mode']))
#            print('ACQ_phase_enc_start = {0}'.
#                                format(hdr['ACQ_phase_enc_start']))
#            print('ACQ_rare_factor = {0}'.format(hdr['ACQ_rare_factor']))
            print('ACQ_spatial_phase_1 = {0}'.
                    format((hdr['ACQ_spatial_phase_1'])))
        except KeyError:
            print('key not found!')
        hdr = BRUKERIO.readJCAMP(fname_method)
    
if __name__ == "__main__":
#    unittest.main()
    stresstest_readfid(skip=[6,9,13,14,15])
