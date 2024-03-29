# -*- coding: utf-8 -*-
"""
Test Routines for lowlevel
"""
import unittest
import doctest
import os

#import matplotlib.pyplot as plt
#import pylab
import numpy

#make the SARlogger features available ??
from .BRUKER_classes import dataroot
from . import lowlevel
#initiate_logging(lowlevel, handler_level=30)

print('test_lowlevel run: %s' % __name__)

fnameroot = os.path.join(dataroot,'stefan','nmr','readfidTest.ix1')

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
    
    def test_doctests(self):
        '''Run the doctests'''
        doctest.testmod(lowlevel)
        
    def test_readfid(self):
        '''lowlevel.io.readfid()'''
#        fname = os.path.expanduser('~/data/readfidTest.ix1/3/fid')   
#        data = lowlevel.readfid(fname)
#        kspace = data['data']
#        self.assertEqual(kspace.shape, (256,105,1,1))
#        # reshape array by stacking them together
#        nslices = kspace.shape[2]
#        kspace_stack = abs(numpy.hstack(
#                        kspace[:,:,i,0] for i in range(0, nslices)))
#        pylab.imshow(numpy.log(kspace_stack+100))
#        pylab.show()
        
    def test_readJCAMP(self):
        '''
        Test readJCAMP for exceptions by reading variety of acqp 
        and method files
        '''
#        logger.info('test_readJCAMP started')
        skip = [] # default 
     
        for k in [keys for keys in list(scan_labels.keys()) if keys not in skip]:
#            lowlevel.logger.info('opening scan {0} which is "{1}"'.
#                         format(k,scan_labels[k]))
            fname_acqp = fnameroot+str(k)+'/acqp'
            fname_method = fnameroot+str(k)+'/method'
    
            acqp = lowlevel.readJCAMP(fname_acqp)
#            lowlevel.logger.info('ACQ_size = {0}'.format(acqp['ACQ_size']))
            self.assertEqual(acqp['ACQ_O2_list'], [0])
            ACQ_size = acqp['ACQ_size']
            self.assertEqual(ACQ_size, ACQ_size_dict[k])

            meth = lowlevel.readJCAMP(fname_method)
            try:
               PVM_EncMatrix = meth['PVM_EncMatrix']
               PVM_EncMatrix[0] = 2*PVM_EncMatrix[0]
               if k not in (12,13,16,17,18):
                   self.assertEqual(ACQ_size, PVM_EncMatrix)
            except KeyError:
                pass
#                lowlevel.logger.info('scan {0} has no PVM_EncMatrix'.format(k))

    def test_read2dseq(self):
        skip = [9, 15, 99]
        for k in [keys for keys in list(scan_labels.keys()) if keys not in skip]:
#            lowlevel.logger.warn('opening procno for scan {0} which is "{1}"'.
#                                format(k,scan_labels[k]))
            fname_pdata = fnameroot+str(k)+'/pdata/1'
            d = lowlevel.read2dseq(fname_pdata)
            print(d['data'].shape)
            self.assertEqual(len(d['dimdesc']), len(d['dimcomment']))
            

def stresstest_readJCAMP(*path_filter):
    import glob
    file_list = glob.glob(os.path.join(dataroot,'stefan',*path_filter))
    for JCAMPname in file_list:
        JCAMPfile=lowlevel.readJCAMP(JCAMPname)
        print(JCAMPname)

def stresstest_read2dseq(*path_filter):
    import glob
    file_list = glob.glob(os.path.join(dataroot,*path_filter))
    for fname in file_list:
        if os.path.isfile(fname+'/visu_pars'):
            print(fname, lowlevel.read2dseq(fname)['data'].shape)
        
        
        
def stresstest_readfid(skip=None):
    '''
    test.test_lowlevel.stresstest_readfid(skip=[9,13,15])
    '''
    skip = skip or [] # default if sentry is None
 
    for k in [keys for keys in list(scan_labels.keys()) if keys not in skip]:
#        lowlevel.logger.info('opening scan {0} which is "{1}"'.
#                     format(k,scan_labels[k]))
        fname = fnameroot+str(k)+'/fid'

        try:
            data = lowlevel.readfid(fname)
            filesize = os.stat(fname).st_size
        except IOError:
            print('File NOT found: '+fname)
        else:            
            kspace = data['data']
            imgspace = lowlevel.fftbruker(kspace, 
                                          encoding=data['header']['encoding'])
            acqp = data['header']['acqp']
#            lowlevel.logger.debug('filesize={0}\ndata shape={1}\nACQ_size = {2}'.
#                         format(filesize, kspace.shape,acqp['ACQ_size']))
                         
            # reshape array by stacking them together
            nslices = kspace.shape[2]
            kspace_stack = numpy.hstack(
                            kspace[:,:,i,0] for i in range(0, nslices))
            img_stack = numpy.hstack(
                            imgspace[:,:,i,0] for i in range(0, nslices))
            panel = numpy.vstack((kspace_stack, img_stack))
            pylab.imshow(numpy.log(abs(panel)+abs(panel).max()/1000))
            pylab.show()

if __name__ == "__main__":
#    unittest.main()
#    stresstest_readfid(skip=[6,9,13,14,15])
#    stresstest_readJCAMP('nmr','*','subject')
#    stresstest_readJCAMP('nmr','*','*','acqp')
#    stresstest_readJCAMP('nmr','*','*','method')
    stresstest_read2dseq('stefan','nmr','*','*','pdata','1')