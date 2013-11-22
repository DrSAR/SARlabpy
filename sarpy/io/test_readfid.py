# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 10:40:00 2013

@author: stefan
"""

import os
import BRUKERIO
import pprint

files =  ['~/bdata/stefan/nmr/readfidTest.ix1/%i/fid' % (i+1) for i 
         in range(18)+[98]]
files = ['~/bdata/stefan/nmr/NecS3Hs02.iJ1/10/fid']         

for fid_file in files:
    print fid_file+'\n'+'-'*50
    try:
        fid = BRUKERIO.readfid(os.path.expanduser(fid_file), untouched=False)
    except Exception, ex:
        print 'EXCEPTION'
        print ex
    else:
        print 'no EXCEPTION'
            
    print fid['data'].shape
    
    pprint.pprint(['%s: %s' % (k,(fid['header']['acqp']).get(k, 'MISSING')) for 
           k in ['ACQ_dim','ACQ_size','GO_block_size','NI','NSLICES','NR',
                 'ACQ_obj_order','ACQ_phase_factor',
                 'ACQ_phase_encoding_mode',
                 'ACQ_phase_enc_start','ACQ_rare_factor',
                 'ACQ_spatial_phase_1','ACQ_spatial_size_1']])

    print '\n'+'='*50