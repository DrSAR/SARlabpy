# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 10:40:00 2013

@author: stefan
"""

import os
import BRUKERIO
import pprint

files =  ['~/bdata/stefan/nmr/readfidTest.ix1/%i' % (i+1) for i 
         in range(18)+[98]]
files = ['~/bdata/stefan/nmr/NecS3Hs02.iJ1/10']         

for fid_file in files:
    print fid_file+'\n'+'-'*50
    try:
        fid = BRUKERIO.readfid(os.path.join(os.path.expanduser(fid_file),'fid'), untouched=True)
    except Exception, ex:
        print 'EXCEPTION'
        print ex
        acqp = BRUKERIO.readJCAMP(os.path.join(os.path.expanduser(fid_file),'acqp'))
    else:
        print 'no EXCEPTION'
        print fid['data'].shape
        acqp = fid['header']['acqp']
            
    pprint.pprint(['%s: %s' % (k,acqp.get(k, 'MISSING')) for 
           k in ['ACQ_dim','ACQ_size','GO_block_size','NI','NSLICES','NR','NREC','NSCANS',
                 'ACQ_obj_order','ACQ_phase_factor','ACQ_ns_list_size','ACQ_ns_list',
                 'ACQ_phase_encoding_mode',
                 'ACQ_phase_enc_start','ACQ_rare_factor',
                 'ACQ_spatial_phase_1','ACQ_spatial_size_1']])

    print '\n'+'='*50