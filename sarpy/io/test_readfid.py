# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 10:40:00 2013

@author: stefan
"""

import os
import BRUKERIO
import pprint
#import sarpy

root = '~/bdata/stefan/nmr/'
files =  ['readfidTest.ix1/%i' % (i+1) for i 
         in range(18)+[98]]
files.append('NecS3Hs02.iJ1/10')
files = ['readfidTest.ix1/8']
files = ['NecS3Hs02.iJ1/10']
for fid_file in files:
        print fid_file+' '+'-'*50
#    try:
        fid = BRUKERIO.readfid(os.path.join(os.path.expanduser(root),fid_file,'fid'), untouched=False)
        print('fid-shape = %s' % (str(fid['data'].shape)))
#    except Exception, ex:
        print 'EXCEPTION ' + fid_file
        print ex
#    try:
#        scn=sarpy.Scan(fid_file)
#        print('2dseq-shape = %s' % (str(scn.pdata[0].data.shape)))
#    except Exception, ex:
#        print 'EXCEPTION reading scan'
#        print ex
        acqp = BRUKERIO.readJCAMP(os.path.join(os.path.expanduser(root),fid_file,'acqp'))        
#    pprint.pprint(['%s: %s' % (k,str(acqp.get(k, 'MISSING'))[0:100]) for 
#        k in [
#        'ACQ_protocol_name',
#        'ACQ_dim','ACQ_size','GO_block_size','NI','NSLICES','NR','NREC','NSCANS',
#        'ACQ_obj_order','ACQ_phase_factor','ACQ_ns_list_size','ACQ_n_movie_frames',
#         'ACQ_ns_list',
#        'ACQ_phase_encoding_mode',
#        'ACQ_phase_enc_start','ACQ_rare_factor',
#        'ACQ_spatial_phase_1',
#        'ACQ_spatial_size_1']])

#    print '\n'+'='*50