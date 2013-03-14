# -*- coding: utf-8 -*-
"""
========================================
Quick check for the presence of all sorts of BRUKER files in the 
entire dataroot. A fully stocked scan should have
acqp, method, fid, 2dseq. This results in a score of
1 +   2 +     4 +  8 = 15. If one is missing, the score is less.
========================================
"""
import os
import glob
from BRUKER_classes import dataroot

print(__doc__)

studies = os.listdir(dataroot)
counter = 0
for study in studies:
    scans = [x for x in os.listdir(os.path.join(dataroot, study))
             if os.path.isdir(os.path.join(dataroot,study,x))]
    for scan in scans:
        if not (scan == 'FieldMap' or scan == 'AdjResult'):
            scanname = os.path.join(dataroot,study,scan)
            testval = 1*os.path.isfile(scanname+'/acqp') + \
                      2*os.path.isfile(scanname+'/method') + \
                      4*os.path.isfile(scanname+'/fid') + \
                      8*(len(glob.glob(scanname+'/pdata/*/2dseq'))>0)
            
            if testval !=15:
                print('{0} / {1} : {2}'.format(study, scan, testval))
            counter += 1
        
print('-'*40 + '\n{0} scans tested'.format(counter))