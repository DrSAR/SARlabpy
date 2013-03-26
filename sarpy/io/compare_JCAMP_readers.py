# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 22:16:02 2013

@author: stefan
"""

import os
from BRUKERIO import readJCAMP
from BRUKER_classes import dataroot
from BRUKERIO_ay import readBrukerParx


if __name__ == '__main__':
    fname = os.path.join(dataroot, 'stefan','nmr',
                         'readfidTest.ix1', '2', 'acqp')
    acqp_sar = readJCAMP(fname)
    acqp_ay =readBrukerParx(fname)

    agreed = disagreed = 0
    for k,v in acqp_sar.iteritems():
        if v != acqp_ay[k]:
            print('{0}: {1}  ====  {2}'.format(
                    k,v,acqp_ay[k]))
            disagreed +=1
        else:
            agreed +=1
            
    print agreed, disagreed