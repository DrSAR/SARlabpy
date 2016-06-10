#!/usr/bin/env python

import sys
import sarpy.io.BRUKER_classes as sib
if __name__ == '__main__':
    scn = sib.Scan('PhantomOrientation.iY1/'+sys.argv[1])
    fname = sys.argv[2]
    print(scn.pdata[0].uid())
    scn.pdata[0].write2nii(sys.argv[2]+'write.nii.gz')
    scn.pdata[0].export2nii(sys.argv[2]+'export.nii.gz')
