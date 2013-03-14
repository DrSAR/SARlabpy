from sarpy import (readJCAMP, dataroot, adataroot)
import os

x=readJCAMP(os.path.join(dataroot,'readfidTest.ix1','2','method'))
print(type(x['ExcPulse']))
print(x['ExcPulse'])

y=readJCAMP(os.path.join(dataroot,'readfidTest.ix1','2','acqp'))
print(type(y['TPQQ']))
print(y['TPQQ'])