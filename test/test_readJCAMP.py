from sarpy import (readJCAMP, dataroot, adataroot)
import os

x=readJCAMP('JCAMPtest.dat')
print(type(x['ExcPulse']))
print(x['ExcPulse'])

print(type(x['TPQQ']))
print(x['TPQQ'])