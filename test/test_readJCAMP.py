from sarpy import (readJCAMP, dataroot, adataroot)
import os

x=readJCAMP('untracked_JCAMPtest.dat')
print(type(x ['ACQ_scan_size']))
print(x ['ACQ_scan_size'])

print(type(x ['ACQ_ns']))
print(x ['ACQ_ns'])

print(type(x ['PULPROG']))
print(x ['PULPROG'])

print(type(x['ExcPulse']))
print(x['ExcPulse'])

print(type(x['TPQQ']))
print(x['TPQQ'])
