#! /bin/sh
""":"
exec python $0 ${1+"$@"}
"""

import sys

while True:
    value = sys.stdin.readline()
    print 'something is coming'
    if value == 'other wixel dead':
        print('it is dead')