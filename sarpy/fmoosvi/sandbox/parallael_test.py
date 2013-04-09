# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 13:09:58 2013

@author: fmoosvi
"""
import numpy
import time

from IPython.parallel import Client

def foo(azimuth, zenith):
    # Do various bits of stuff
    # Eventually get a result
    
    for i in range(1000):
        for j in range(1000):
            pass
    return 5e8

azimuths = numpy.linspace(0,100,50)
zeniths = numpy.linspace(1000,0,50)

####### Without Parallelization

tasks1 = []
st2 = time.time()
for azimuth in azimuths:
    for zenith in zeniths:
        tasks1.append(foo(azimuth, zenith))

et2 = time.time()
print  'Without parallelization : {0} s'.format(et2 - st2)


####### With Parallelization
c = Client()   # here is where the client establishes the connection
lv = c.load_balanced_view()   # this object represents the engines (workers)

tasks = []

st1 = time.time()
for azimuth in azimuths:
    for zenith in zeniths:
        tasks.append(lv.apply(foo, azimuth, zenith))

et1 = time.time()

print 'With parallelization : {0} s'.format(et1 - st1)

