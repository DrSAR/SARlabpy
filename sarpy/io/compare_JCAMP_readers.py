# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 22:16:02 2013

@author: stefan
"""

import os
import timeit
import cProfile
from BRUKER_classes import dataroot
from BRUKERIO import readJCAMP
from BRUKERIO_ay import readBrukerParx
from parseJCAMP import parse_JCAMP, join_JCAMP_lines


def sarcode(fname):
    from BRUKERIO import readJCAMP
    acqp =readJCAMP(fname)
    
def acycode(fname):
    from BRUKERIO_ay import readBrukerParx
    acqp =readBrukerParx(fname)
    
def parsecode(fname):
    #from parseJCAMP import parse_JCAMP, join_JCAMP_lines
    acqp =parse_JCAMP(join_JCAMP_lines(fname))

def show_dict_comparison(d1, d2):
    agreed = disagreed = 0
    for k,v in d1.iteritems():
        try:
            if v != d2[k]:
                print('{0}:\n{1}\n{2}'.format(
                        k,v,d2[k]))
                disagreed += 1
            else:
                agreed += 1
        except ValueError:
            disagreed += 1
            print('{0}: {1}  ----  {2}'.format(
                        k,v,d2[k])) 
            
    print('agreeing on {0} and disagreeing on {1} lines'.
            format(agreed, disagreed))

if __name__ == '__main__':
    fname = os.path.join(dataroot, 'stefan','nmr',
                         'readfidTest.ix1', '2', 'acqp')

    # time trials                         
    t1 = timeit.Timer('sarcode("'+fname+'")',
                     'from __main__ import sarcode')
    print('sarcode ran for %s s' % t1.timeit(100))
    t2 = timeit.Timer('acycode("'+fname+'")',
                     'from __main__ import acycode')
    print('acycode ran for %s s' % t2.timeit(100))
    t3 = timeit.Timer('acqp =parse_JCAMP(join_JCAMP_lines("'+fname+'"))',
                     'from parseJCAMP import parse_JCAMP, join_JCAMP_lines')
    print('pyparsing code ran for %s s' % t3.timeit(100))
    
    #profiler
    cProfile.run('sarcode("'+fname+'")')
    cProfile.run('acycode("'+fname+'")')
    cProfile.run('parsecode("'+fname+'")')
   
 
    print('='*80+'\nvalue comparison for acqp\n'+'-'*80)
    JCAMP_sar = readJCAMP(fname)
    JCAMP_ay = readBrukerParx(fname)
    show_dict_comparison(JCAMP_sar, JCAMP_ay)

            
    print('='*80+'\nvalue comparison for visu_pars\n'+'-'*80)
    fname = os.path.join(dataroot, 'stefan','nmr',
                         'readfidTest.ix1', '2', 'pdata', '1', 'visu_pars')
    JCAMP_sar = readJCAMP(fname)
    JCAMP_ay = readBrukerParx(fname)
    show_dict_comparison(JCAMP_sar, JCAMP_ay)
