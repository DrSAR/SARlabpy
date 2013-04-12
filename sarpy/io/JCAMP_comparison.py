# -*- coding: utf-8 -*-
"""
Copyright: SARlab members, UBC, Vancouver, 2013
"""
from __future__ import print_function

import re
import numpy
import sarpy

def regex_comp(a,b): return re.match(a,b)
def arr_comp(a,b): return all(a==b)
def int_comp(a,b): return a==b

dictionary = {}




def list_param_by_type(type_=str, exp_names=['NecS3'], fname = None):
    '''
    create part of the soure code in this module:
    sort JCAMP parameter names into types (str, array and other) and create a
    dictionary that associate frozensets of these parameters with the
    corresponding compare function
    '''

    set_dict = {}
    for typ in [float, numpy.ndarray, str, list, int]:
        set_dict[typ]=set([])

    exp_all = None
    for x in exp_names:
        exp_i = sarpy.Experiment(x)
        if exp_all is None:
            exp_all = exp_i
        else:
            for stdy in exp_i.studies:
                exp_all.add_study(study=stdy)

    for stdy in exp_all.studies:
        for scn in stdy.scans:
            print(scn)
            for t in scn.acqp.__dict__:
                set_dict[type(scn.acqp[t])].add(t)
            for t in scn.method.__dict__:
                set_dict[type(scn.method[t])].add(t)

    if fname:
        with open(fname, 'w') as f:
            print('{frozenset([', file=f)
            for a in sarpy.natural_sort(set_dict[str]): #natural_sort
                print("'%s', " % a, file=f, end='')
            print("]): regex_comp,", file=f)
            print('frozenset([', file=f)
            for a in sarpy.natural_sort(set_dict[numpy.ndarray]):
                print("'%s', " % a, file=f, end='')
            print("]): arr_comp,", file=f)
            print("frozenset(['default']): int_comp}", file=f)
    else:
        return {frozenset(set_dict[str]): regex_comp,
                frozenset(set_dict[numpy.ndarray]): arr_comp,
                frozenset(['default']): int_comp,
               }

if __name__ == "__main__":
    list_param_by_type(exp_names=['readfid', 'NecS3Hs15', 'DiLL', 'HPG'],
                             fname = 'JCAMP_dict.py')