# -*- coding: utf-8 -*-
"""
Copyright: SARlab members, UBC, Vancouver, 2013
"""
import re
list_of_dicts = [{'s1':'abcd', 's2':'ABC', 'i':42, 'f':4.2},
                 {'s2':'xyz', 'i':84, 'f':8.4, 'notmentioned': 10}]

 # just making up some comparison strategies
def regex_comp(a,b): return re.match(a,b)
def int_comp(a,b): return a==b
def float_comp(a,b): return round(a,-1) == round (b,-1)

pre_specified_comp_dict = {frozenset(['s1','s2']) : regex_comp,
                           frozenset(['i']): int_comp,
                           frozenset(['f']): float_comp,
                           frozenset(['default']): int_comp,}

def fle_new(**kwargs):
    chosen_comps={}
    for key in kwargs.keys():
        # remember, the keys here are frozensets
        cand_comp = [frozenset(['default'])] or [x for x in
                                           pre_specified_comp_dict if key in x]
        chosen_comps[key] = pre_specified_comp_dict[cand_comp[0]]

    matches = lambda d: all(k in d and chosen_comps[k](v, d[k])
                            for k, v in kwargs.items())

    return filter(matches, list_of_dicts)

def fle_vivekpoddar(**kwargs):
    for data in list_of_dicts:
        for key, val in kwargs.iteritems():
            if not data.get(key, False) == val:
                return []
            else:
                return data


def fle_orig(**kwargs):
    for s in list_of_dicts:
        for criterion, criterion_val in kwargs.iteritems():
            if type(criterion_val) is str:
                if re.match(criterion_val, s.get(criterion, 'unlikely_return_val')):
                    yield s
                    continue
            if s.get(criterion, None) == criterion_val:
                yield s


def fle_srgerg(**kwargs):
    compare = lambda e, a: re.match(e, a) is not None if isinstance(e, basestring) else e == a
    matches = lambda d: all(k in d and compare(v, d[k]) for k, v in kwargs.iteritems())
    return filter(matches, list_of_dicts)

for method in [fle_orig,
               fle_new,
#               fle_vivekpoddar,
               fle_srgerg]:
    print '-'*49
    print method.__name__
    print '-'*49

    print 'no key        {0}'.format(list(method()))
    print 'key missing   {0}'.format(list(method(i=41)))
    print 'key there     {0}'.format(list(method(i=42)))
    print 'string search {0}'.format(list(method(s2='xyz')))
    print 'regex search  {0}'.format(list(method(s2='[a-z]')))
    print 'keys combined {0}'.format(list(method(i=42, f=4.2)))
    print 'key not known {0}'.format(list(method(notmentioned=1)))
    print 'key not known {0}'.format(list(method(notmentioned=10)))
