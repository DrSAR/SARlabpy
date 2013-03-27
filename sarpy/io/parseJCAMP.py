# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 21:10:04 2013

@author: stefan
"""

import re, os
from pyparsing import (Word, Group, Optional, SkipTo, alphanums, nums,
                       StringEnd, delimitedList, Forward)

from BRUKERIO import pairwise
from BRUKER_classes import dataroot

def parse_JCAMP(string_list):             
    
    EOL = SkipTo(StringEnd())       
    integer = Word(nums)
    shape = '(' + delimitedList(integer)('shape') + ')'
    paramname=Word(alphanums+"_")('pn')
    paramvalue=Word(alphanums+"_"+'.'+'<'+'>'+'('+')' +'+'+'-')
    param = ('##' + Optional('$') + paramname('paramname') + '=' 
            + Optional(shape) + paramvalue('paramvalue'))
    
    comment = "$$ @vis" + EOL
    
    date = "$$ " + Word(alphanums)('datestr') + EOL
    location = "$$ " + Word(alphanums+'/'+'.') + EOL
    
    param_line = Group(param | comment | date| location)

    LDRdict = {}
    for line in string_list:
#        print line
        parsed_data = param_line.parseString(line)
        if parsed_data[0].paramname:
#            print('{0} == ## {2} ## {1}'.format(parsed_data[0].paramname,
#                          parsed_data[0].paramvalue, parsed_data[0].shape))
            LDRdict[parsed_data[0].paramname] = parsed_data[0].paramvalue
    return LDRdict
        
def join_JCAMP_lines(fname):
    with open(fname, "r") as JCAMPfile:
        JCAMPdata = JCAMPfile.read().splitlines() # read and lose the "\n"
    LDRlist = [] # start with empty list
    currentline = ''
    JCAMPdata_iterator = pairwise(JCAMPdata)
    for line, next_line in JCAMPdata_iterator:
        currentline += line + ' '
        if re.match('\\$\\$|##', next_line):
            LDRlist.append(currentline)
            currentline = ''
            
    return LDRlist
            
    

if __name__ == '__main__':
    fname = os.path.join(dataroot, 'stefan','nmr',
                         'readfidTest.ix1', '2', 'acqp')
    joined = join_JCAMP_lines(fname)
    parse_JCAMP(joined)
    