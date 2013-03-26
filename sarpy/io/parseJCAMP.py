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
    param_line = Forward()
    EOL = SkipTo(StringEnd())       
    integer = Word(nums)
    shape = '(' + delimitedList(integer) + ')'
    paramname=Word(alphanums+"_")('pn')
    paramvalue=Word(alphanums+"_"+'.'+'<'+'>'+'('+')' +'+'+'-')
    param = '##' + Optional('$') + paramname + '=' + Optional(shape) + paramvalue
    
    comment = "$$ @vis" + EOL
    
    date = "$$ " + Word(alphanums)('datestr') + EOL
    location = "$$ " + Word(alphanums+'/'+'.') + EOL
    
    param_line << Group(param("pp") | comment | date| location)

    for line in string_list:
        print line
        parsed_data = param_line.parseString(line)
        print parsed_data.dump()
        print parsed_data.pn
        print parsed_data.datestr

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
    