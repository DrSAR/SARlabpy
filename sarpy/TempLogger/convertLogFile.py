# -*- coding: utf-8 -*-
"""
Copyright: SARlab members, UBC, Vancouver, 2013
"""

import sys
import os
from pyparsing import (Word, Group, Optional, SkipTo, alphanums, nums,
                       StringEnd, delimitedList, OnlyOnce)
import datetime

# 2013-04-04 23:59:55 ch1(V/3.3V,Temp)=724.0,22.4 ch2=248.0,-25.2. Temp norm, Wixel Alive, Server Up. (8)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        fname = os.path.expanduser('~/Desktop/TemperatureLogFile')
    else:
        fname = sys.argv[1]


    date = Word(nums+'-')
    time = Word(nums+':')
    float_word = delimitedList(Word('-' + nums + '.') | "'error'") 
    EOL = SkipTo(StringEnd())       

    log_line = (date + time + 'ch1(V/3.3V,Temp)=' + float_word +
                'ch2=' + float_word + EOL)


    with open(fname) as f:
        for line in f:
            parsed_line = log_line.parseString(line)
            new_line={
                'date': parsed_line[0],
                'time': parsed_line[1],
                'T1' : parsed_line[4],
                'T2' : parsed_line[7].strip('.'),
                'EOL': parsed_line[8].strip()}
#            time = datetime.datetime.strptime('2013-04-04 23:59:55',
#                                              '%Y-%m-%d %H:%M:%S')
            print('{date} {time} T1={T1} T2={T2} {EOL}'.format(**new_line))
