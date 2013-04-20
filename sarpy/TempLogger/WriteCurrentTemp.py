# -*- coding: utf-8 -*-
"""
Copyright: SARlab members, UBC, Vancouver, 2013
"""

# WriteCurrentTemp.py performs the following operations:
# 1) Open T2.json (written by cronjob)
# 2) Extract last temperature reading (last line)
# 3) Compares last read in time (over ten minues -> email) (still need to code)
# 4) Writes last temp to currenttemp.txt (read by html -> updates website)

import sys
import os
from pyparsing import (Word, Group, Optional, SkipTo, alphanums, nums,
                       StringEnd, delimitedList, OnlyOnce)
import datetime

if len(sys.argv) < 2:
    fname = 'T2.json'
else:
    fname = sys.argv[1]


# Extract last line from T2.json
with open(fname, 'rb') as fh:   # JSON = JavaScript Object Notation
    #first = next(fh)
    offset = -100 # 100 characters
    
    while True:
        fh.seek(offset, 2)          # lines = last 100 characters of file
        lines = fh.readlines()     
        #print lines
        #print len(lines)
        if len(lines)>1:        # len = length of line function
            #print 'In if statement'
            last = lines[-1]    # last = last line of "lines"
            #print lines[-3]
            break
        # Increase offset b/c last line is longer than 100 characters
        offset = offset*2
    #print first    
    #print last



# PARSING: Extract temperature from last (last line)
print last
# 2013-04-04 23:59:55 ch1(V/3.3V,Temp)=724.0,22.4 ch2=248.0,-25.2. Temp norm, Wixel Alive, Server Up. (8)
#["2013-04-12 14:29:33",-8.1]]]     # Our line to parse (extract -8.1)

date = Word(nums+'- :')
#time = Word(nums+':')
float_word = delimitedList(Word('-' + nums + '.') | "'error'") 
EOL = SkipTo(StringEnd())       

log_line = ('["' + date + '",' + float_word + EOL)

parsed_line = log_line.parseString(last)
print parsed_line
new_line={
    # 'date': parsed_line[0], 
    # 'time': parsed_line[1], 
    'T2': parsed_line[3] }
    # 'T2' : parsed_line[7].strip('.'),
    # 'EOL': parsed_line[8].strip()}
#time = datetime.datetime.strptime('2013-04-04 23:59:55',
#                                     '%Y-%m-%d %H:%M:%S')

#print( '[[{T1}]]'.format(**new_line) )
T2var='[[{T2}]]'.format(**new_line) # Extracted T2 is now a variable

# Calculate time since last ReadIn
dt = (datetime.datetime.now() - datetime.datetime.strptime(parsed_line[1],'%Y-%m-%d %H:%M:%S') 
      ).total_seconds()
print T2var, dt


# WRITE TO 'currenttemp.txt'
fo = open(fname+'.current', 'w')    # JSON = JavaScript Object Notation  #Opens a file for writing only. Overwrites the file if the file exists. If the file does not exist, creates a new file for writing.
fo.write(str(T2var));     # \n makes new line
#fo.write(repr( x));
#fo.write(repr( T));
fo.close() # close file









