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

if len(sys.argv) < 2:   # Command in terminal: python WriteCurrentTemp.py T1.json
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



# PARSING: Extract temperature from (last line)
# print last
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



# SEND email w msg ('No update T2.json file in over 10 mins  ->  Check wixel')
EmailCount=0;
EmailThreshold=1;

if dt > 60*10: # Send email if dt (time diff b/w now and last update) > 10mins
    EmailCount=EmailCount+1;

    if EmailCount > EmailThreshold:     # Send an email if we are above the Threshold
        EmailThreshold = EmailThreshold*2 # this should lead to message every 1min, 2min, 4min, 8min, 16min etc after disaster...
        #logger.error('Email Actually sent!')
        # Open .txt file for sending latest temperature values  
        #fp = open(LogFileName, 'rb')
        # Create a text/plain message
        #msg = MIMEText(fp.read())
        msg = 'No update to the T2.json file in over 10 mins  ->  Check wixel'
        #fp.close()

        me = 'ubcmri7t@gmail.com'
        #you = ['stefan@phas.ubc.ca'] #,'stefan@phas.ubc.ca','clayton.wong1x@gmail.com']   # Can send to 2 recipients
        you = ['clw9@sfu.ca']
        #msg = MIMEMultipart('alternative')
        msg['Subject'] = 'No update to T2.json'  # "Hello %s, my name is %s" % ('john', 'mike') # Hello john, my name is mike".
        msg['From'] = me
        msg['To'] = ", ".join(you)
        #text = 'The temperature of the TMP36 in the 7T Freezer is hot (greater than 23 degCelius)!'



        #part1 = MIMEText(text, 'plain')
        #msg.attach(part1)
        #msg.attach(part2)                

        # Credentials
        username = 'ubcmri7t@gmail.com'  
        password = 'xwinnmr123'  

        # Send the message via our own SMTP server, but don't include the envelope header
        s = smtplib.SMTP('smtp.gmail.com:587')
        s.starttls() # SMTP connection in TLS (Transport Layer Security) mode. All SMTP commands that follow will be encrypted 
        s.login(username,password) 
        s.sendmail(me, you, msg.as_string())
        s.quit()



# WRITE TO 'T2.json.current'
fo = open(fname+'.current', 'w')    # JSON = JavaScript Object Notation  #Opens a file for writing only. Overwrites the file if the file exists. If the file does not exist, creates a new file for writing.
fo.write(str(T2var));     # \n makes new line
#fo.write(repr( x));
#fo.write(repr( T));
fo.close() # close file


# html then reads 'T2.json.current' and updates the gauge 









