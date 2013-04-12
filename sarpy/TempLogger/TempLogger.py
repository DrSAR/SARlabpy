#!/usr/bin/env python
# Temperature Logger (for monitoring Brain Tumour Freezer Temperature)
# Clayton Wong
# April 5th, 2013

# Main loop calls the functions:
# 1) CheckServer
# 2) ReadSerial
# 3) Logger
# 4) Email

# STATE TABLE:
#----------------------------------------------------------------------
#             INPUTS             ||       OUTPUTS      || 
# Tnorm   WixelAlive   ServerUp  ||   Log      Email   || Case Number
#   x          0           0     || Tn,WA,SU     0     ||     1
#   x          0           1     || Tn,WA,SU     1     ||     2
#   0          1           0     || Tn,WA,SU     0     ||     3
#   0          1           1     || Tn,WA,SU     1     ||     4
#   x          0           0     || Tn,WA,SU     0     ||     5 (Same as 1)
#   x          0           1     || Tn,WA,SU     1     ||     6 (Same as 2)
#   1          1           0     || Tn,WA,SU     0     ||     7
#   1          1           1     || Tn,WA,SU     0     ||     8 (WANT THIS! Everything Working)

import serial # import serial for serial input
#import datetime # import datetime for getting localtime
import smtplib # import smtplib simple mail transfer protocol for EMAIL
# Import the email modules we'll need 
from email.mime.text import MIMEText
#from email.mime.multipart import MIMEMultipart

T1HIGH = 32  # room temp getting too hot
T2HIGH = -20 # freezer is getting too warm
LogFileName = '/var/log/SystemsHealth/TemperatureLogFile'

import logging, logging.handlers    # Import logger modules
logger = logging.getLogger('TemperatureLogFile') # create logger instance
hdlr = logging.handlers.TimedRotatingFileHandler(LogFileName, when="midnight")
formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')    # create formatter (Layout of log)

hdlr.setFormatter(formatter)    # attach formatter to handler
logger.addHandler(hdlr)         # attach handler to logger
logger.setLevel(logging.WARNING)    # sets debug level for logger


import sys

import urllib   # Import module to check if website is up
import time     # Import time for time.sleep()

import numpy

#----------------------------------------------------------------------------------------
## First CHECK IF SERVER IS UP 
def CheckServer():      # CodeFolding Shortcut: option+command+[  
    ServerUp=0; # ServerUp=1 if server is up, 0 if server is down
    try:
        a=urllib.urlopen("http://www.gmail.com").getcode() # a=200 when website can be opened
        if a==200:
            ServerUp=1;    # Server is up (define a=1)
    except IOError: # If get an IOError (Can't connect to google)
        #print('IOError, internet up?')
        ServerUp=0; 
        # No msg when python can connect to google
    return ServerUp # Returns Var. ServerUp


# SERIAL READ IN FUNCTION
def ReadSerial():       # Inputs: None. Outputs: x,T,x2,T2,WixelAlive
    try:
        # Setup SerialData Read
        port='/dev/ttyACM0'
        ser = serial.Serial(port,9600,timeout=50)    # When timeout occurs, nothing is read in and code goes to EXCEPTION: WixelAlive=0 (Dead)
        value = ser.readline()

        x=float(value.split()[2])   # [2] for ch1, [3] for ch2, [4] for ch3 monitor battery
        x2=float(value.split()[3])
        T=(x-500)/10 # T is degrees in celsius
        T2=(x2-500)/10

        # FAKING DATA (Uncomment while Wixel connected and ACTUALLY reading data in)
        # x=700; # x is raw voltage from TMP36 (max Vwixel=3400mV)
        # T=140; # T is temperature
        # x2=732;
        # T2=45;

        WixelAlive=1;
        return (WixelAlive, x, T, x2, T2)
    except Exception:
        WixelAlive=0; # WIXEL PROBLEM (CAN'T READ IN)
        
        x=numpy.NaN;    # numpy.NaN are float numbers
        T=numpy.NaN;
        x2=numpy.NaN;
        T2=numpy.NaN;
        return (WixelAlive, x, T, x2, T2)


# LOGGER FUNCTION
def Logger(Tnorm,WixelAlive,ServerUp,x,T,x2,T2):    # Outputs: Logs Tnorm,WixelAlive,ServerUp   
    # LOGGER (31 bytes per line * 1440 lines per day * 30 days = 1.3 Mb per month)
    # logger = logging.getLogger('TemperatureLogFile') # create logger instance
    # hdlr = logging.handlers.TimedRotatingFileHandler("TemperatureLogFile",when="midnight") # create filehandler (handler sends log to specified destination)
    # formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')    # create formatter (Layout of log)

    # hdlr.setFormatter(formatter)    # attach formatter to handler
    # logger.addHandler(hdlr)         # attach handler to logger
    # logger.setLevel(logging.WARNING)    # sets debug level for logger

    if Tnorm==0:    # T is NOT normal => TEMP IS HIGH !!
        if WixelAlive==0:   # Wixel is DEAD !!
            if ServerUp==0: # If Server is down..
                logger.error('T1=' + repr(T) + ' T2='+ repr(T2) 
                    + ' Temp DUNNO, Wixel DEAD, Server DOWN (1) ') #%(Tnorm,WixelAlive,ServerUp) )      # This is the message 
                print 'Tnorm, WixelAlive, ServerUp = %i%i%i' %(Tnorm,WixelAlive,ServerUp) # Prints action
                print 'Temp DUNNO, Wixel DEAD, Server DOWN -> Log & Cant Email (1)'   
            
            elif ServerUp==1:   # Server is UP
                logger.error('T1=' + repr(T) + ' T2='+ T2
                    + ' Temp DUNNO, Wixel DEAD, Server Up -> Email (2)') #%(Tnorm,WixelAlive,ServerUp) )      # This is the message 
                print 'Tnorm, WixelAlive, ServerUp = %i%i%i' %(Tnorm,WixelAlive,ServerUp) # Prints action
                print 'Temp DUNNO, Wixel DEAD, Server is Up -> Log & Email (2)'
        
        
        elif WixelAlive==1:   # Wixel is Alive
            if ServerUp==0:     # Server is DOWN
                logger.error('T1=' + repr(T) + ' T2='+ repr(T2)
                    + ' Temp HIGH, Wixel Alive, Server DOWN (3) ') #%(Tnorm,WixelAlive,ServerUp) )      # This is the message 
                print 'Tnorm, WixelAlive, ServerUp = %i%i%i' %(Tnorm,WixelAlive,ServerUp) # Prints action
                print 'Temp HIGH, Wixel Alive, Server DOWN -> Log & Cant Email. (3)'
                    
            elif ServerUp==1: # Server is UP
                logger.error('T1=' + repr(T) + ' T2='+ repr(T2)
                    + ' Temp HIGH, Wixel Alive, Server Up -> Email (4)') #%(Tnorm,WixelAlive,ServerUp) )      # This is the message 
                print 'Tnorm, WixelAlive, ServerUp = %i%i%i' %(Tnorm,WixelAlive,ServerUp) # Prints action
                print 'Temp HIGH, Wixel Alive, Server is Up -> Log & Email. (4)'

    
    elif Tnorm==1:  # Temp is normal
        if WixelAlive==0: # Wixel is Dead !!
            if ServerUp==0: # Server is DOWN !!
                logger.error('T1=' + repr(T) + ' T2='+ repr(T2)
                    + ' Temp DUNNO, Wixel DEAD, Server DOWN (5)') #%(Tnorm,WixelAlive,ServerUp) )      # This is the message 
                print 'Tnorm, WixelAlive, ServerUp = %i%i%i' %(Tnorm,WixelAlive,ServerUp) # Prints action
                print 'Temp DUNNO, Wixel DEAD, Server DOWN -> Log & Cant Email (5)'   
            
            elif ServerUp==1:   # Server is UP
                logger.error('T1=' + repr(T) + ' T2='+ repr(T2)
                    + ' Temp DUNNO, Wixel DEAD, Server Up -> Email (6)') #%(Tnorm,WixelAlive,ServerUp) )      # This is the message 
                print 'Tnorm, WixelAlive, ServerUp = %i%i%i' %(Tnorm,WixelAlive,ServerUp) # Prints action
                print 'Temp DUNNO, Wixel DEAD, Server is Up -> Log & Email (6)'

        elif WixelAlive==1: # Wixel Alive
            if ServerUp==0: # Server is DOWN !!
                logger.error('T1=' + repr(T) + ' T2='+ repr(T2)
                    + ' Temp norm, Wixel Alive, Server DOWN (7)') #%(Tnorm,WixelAlive,ServerUp) )      # This is the message 
                print 'Tnorm, WixelAlive, ServerUp = %i%i%i' %(Tnorm,WixelAlive,ServerUp) # Prints action
                print 'Temp norm, Wixel Alive, Server DOWN -> Log & Cant Email (7)'   
            
            elif ServerUp==1:   # Server is UP
                logger.error('T1=' + repr(T) + ' T2='+ repr(T2)
                    + ' Temp norm, Wixel Alive, Server Up (8)') #%(Tnorm,WixelAlive,ServerUp) )      # This is the message 
                print 'Tnorm, WixelAlive, ServerUp = %i%i%i' %(Tnorm,WixelAlive,ServerUp) # Prints action
                print 'Temp norm, Wixel Alive, Server is Up -> Log (8)'
		EmailCount = 0
		EmailThreshold = 1 # back to normal email threshold


#----------------------------------------------------------------------------------------
## SEND EMAIL FUNCTION
def Email(ServerUp,EmailCount, EmailThreshold):
    if (ServerUp==1 and Tnorm==0) or (ServerUp==1 and WixelAlive==0):     # 'and' so that don't send email on case (8) (Normal Operation)
        print "Email: Request"
        EmailCount=EmailCount+1;

        if EmailCount > EmailThreshold:     # Send an email if we are above the Threshold
	    EmailThreshold = EmailThreshold*2 # this should lead to message every 1min, 2min, 4min, 8min, 16min etc after disaster...
            logger.error('Email Actually sent!')
            # Open .txt file for sending latest temperature values  
            fp = open(LogFileName, 'rb')
            # Create a text/plain message
            msg = MIMEText(fp.read())
            fp.close()

            me = 'ubcmri7t@gmail.com'
            you = ['stefan@phas.ubc.ca'] #,'stefan@phas.ubc.ca','clayton.wong1x@gmail.com']   # Can send to 2 recipients
            #msg = MIMEMultipart('alternative')
            msg['Subject'] = 'Alert: 7T Freezer Temperature HOTT'  # "Hello %s, my name is %s" % ('john', 'mike') # Hello john, my name is mike".
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
    elif ServerUp==0:   #Server is down => Don't send email
        print ('Email: No email because server is DOWN')

    return EmailCount, EmailThreshold




##----------------------------------------------------------------------------------------
## MAIN LOOP
EmailCount=0;
EmailThreshold=1;
print EmailCount

while True: # Infinite loop (ctr-c to break)

    
    # Check Server
    ServerUp=CheckServer();         # Calls Function 'CheckServer'
    #print "ServerUp is : ", ServerUp


    #while one==1: # Infinite loop (ctr-c to break)

    # Read Serial in
    [WixelAlive, x, T, x2, T2]=ReadSerial();    # Calls Function ReadSerial
    #print "WixelAlive is : ", WixelAlive
    print ("x, T, x2, T2 are : %s, %s, %s, %s"  % (x, T, x2, T2))


    # Logger
    if (T > T1HIGH) or (T2 > T2HIGH):   # High Temperature
        Tnorm=0;
    else:         # Normal Temperature
        Tnorm=1;
        
    Logger(Tnorm,WixelAlive,ServerUp,x,T,x2,T2)    # Calls Function 'Logger'
    
    # Call 'Email' Function  # Actually send an email if we are above EmailThreshold
    (EmailCount, EmailThreshold) = Email(ServerUp, EmailCount, EmailThreshold)     
    print EmailCount

    del x, T, x2, T2
    time.sleep(60)


# To do: 
# - Alert when Temp too COLD?
