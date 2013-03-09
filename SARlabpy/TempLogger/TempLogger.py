# Temperature Logger (TempLogger2 uses logger)
# Clayton Wong
# Feb 10th, 2013

import serial # import serial for serial input
#import datetime # import datetime for getting localtime
import smtplib # import smtplib simple mail transfer protocol for EMAIL
# Import the email modules we'll need 
from email.mime.text import MIMEText
#from email.mime.multipart import MIMEMultipart

import logging, logging.handlers    # Import logger modules
import sys

import urllib   # Import module to check if website is up

#----------------------------------------------------------------------------------------
## First CHECK IF SERVER IS UP 
def CheckServer():      # CodeFolding Shortcut: option+command+[  
    ServerUp=0; # ServerUp=1 if server is up, 0 if server is down
    try:
        a=urllib.urlopen("http://www.google.com").getcode() # a=200 when website can be opened
        if a==200:
            ServerUp=1;    # Server is up (define a=1)
    except IOError: # If get an IOError (Can't connect to google)
        #print('IOError, internet up?')
        ServerUp=0; 
        # No msg when python can connect to google
    return ServerUp # Returns Var. ServerUp
     
    # Setup SerialData Read     # UNCOMMENT WHEN READING SERIAL IN
    #port='/dev/tty.usbmodemfa131'
    #ser = serial.Serial("/dev/tty.usbmodemfa131",9600)
    #print ser


def ReadSerial():
    try:
        #value = ser.readline();
        #print value
        #print value.split()[2] # This reads channel 1 (the [2] reads the 2nd col sep by white space)

        #x=float(value.split()[2])   # [2] for ch1, [3] for ch2, [4] for ch3 monitor battery
        #x2=float(value.split()[3])
        #T=(x-500)/10 # T is degrees in celsius
        #T2=(x2-500)/10
        #print x , T , x2 , T2

        x=700; # x is raw voltage from TMP36 (max Vwixel=3400mV)
        T=120; # T is temperature
        x2=732;
        T2=45;
        WixelAlive=1;
        return (WixelAlive, x, T, x2, T2)
    except IOError:
        WixelAlive=0; # WIXEL PROBLEM (CAN'T READ IN)
        x="error";
        T="error";
        x2="error";
        T2="error";
        return (WixelAlive, x, T, x2, T2)


def Logger(Thigh,WixelAlive,ServerUp,x,T,x2,T2):        
    # LOGGER (31 bytes per line * 1440 lines per day * 30 days = 1.3 Mb per month)
    logger = logging.getLogger('TemperatureLogFile') # create logger instance
    hdlr = logging.handlers.TimedRotatingFileHandler("TemperatureLogFile",when="midnight") # create filehandler (handler sends log to specified destination)
    formatter = logging.Formatter('%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')    # create formatter (Layout of log)

    hdlr.setFormatter(formatter)    # attach formatter to handler
    logger.addHandler(hdlr)         # attach handler to logger
    logger.setLevel(logging.WARNING)    # sets debug level for logger

    if Thigh==1:    # TEMP IS HIGH
        if WixelAlive==1:   # Wixel is ALIVE
            if ServerUp==1: # If Server is up..
                logger.error('ch1(V/3.3V,Temp)=' + repr(x) + ',' + repr(T) + ' ch2(V/3.3V,Temp)=' + repr(x2) + ',' + repr(T2) + ' Emailed')      # This is the message 
                print ("Thigh, WixelAlive, ServerUp = %i%i%i => Log,Email" %(Thigh,WixelAlive,ServerUp))
            elif ServerUp==0:   # Server is DOWN!!!!
                # Log "server down" in logfile, but can't send email (since server is down)
                logger.error('ch1(V/3.3V,Temp)=' + repr(x) + ',' + repr(T) + ' ch2(V/3.3V,Temp)=' + repr(x2) + ',' + repr(T2) + ' SERVER DOWN')      # This is the message 
                print "WixelAlive, ServerDown => Log" 
        
        elif WixelAlive==0:   # Wixel is DEADD
            if ServerUp==1:
                logger.error('ch1(V/3.3V,Temp)=' + repr(x) + ',' + repr(T) + ' ch2(V/3.3V,Temp)=' + repr(x2) + ',' + repr(T2) + ' WIXEL DEAD')      # This is the message 
                print "WixelDead, ServerUp => Log"        
            elif ServerUp==0:
                logger.error('ch1(V/3.3V,Temp)=' + repr(x) + ',' + repr(T) + ' ch2(V/3.3V,Temp)=' + repr(x2) + ',' + repr(T2) + ' WIXEL DEAD AND SERVER DOWN')
                print "WixelDead, ServerDown => Log"

    elif Thigh==0:
        print 'Thigh == 0'


    # if T < 0:
    #     logger.error('THIS IS SECOND LINE OF MSG (TEMP IS TOO LOW)')
    #     #logger.error('3rd line of log')


#----------------------------------------------------------------------------------------
## SEND EMAIL 
def Email(ServerUp):
    if ServerUp==1:
        print "Email Sent"
        # Open .txt file for sending latest temperature values  
        fp = open("TemperatureLogFile", 'rb')
        # Create a text/plain message
        msg = MIMEText(fp.read())
        fp.close()

        me = 'ubcmri7t@gmail.com'
        you = ['clw9@sfu.ca'] #,'clayton.wong1x@gmail.com']   # Can send to 2 recipients
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
    else:   #Server is down => Don't send email
        print ('No email because server is DOWN')




# Other STuff to implement:
# 1. DONE Cannot connect to gmail server (ev Hour?) LOG THIS
# 2. Stop Reading in voltage (ev 30 mins?) SEND EMAIL
# 3. DONE Temperature too cold   



##----------------------------------------------------------------------------------------
## MAIN LOOP

# Check Server
ServerUp=CheckServer();         # Calls Function CheckServer
print "ServerUp is : ", ServerUp


#while one==1: # Infinite loop (ctr-c to break)

# Read Serial in
[WixelAlive, x, T, x2, T2]=ReadSerial();    # Calls Function ReadSerial
if WixelAlive==1:   # WixelAlive
    print "WixelAlive is : ", WixelAlive
    print ("x, T, x2, T2 are : %i, %i, %i, %i"  % (x, T, x2, T2))
    #print ('%s blahblahblah %s' % ("string1", "string2"))
elif WixelAlive==0:    # Wixel is DEAD
    print "WixelAlive is : ", WixelAlive
    print ("x, T, x2, T2 are : %s, %s, %s, %s"  % (x, T, x2, T2))

# Logger
if T > 100:   # High Temperature
    Thigh=1;
    Logger(Thigh,WixelAlive,ServerUp,x,T,x2,T2)    # Calls Function Logger
    Email(ServerUp)     # Calls Function Email
    #print ("TempHigh => Log, Email (if can)")

else:   # Normal Temperature
    Thigh==0;
    Logger(Thigh,WixelAlive,ServerUp,x,T,x2,T2)    # Calls Function Logger
    #print ("TempNorm => Log")







