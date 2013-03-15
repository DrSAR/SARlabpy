#! /bin/sh
""":"
exec python $0 ${1+"$@"}
"""
# Temperature Logger (TempLogger2 uses logger)
# Clayton Wong
# Feb 10th, 2013

import serial # import serial for serial input
#import datetime # import datetime for getting localtime

import logging, logging.handlers    # Import logger modules
import sys

#----------------------------------------------------------------------------------------
# Setup SerialData Read     # UNCOMMENT WHEN READING SERIAL IN
port='/dev/tty.usbmodemfa131'
try:
    ser = serial.Serial("/dev/tty.usbmodemfa131",9600, timeout=2)
except:
    print("problem seeing the wixel")
else:
    print ser

    # MAIN LOOP
    one=1
    while True:
        try:
            value = ser.readline();
            if value:
                print value
            else:
                print ('other wixel dead')
        except OSError:
            print("problem reading from serial")
        sys.stdout.flush()
        #print value.split()[2] # This reads channel 1 (the [2] reads the 2nd col sep by white space)

        #x=float(value.split()[2])   # [2] for ch1, [3] for ch2, [4] for ch3 monitor battery
        #x2=float(value.split()[3])
        #T=(x-500)/10 # T is degrees in celsius
        #T2=(x2-500)/10
        #print x , T , x2 , T2

