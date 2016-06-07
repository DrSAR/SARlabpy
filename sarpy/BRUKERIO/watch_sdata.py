# -*- coding: utf-8 -*-
"""
Watch directory and trigger processing

@author: stefan a reinsberg
Created on Thu Jul 25 16:36:34 2013
"""

import os
import pyinotify
import datetime
import time
import hashlib
import sarpy.io.mriBoards

class PTmp(pyinotify.ProcessEvent):
    def process_IN_CREATE(self, event):
        print "Create: %s" %  os.path.join(event.path, event.name)
    def process_IN_DELETE(self, event):
        print "Remove: %s" %  os.path.join(event.path, event.name)
    def process_IN_MODIFY(self, event):
        fname = os.path.join(event.path, event.name)
        print "Modify: {0} ({1})".format(fname,
                                         datetime.datetime.now().strftime('%X'))
        try:
            sarpy.io.mriBoards.generate(conf_file=fname)        
        except Exception as e:
            print e
       
wm = pyinotify.WatchManager()
# watched events
mask = pyinotify.ALL_EVENTS

pm = PTmp()

notifier = pyinotify.Notifier(wm, pm)

wm.add_watch(os.path.join(os.path.expanduser('~'),'sdata','mriBoards'), 
                                mask, rec=True, quiet=False)

while True:
    try:
        # process the queue of events as explained above
        notifier.process_events()
        if notifier.check_events():
            # read notified events and enqueue them
            notifier.read_events()
           
    except KeyboardInterrupt:
        # destroy the inotify's instance on this interrupt (stop monitoring)
        notifier.stop()
        print('Keyboard interrupt detected')
        break