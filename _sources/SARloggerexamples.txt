How to use Module Logging
=========================

Modules in this package (should) provide a thin layer to the python
logging module. The idea is to emit messages inside the library and
let the user of said library decide on whether they want to see it. 

Module setup
-------------
Any module interested in logging should include the following header:

.. code-block:: py

   import logging
   logger=logging.getLogger(__name__)
   logger.addHandler(logging.NullHandler())

Inside the module being watch, statements of the kind:

.. code-block:: py

   logger.debug('Variable x={0}'.format(x))
   logger.warning('File not found, using random data instead. Good luck')

will lead to debug messages. With the default setting the first will be swallowed
(since it is not severe enough, i.e. low enough log level), while the second should
be visible.

Configuration by the Module User
--------------------------------
The importer of the module (e.g. BRUKERIO) needs to issue the following 
statement in addition to the actual import statement:

.. code-block:: py

   #make the SARlogger features available
   from SARlabpy import SARlogger
   #test logging for BRUKERIO
   import BRUKERIO
   SARlogger.initiate_logging(BRUKERIO)

At this point, a console handler that issues anything that's a WARNING or 
more sever will be written to to the console. 
Without the call to py:func::SARlogger.initiate_logging , 
logging messages will be ignored alltogether. 

If you desire output to a file, another convenience function in SARlogger can be used
to attach another handler (with another log level and formatter): 

.. code-block:: py

   SARlogger.add_file_handler(BRUKERIO)

By default, this will lead to logging at the DEBUG level to file /tmp/SARlabpy.log.

Log levels
----------
Log levels can be changed both at the logger level and the handler level. Note that no
log messages are handed on to the handlers if they are above the logger log level. 
Each individual handler in turn decides whether to publish a received log message or
not. This way it is possible to reduce console output to the more severe messages
and log every little fart in a log file.


More sophisticated log control
-------------------------------
Since the module will expose access to the logger through module attribute `logger`,
logging is as configurable as the python logging module itself. You can attach
further handlers, change the levels, come up with new formatters and filters
as you desire....
