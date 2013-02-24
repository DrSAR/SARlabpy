# -*- coding: utf-8 -*-
"""
Custom Logger to be used across modules

To setup logging in your code. Your module needs to initiate logging using
the custom SARlogger.Good info to be found at the 
`Logging Tutorial 
<http://docs.python.org/2/howto/logging.html#logging-basic-tutorial>`_.
    
    
.. code-block:: py

   # setup logging
   import SARlogger
   logger=SARlogger.setup_custom_logger('root')

Currently BRUKERIO instantiates such a logger. You can get access to it through 
the module variable: BRUKERIO.logger or you can access it by finding where the
'root' logger is:

.. code-block:: py

   import logging
   logger = logging.getLogger('root')
   logger.setLevel(logging.DEBUG)

In either case you can then emit messages in your code 
(whereever logger is accessible):

.. code-block:: py
   
   logger.debug('this is a debug message')
"""

print('SARlogger imported: %s' % __name__)
import logging

def setup_custom_logger(name, level=None):
    level = level or logging.DEBUG
    formatter = logging.Formatter(fmt=
            '%(asctime)s - %(levelname)s - %(module)s - %(message)s')

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    while len(logger.handlers)>1:
        logger.removeHandler(logger.handlers[0])
    return logger