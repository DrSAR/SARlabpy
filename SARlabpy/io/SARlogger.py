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
   logger=SARlogger.setup_custom_logger(level=SARlogger.WARNING)

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
from logging import (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    
def initiate_logging(module, formatter=None, handler=None, level=None):
    '''
    Initiate logging for a module (in the SARlabpy package). This assumes 
    a module has been imported.
    
    :param module module:
        Module object that points to module in the SARlabpy package 
        framework (collection of libraries? modules?).
    :param int level: 
        log level, default logging.DEBUG(=10)
        
    The module in question will have set up minimal logging. A NullHandler 
    has been set so the default logging (methods in the libraries should
    emit appropriate log message at opportune times) is swallowed. If this
    handlers is overwritten, logging can commence.
    '''
    if formatter is None:
        formatter = logging.Formatter(fmt=
            '%(asctime)s - %(levelname)s - %(module)s - %(message)s')
    if handler is None:
        handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    if level is None:
        level = DEBUG
    try:
        module.logger.setLevel(level)
        print module.logger.handlers
        if len(module.logger.handlers)<1:
            module.logger.addHandler(handler)
    except AttributeError, exc:
        print('Initiation of logging for module {0} failed'.
                format(module))
        print('Presumably logging has not been prepared.')
        raise AttributeError(exc)
        
def add_file_logging_handler(module, formatter=None, handler=None, level=None):
    '''
    add handler for logging to file. typically with lower level and timestamps
    '''
    if formatter is None:
        formatter = logging.Formatter(fmt=
                        '%(asctime)s %(name)-12s %(levelname)-8s '+
                        '%(module)s - %(message)s',
                        datefmt='%m-%d %H:%M')
    
    if handler is None:
        handler = logging.FileHandler('/tmp/SARlabpy.log')
    if level is None:
        level = DEBUG
    logging.getLogger(module.__name__).addHandler(handler)    

def change_logging_level(module, level=None):
    '''
    Since the :py:func: initiate_logging method is careful not toadd a handler if
    one already exists, changing log levels can de done directly there ...
    initiate_loggin(module, level=level)
    '''
    initiate_logging(module, level=level)