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

import logging
from logging import (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    
def initiate_logging(module, formatter=None, 
                     handler=None, handler_level=None,
                     logger_level=None):
    '''
    Initiate logging for a module (in the SARlabpy package). This assumes 
    a module has been imported.
    
    :param module module:
        Module object that points to module in the SARlabpy package 
        framework (collection of libraries? modules?).
    :param int handler_level: 
        log level, default logging.DEBUG(=10)
        this will be ignored if there is no handler provided and one or 
        more handlers pre-exist. in other words, the level can only be set 
        for the default (console) handler. Use :py:func: change_logger_level:
        to change the overall level
        
    The module in question will have set up minimal logging. A NullHandler 
    has been set so the default logging (methods in the libraries should
    emit appropriate log message at opportune times) is swallowed. If this
    handlers is overwritten, logging can commence.
    '''
    if formatter is None:
        formatter = logging.Formatter(
                fmt='%(asctime)s %(levelname)s: %(module)s  - %(message)s',
                datefmt='%H:%M')

    handler_level = None or WARNING
    logger_level = None or DEBUG
            
    # remove all NullHandlers:
    handlers_to_be_removed = []
    for hdlr in module.logger.handlers:
        if isinstance(hdlr, logging.NullHandler):
            handlers_to_be_removed.append(hdlr)
    for hdlr in handlers_to_be_removed[:]:
        module.logger.removeHandler(hdlr)

    if (handler is None and len(module.logger.handlers) == 0):
        # no handler provided and the previously defined
        # handlers are NullHandler -> adding default console handler
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        handler.setLevel(handler_level)
    if handler is not None:
        handler.setFormatter(formatter)
        handler.setLevel(handler_level)

    try:
        module.logger.addHandler(handler)
        module.logger.setLevel(logger_level)
    except AttributeError, exc:
        # TODO reformat output
        print('Initiation of logging for module {0} failed'.
                format(module))
        print('Presumably logging has not been prepared.')
        raise AttributeError(exc)
    
        
           
def add_file_handler(module, formatter=None, handler=None, level=None):
    '''
    add handler for logging to file. typically with lower level and timestamps
    '''
    if formatter is None:
        formatter = logging.Formatter(fmt=
                        '%(asctime)s %(levelname)s: %(name)s  - %(message)s',
                        datefmt='%y-%m-%d %H:%M')
    
    if handler is None:
        handler = logging.FileHandler('/tmp/SARlabpy.log')
    if level is None:
        level = DEBUG
    handler.setFormatter(formatter)
    handler.setLevel(level)
    module.logger.addHandler(handler)    

def change_logger_level(module, level=None):
    '''
    This changes the overall logger level. Above in :py:func: initiate_logging,
    the level for the individual handler was being set.
    '''
    module.logger.setLevel(level=level)
