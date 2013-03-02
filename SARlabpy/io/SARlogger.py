# -*- coding: utf-8 -*-
"""
Custom Logger to be used across modules

To setup logging in your code. Your module needs to initiate logging using
the custom SARlogger.Good info to be found at the 
`Logging Tutorial 
<http://docs.python.org/2/howto/logging.html#logging-basic-tutorial>`_.
    
    
.. code-block:: py

   #make the SARlogger features available
   from SARlabpy import SARlogger
   #test logging for BRUKERIO
   import BRUKERIO
   SARlogger.initiate_logging(BRUKERIO)

At this point, a console handler that issues anything that's a WARNING or 
more severe (higher level) will be written to to the console. 
If you  also desire output to 
a file, another convenience function in SARlogger can be used
to attach another handler (with another log level and formatter): 

.. code-block:: py

   SARlogger.add_file_handler(BRUKERIO)
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
        log level for the handler, default logging.DEBUG(=10)
    :param int logger_level: 
        log level for the entire logger, default logging.DEBUG(=10)
        Use :py:meth:`change_logger_level` to change the overall level
    :param formatter:
        default provides time, level, module and message
    :param handler:
        default provides a StreamHandler (that prints to the console)
        
    The module in question will have set up minimal logging. A NullHandler 
    has been set so the default logging (methods in the libraries should
    emit appropriate log message at opportune times) is swallowed. If this
    handlers is overwritten, logging can commence.
    '''
    if formatter is None:
        formatter = logging.Formatter(
                fmt='%(asctime)s %(levelname)s: %(module)s  - %(message)s',
                datefmt='%H:%M')

    handler_level = handler_level or WARNING
    logger_level = handler_level or DEBUG
            
    # remove all NullHandlers:
    handlers_to_be_removed = []
    try:
        for hdlr in module.logger.handlers:
            if isinstance(hdlr, logging.NullHandler):
                handlers_to_be_removed.append(hdlr)
    except AttributeError, exc:
        # TODO reformat output
        print('Initiation of logging for module {0} failed'.
                format(module))
        print('Presumably logging has not been prepared.')
        raise AttributeError(exc)    
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

    module.logger.addHandler(handler)
    module.logger.setLevel(logger_level)
    
        
           
def add_file_handler(module, formatter=None, handler=None, level=None):
    '''
    Add handler for logging to file. typically with lower level and timestamps
    
    :param module module: module whose logger to add the handler to
    :param int level: 
        log level for the entire logger, default logging.DEBUG(=10)
        This is equivalent to handler_level option in
        :func:`initiate_logging`
    :param formatter:
        default provides date, time, level, module and message
    :param handler:
        default provides a FileHandler (that prints to /tmp/SARlabpy.log)
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

    :param module module: module whose logger to change the logger level for
    :param int level: 
        log level for the entire logger, default logging.DEBUG(=10)
        This is equivalent to logger_level option in
        :func:`initiate_logging`
    '''
    module.logger.setLevel(level=level)
