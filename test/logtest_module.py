# -*- coding: utf-8 -*-
"""
Test module for SARlogger features

The logging is set up so that if the library user does nothing,
all will be silent. Details in :py:mod: SARlogger
"""
import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def foo(msg='message from foo'):
    logger.debug(msg)
    logger.info(msg)    
    logger.warning(msg)
    
