# -*- coding: utf-8 -*-
"""
Testing the Logger
"""
import unittest
from datetime import datetime

class Test(unittest.TestCase):
    '''
    Unit tests for SARlogger.
    '''
    
    def test_SARlogger(self):

        #setup this module should be using a logger
        from SARlabpy.io import BRUKERIO
        self.assertIsNotNone(BRUKERIO.logger)
        
        #or we should be able to query it
        import logging
        logger = logging.getLogger('root')
        self.assertIsNotNone(logger)
        
        #let us do some logging
        
        logger.debug('debug message on {0}'.
                      format(datetime.now().strftime('%c')))
        logger.info('confirming that all is working well')
        logger.warning('there are no warnings yet')
        logger.error('nor are there any errors')
        logger.critical('not to mention something critical that would make us wheep')

        logger.setLevel(logging.CRITICAL)
        logger.debug('this message goes unheeded')    
        logger.critical('whereas this one is loud enough')


if __name__ == "__main__":
    unittest.main()