# -*- coding: utf-8 -*-
"""
Testing the Logger
"""
import unittest

#make the SARlogger features available
from sarpy import SARlogger

#test logging for BRUKERIO
import test.logtest_module
SARlogger.initiate_logging(test.logtest_module)


class Test(unittest.TestCase):
    '''
    Unit tests for SARlogger.
    '''
    
    def test_SARlogger(self):
        '''
        This should setup the logger for a test module.
        You will notive how the levels for the console and the file logging
        can be independent.
        '''
        test.logtest_module.foo()
        #switch to INFO only
        print('\nswitch to INFO level')
        SARlogger.change_logger_level(
                test.logtest_module,
                level=SARlogger.INFO)
        test.logtest_module.foo()
        #switch to WARNING only
        print('\nswitch to WARNING level')
        SARlogger.change_logger_level(
                test.logtest_module,
                level=SARlogger.WARNING)
        test.logtest_module.foo()
                
        # Now, let's see some file action
        print('\nadd file logger and switch to DEBUG level (for logger)')
        SARlogger.change_logger_level(
                test.logtest_module,
                level=SARlogger.DEBUG)
    
        SARlogger.add_file_handler(test.logtest_module)
        test.logtest_module.foo(msg='file logging?')

if __name__ == "__main__":
    unittest.main()