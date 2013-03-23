'''
Central __init__.py file that pulls various typically used functions from the
various subpackes into the main sarpy namespace. This way, e.g., instead of using the
cumbersome :py:class:`sarpy.io.BRUKER_classes.Scan()` users can use the rather 
snappy :py:class:`sarpy.Scan()`.
'''
version = '0.1'
release = version + 'alpha'

import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from io.BRUKERIO import (readJCAMP, readfid, readfidspectro, read2dseq,
                        dict2string, fftbruker, readRFshape)

from io.BRUKER_classes import (natural_sort, 
                               Scan, Study, StudyCollection,
                               Patient, Experiment,
                               dataroot)
from io.AData_classes import (adataroot, AData)
                        
from io.SARlogger import initiate_logging
