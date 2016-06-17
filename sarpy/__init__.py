version = '5.0'
release = version + 'beta'
import os
import json

from .helpers import natural_sort, smooth_SG, git_repo_state
repo_state = git_repo_state()

import logging
import logging.config
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())
# the logging configuration is kept in a json file at the
# root directory of the package: sarpy/logging.json
conf_name = os.path.join(__path__[0], 'logging.json')
with open(conf_name, 'rt') as f:
    logging.config.dictConfig(json.load(f))
logger.info('restarting sarpy application')

#from . import BRUKERIO
from .BRUKERIO.BRUKER_classes import (Scan, Study, StudyCollection,
                            Patient, Experiment, dataroot)

from .BRUKERIO.AData_classes import (adataroot, AData)

print(('Imported sarpy v{0}\n'+
      'git branch: {1} ({2}) -- {3}').format(
      release,
      repo_state['branch'],
      'dirty' if repo_state['dirty'] else 'clean',
      repo_state['sha1'][0:10]))
#from inputoutput import (rois, write_csv, bulk_analyze, mriBoards)

#from inputoutput.visu_pars_2_Nifti1Header import (visu_pars_2_Nifti1Header,
#                                      visu_pars_2_matrix_size)
                                      
#from fmoosvi import roipoly
