version = '5.0'
release = version + 'beta'
from .helpers import natural_sort, smooth_SG, git_repo_state
repo_state = git_repo_state()

import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

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
                                      
#from inputoutput.SARlogger import initiate_logging

#from fmoosvi import roipoly
