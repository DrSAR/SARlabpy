version = '4.0'
release = version + 'beta'
from helpers import git_repo_state
repo_state = git_repo_state()

import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

from inputoutput.BRUKERIO import (readJCAMP, readfid, readfidspectro, read2dseq,
                        dict2string, fftbruker, readRFshape)

from inputoutput.BRUKER_classes import (Scan, Study, StudyCollection,
                               Patient, Experiment,
                               dataroot)

from inputoutput.AData_classes import (adataroot, AData)

from inputoutput import (rois, write_csv, bulk_analyze, mriBoards)

from inputoutput.visu_pars_2_Nifti1Header import (visu_pars_2_Nifti1Header,
                                      visu_pars_2_matrix_size)
                                      
from inputoutput.SARlogger import initiate_logging

from helpers import natural_sort, generate_summary, smooth_SG

from fmoosvi import roipoly
