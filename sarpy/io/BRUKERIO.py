#!/usr/bin/env python

# Copyright (C) 2012-2013 Stefan A Reinsberg and SARlab members
# full license details see LICENSE.txt
"""Collection of BRUKER input routines

Handy functions to read BRUKER data and header files.

The logging is set up so that if the library user does nothing,
all will be silent. Details in :py:mod: SARlogger
"""
from __future__ import division
import logging
logger=logging.getLogger('sarpy.io.BRUKERIO')

import numpy
import os.path
import re
from types import StringType, FileType, UnicodeType
from itertools import tee, izip

def pairwise(iterable):
    """
    This is a solution to the problem of looking ahead in a for loop
    mentioned in on `stackoverflow: <http://stackoverflow.com/questions/4197805/python-for-loop-look-ahead>`_

    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    """
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

def convert_int_float_string(param):
    try:
        y = int(param)
    except ValueError:
        try:
            y = float(param)
        except ValueError:
            y = param.strip()
    return y

def readJCAMP(filename):
    """
    Parse text file in JCAMP format

    :param string filename: filename of fid file
    :return:
        Dictionary of labelled data records (LDR) with LDR-names as keys
        and their content as values.
    :rtype: dict
    :raises: IOERROR if opening file fails or passes on any other error

    The *JCAMP format* is a self-documenting, ASCII text file format
    that is maintained by IUPAC (`see a report
    here <http://iupac.org/publications/pac/78/3/0613/>`_).
    It consists of labelled data records (LDR) that start with a ##
    and end when the next record begins.
    They can span several lines. The data label is
    enclosed between '##' and '='. If it starts with a $, we are
    dealing with a private LDR.

    The issue of reading these is complicated due to the various
    types of data (integers, floats, strings, arrays and nested
    structures) that can be present.

    ~/bdata/readfidTest.ix1/subject:
    ::

        ##$SUBJECT_name_string=( 64 )
        <Moosvi, readfidTest>
        ##$SUBJECT_name=(<Moosvi>, <readfidTest>)

        >>> import os
        >>> a=readJCAMP(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/subject'))
        >>> a['SUBJECT_name']
        ['<Moosvi>', '<readfidTest>']
        >>> a['SUBJECT_name_string']
        'Moosvi, readfidTest'

    ~/bdata/readfidTest.ix1/1/acqp:
    ::

        ##$ACQ_user_filter=No
        ##$ACQ_dim_desc=( 2 )
        Spatial Spatial
        ##$NR=1
        ##$D=( 64 )
        0.00502257333333333 0 0.00066296 0.000114 0.000114 0 0.001886 0 2.5e-05 0
        0.000886 0.000886 0.000368 0 0 0 0 0 0 0 1e-05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        ##$TPQQ=( 16 )
        (<hermite.exc>, 17.9862515696541, 0) (<>, 30, 0) (<>, 30, 0) (<>, 30, 0) (<>,
        30, 0) (<>, 30, 0) (<>, 30, 0) (<>, 30, 0) (<>, 30, 0) (<>, 30, 0) (<>, 30, 0)
         (<>, 30, 0) (<>, 30, 0) (<>, 30, 0) (<>, 30, 0) (<>, 30, 0)
        ##$ACQ_grad_matrix=( 15, 3, 3 )
        1 0 -0 0 1 0 0 -0 1 1 0 -0 0 1 0 0 -0 1 1 0 -0 0 1 0 0 -0 1 -0 0 1 -0 1 -0 1
        0 0 -0 0 1 -0 1 -0 1 0 0 0 0 1 1 -0 0 0 1 -0 0 0 1 1 -0 0 0 1 -0 0 0 1 1 -0 0
        0 1 -0 1 0 -0 0 1 0 0 -0 1 1 0 -0 0 1 0 0 -0 1 -0 0 1 -0 1 -0 1 0 0 -0 0 1 -0
        1 -0 1 0 0 -0 0 1 -0 1 -0 1 0 0 0 0 1 1 -0 0 0 1 -0 0 0 1 1 -0 0 0 1 -0

        >>> import os
        >>> a=readJCAMP(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/1/acqp'))
        >>> a['ACQ_user_filter']
        'No'
        >>> a['ACQ_dim_desc']
        ['Spatial', 'Spatial']
        >>> a['NR']
        1
        >>> a['D'][0:4]
        [0.00502257333333333, 0, 0.00066296, 0.000114]
        >>> a['TPQQ'][0:4]
        [['<hermite.exc>', 17.9862515696541, 0], ['<>', 30, 0], ['<>', 30, 0], ['<>', 30, 0]]
        >>> a['TPQQ'][0]
        ['<hermite.exc>', 17.9862515696541, 0]
        >>> a['TPQQ'][0][1]
        17.9862515696541
        >>> a['ACQ_grad_matrix'][0]
        array([[1, 0, 0],
               [0, 1, 0],
               [0, 0, 1]])


    """

    logger.info("opening {0}".format(filename))
    JCAMPfile = open(filename, "r")
    JCAMPdata = JCAMPfile.read().splitlines() # read and lose the "\n"
    JCAMPfile.close()

    # let's loop through the file, remove comments and put each
    # LDR on its own line
    LDRlist = [] # start with empty list
#    LDRmetadict = {} # meta info about array dimensions for each entry
    JCAMPdata_iterator = pairwise(JCAMPdata)
    for line, next_line in JCAMPdata_iterator:
        # match location comment (of the form "$$ /*")
        if re.match("\\$\\$ \\/", line):
            line = [re.sub("\\$\\$ ", "##$FILE_LOCATION=", line)]
        # match date comment (assumes there is nothing else in the comments)
        elif re.match("\\$\\$ [^@]", line):
            line = re.sub("\\$\\$ ", "##$DATE=", line)
            uname = "##$USERNAME="+(line.split(' '))[-1]
            line = re.sub("[^ ]+$", "", line) # remove username
            line = [line, uname]
        # match all other comments
        elif re.match("\\$\\$ @", line):
            line = []
        # match lists that are normal LDRs
        elif re.match("##", line):
            # we also need to append the next line(s) if they are part of
            # the LDR
            while not re.match(r'[#\$]{2}', next_line):
                line=line+' '+ next_line
                newline, next_line = JCAMPdata_iterator.next()
            # we should keep information about array dimensions in a
            # separate meta storage list
            line = [line]
        # this must be a line that belongs to the preceeding LDR
        # attach this to the line that was previously appended to
        # the LDRlist
        else:
            raise IOError('encountered line that cannot be here:\n{0}'.
                    format(line,next_line))
#            LDRlist[-1] = LDRlist[-1] + " " + line
#            line = []
        #add this to the list of LDRs
        LDRlist.extend(line)

    # strip every labelled data record of its preceding ## or ##$.
    LDRlist = [re.sub('##[\\$]*', '', LDR) for LDR in LDRlist]

    # split every LDR at the " = " sign and turn it into a dictionary entry
    # make sure to specify maxsplit=1 to the str.split function in case
    # some parameter values contain = as a sign...
    LDRdict = dict([LDR.split("=", 1) for LDR in LDRlist])

    for k, v in LDRdict.iteritems():
        # is it an array or struct? (signified by the size indicator)
        if not re.match(r'\(', v):
            # we have an int/float/simple string
            val = convert_int_float_string(v)
            LDRdict[k] = val
            logger.debug('found {2}: {0}={1}'.format(k,LDRdict[k], type(val)))
        else:
            # this isunfortunately harder. Let's get the array dimensions:
            # is there something following the 'array' definition?
            struc_match = re.match(r'\(([^\)]*)\)$', v)
            if struc_match:
                # it starts and ends with a round bracket.
                # this is actually a structure (without dimension instructions)
                v = struc_match.group(1)
                shape = [1]
                LDRdict[k] = [convert_int_float_string(val)
                                for val in v.split(',')]
                logger.debug('found dim-less structure({2})): {0}={1}'.
                            format(k,v,shape))
            else:
                string_match = re.match(r'\(([^\)]*)\) *([^ ].*)', v)
                assert string_match is not None, 'Cannot parse parameter file'
                v = string_match.group(2)
                shape = [int(s) for s in string_match.group(1).split(',')]
                if re.match('<', v):
                    # this is a string, possibly a list of strings (all 
                    # demarcated  with <>)
                    ll = re.compile("\s*>\s+<\s*").split(v.strip('<>'))
                    LDRdict[k] = inner_value(ll)
                    logger.debug('found string): {0}={1}'.
                                format(k,LDRdict[k]))
                    continue
                elif re.match('\(', v):
                    # this is a (nested?) structure
                    list_level_one = [s.strip(' ()')
                                      for s in re.split('\) *\(',
                                      v)]
                    split_list = [list_element.
                            split(',') for list_element in list_level_one]

                    LDRdict[k] = [[convert_int_float_string(x) for x in listA]
                                 for listA in split_list]
                    logger.debug('found comlex structure): {0}={1}'.
                                format(k,LDRdict[k]))
                else:
                    # this is a proper array
                    split_string = [s for s in re.split(' ',v)]
                    if len(shape) > 1:
                        LDRdict[k] = numpy.array([convert_int_float_string(s) for
                                    s in split_string if s]).reshape(shape)
                    else:
                        LDRdict[k] = [convert_int_float_string(s) for
                                    s in split_string if s]
                    logger.debug('found array({1}): {0}={2}'.
                                format(k,shape,LDRdict[k]))
    return LDRdict

def inner_value(somelist):
    '''
    Return somelist[0] a one-element list or the whole list otherwise

    :param list somelist: list that might contain only one element
    :returns: element of somelist or somelist

    >>> inner_value([42])
    42
    >>> inner_value([42,43])
    [42, 43]
    >>> inner_value([[42]])
    42
    >>> inner_value([[42,43]])
    [42, 43]
    >>> inner_value('spam')
    'spam'
    >>> inner_value(['spam'])
    'spam'
    >>> inner_value(['s'])
    's'
    >>> inner_value(['spam','eggs','bacon'])
    ['spam', 'eggs', 'bacon']

    .. warning::
       This method does not work for dictionaries (KeyError)

    '''
    if isinstance(somelist,list):
        if len(somelist) == 1:
            return inner_value(somelist[0])
        else:
            return somelist
    else:
        return somelist


def readfid(fptr=None,
            acqp=None,
            method=None,
            untouched=False,
            squeezed=True,
            resetNR=False):
    """
    Returns BRUKER's fid file as a properly dimensioned & rearranged array.

    :param FileType,StringType fptr:
        filename of fid file or filehandle to open fid file
    :param dict acqp: dictionary of acqp parameters
        (default None: parameter file will be loaded)
    :param dict method:  dictionary of method parameters
        (default None: parameter file will be loaded)
    :param boolean untouched: do not rearrange/reshape data (default: False)
    :param boolean squeezed: squeeze dimension of length 1 (default: True)
    :param boolean resetNR: for incomplete scans, reset NR if needed (default: False)
    :return: Flat (untouched = True) or Rearranged and assembled array of the acquire k-space
    :rtype: numpy array
    :raises: IOERROR if filesize and matrix description appear to be inconsistent

    **BRUKER manual D.14 File Formats (Raw data files)**

        The reconstruction assumes that the raw data to be reconstructed will
        be located in a file named either fid or ser in the EXPNO directory.
        Any averaging of the raw data must be done before the data is written
        to the file, no post-acquisition averag- ing is done by the
        reconstruction.

        K-format
          The raw data file must be written in K-format to be processed by the TOPSPIN
          processing software.  This means that every profile must be written
          to the file at a posi- tion which is a multiple of 1K-bytes (1024 bytes).
          Datasets which are imported may need to be reformatted. Both little-endian
          and big-endian word formats are sup- ported, the BYTORDA parameter in the
          acqp parameter file is assumed to correctly specify the format of the
          raw data file. Finally, the raw data file is assumed to contain data only,
          without any headers.

        Packed format
          Data to be processed with ParaVision can also be stored in packed format.
          For this purpose the parameter GO_block_size must be set to continuous.

        Integer formats
          By default, raw data is stored in 32-bit signed integer format, real and
          imaginary points of the complex signal interleaved. For applications where
          the data range may not exceed 16-bit, e.g. when accumulation and oversampling
          are disabled in analog mode, data can be stored as 16-bit short integers,
          setting the parameter GO_raw_data_format to GO_16BIT_SGN_INT. A warning
          message is displayed if an overflow occurs during data acquisition.

        Floating point format
          Alternatively, data can be stored in floating point format:
          GO_raw_data_format = GO_32BIT_FLOAT. Note, that data not acquired in the
          default 32-bit signed integer format cannot be processed with TOPSPIN but
          only by the ParaVision reconstruction.

        Order of data storage
          The raw data file usually contains the data in the order of acquisition.
          There is only one exception: in case a pipeline filter (AU) is used to
          process the raw data, the order of data storage is defined by the output
          order of the filter.

          For raw datasets acquired by ParaVision methods, ordering and size of the
          data are described by the ACQP parameters contained in the file acqp:

           * ACQ_dim - Spatial/spectroscopic dimension of the experiment.
           * ACQ_size[ ] - Array of length ACQ_dim with length of the given dimensions.
             ACQ_size[0] counts real valued data and for this it must be an even number.
             Real and imaginary data is stored shuffled. For experiments acquired with
             multiple receivers (ACQ_experiment_mode == ParallelExperiment) the scans
             from different receivers are appended. The number of active receivers is
             derived from the number of selected receivers in the parameter array
             GO_ReceiverSelect.
           * NI - Number of objects.
             e.g. slices, but not purely!, could be slices*echoes
           * ACQ_obj_order - Permutation of the order of the acquired objects.
           * ACQ_phase_factor - Number of subsequent scans in the raw dataset belonging
             to the same Object. Typically, this parameter is set to 1 which will keep
             the cod- ing constant within a Multiplex step. If ACQ_phase_factor > 1,
             this parameter will give the number of consecutively acquired phase encoding
             steps that belong to a single Object. Concerning phase increment:
             see ACQ_rare_factor.
           * ACQ_phase_enc_mode[ ] - Ordering scheme for scans within the k-space.

                .. note::  parameter seems to be renamed to ACQ_phase_encoding_mode

           * ACQ_phase_enc_start[ ] - For positioning of first scan in phase encoding
             scheme.
           * ACQ_rare_factor - For positioning of the ACQ_phase_factor scans within the
             k-space. In the case of ACQ_phase_factor > 1 the phase encoding increment
             is determined by ACQ_size[1] / ACQ_rare_factor.
             ACQ_rare_factor = ACQ_phase_factor is used for RARE experiments.
           * NR - Number of repeated experiments within the dataset.
           
    These salient parameters are printed using the following code:
    ['%s: %s' % (k,acqp[k]) for k in ['ACQ_dim','ACQ_size','GO_block_size','NI','ACQ_obj_order','ACQ_phase_factor','ACQ_phase_encoding_mode','ACQ_phase_enc_start','ACQ_rare_factor','NR','NSLICES','ACQ_spatial_phase_1','ACQ_spatial_size_1']]

    Examples:

        >>> import os
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/1/fid'))
        >>> fid['data'].shape   # TriPilot multi
        (128, 128, 15)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/2/fid'))
        >>> fid['data'].shape   # FLASH 2D
        (133, 105, 5)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/3/fid'))
        >>> fid['data'].shape   # FLASH 3D
        (133, 105, 25)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/4/fid'))
        >>> fid['data'].shape   # MSME 2D
        (133, 105, 5)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/5/fid'))
        >>> fid['data'].shape   # MSME 3D
        (133, 105, 25)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/6/fid'))
        >>> fid['data'].shape   # MSME 2D-TURBO
        (256, 256, 5, 3)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/7/fid'))
        >>> fid['data'].shape  # FLASH 2D (NR=25, NI=5, NSLICES=5)
        (133, 105, 5, 25)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/8/fid'))
        >>> fid['data'].shape  # FLASH 2D, partial acq. NR auto reset to 5
        (133, 105, 5, 5)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/9/fid')) # doctest:+ELLIPSIS
        Traceback (most recent call last):
        ...
        OSError: ...
        >>> # fid file 9 was missing due to incomplete scans
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/10/fid'))
        >>> fid['data'].shape # FLASH 2D (MATRIX 32 X 32)
        (32, 32, 5)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/11/fid'))
        >>> fid['data'].shape # FLASH 3D (MATRIX 32 X 32)
        (32, 32, 5)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/12/fid'))
        ... # doctest: +ELLIPSIS
        ... # 1-segment EPI - FIXME, this should be easy but somehow ACQ_size=(8192,1) 
        ... # and ACQ_scan_size=ACQ_phase_factor_scans (instead of One_Scan)
        ... # this should be an easy 64x64x5slice single shot EPI...
        Traceback (most recent call last):
        ...
        ValueError: ...
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/13/fid'))
        ... # doctest: +ELLIPSIS
        ... # 16 segment EPI - FIXME
        Traceback (most recent call last):
        ...
        ValueError: ...
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/14/fid'))
        >>> fid['data'].shape # DTI Standard
        (133, 105, 5, 2)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/15/fid'))
        ... # doctest: +ELLIPSIS
        ... # DTI SPIRAL - FIXME
        Traceback (most recent call last):
        ...
        TypeError: ...
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/16/fid'))
        >>> fid['data'].shape # UTE 2D
        (64, 402, 5)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/17/fid'))
        >>> fid['data'].shape # UTE 3D
        (64, 51360)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/18/fid'))
        >>> fid['data'].shape # ZTE 3D
        (512, 51896)
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/99/fid'))
        ... # doctest: +ELLIPSIS
        ... # interrupted FLASH DCE, mismatch in filesize causes IOError
        ... # this is desired (and expected) behaviour
        Traceback (most recent call last):
        ...
        IOError: ...
        >>> fid = readfid(os.path.expanduser('~/bdata/stefan/nmr/readfidTest.ix1/99/fid'), resetNR=True)
        >>> fid['data'].shape # interrupted FLASH DCE (NR=7 from formerly 25)
        (133, 105, 5, 7)
    """
    if isinstance(fptr, FileType):
        fidname = fptr.name
    if isinstance(fptr, (StringType,UnicodeType)):
        fidname = fptr
    dirname = os.path.abspath(os.path.dirname(fidname))

    # use parameter files provided by caller or load if needed
    acqp = acqp or readJCAMP(os.path.join(dirname,'acqp'))
    method = method or readJCAMP(os.path.join(dirname,'method'))

    if "Spectroscopic" in acqp['ACQ_dim_desc']:
        raise TypeError(
            "Problem: Could this be a spectro scan instead of an image?\n"+
            "If so, use readfidspectro()")

    # determine data type
    if acqp['GO_raw_data_format'] == 'GO_32BIT_SGN_INT':
        datatype = 'i4'
    elif acqp['GO_raw_data_format'] == 'GO_16BIT_SGN_INT':
        datatype = 'i2'
    elif acqp['GO_raw_data_format'] == 'GO_32BIT_FLOAT':
        datatype = 'f4'
    else:
        raise IOError('Unknown ##$GO_raw_data_format = '\
                          + acqp['GO_raw_data_format'])
    dtype = numpy.dtype(datatype)

    # Spatial/spectroscopic dimension of the experiment.
    logger.debug('ACQ_dim={0}'.format(acqp['ACQ_dim']))
    # only the following parameters are retrieved for easier access AND
    # modification from their value in the acqp parameter file
    #ACQ_size[0] should be even (real valued counts)
    ACQ_size = acqp['ACQ_size'][:]
    NR = acqp['NR']
    logger.debug('ACQ_size={0}'.format(ACQ_size))
#    assert acqp['ACQ_experiment_mode'] == 'SingleExperiment',(
#            'I am not clever enough to read Parallel acquired data, yet')
    assert acqp['ACQ_dim'] == len(ACQ_size),(
            'ACQ_dim = {0} != len(ACQ_size={1}) ??'.format(
            acqp['ACQ_dim'], ACQ_size))

    # There is the possibility of 'zero filling' since objects (aka kspace
    # lines) are written n 1Kb blocks if GO_block_size != continuous
    logger.debug('GO_block_size={0}'.format(acqp['GO_block_size']))
    if acqp['GO_block_size'] == 'continuous':
        # we acquire complex data which requires numbers in the read direction
        ACQ_size[0] /= 2
    elif acqp['GO_block_size'] == 'Standard_KBlock_Format':
        true_obj_blocksize = (int(datatype[1])*ACQ_size[0])
        obj_blocksize = 1024 * ((true_obj_blocksize-1) // 1024 + 1)
        ACQ_size[0] = obj_blocksize//(2*int(datatype[1]))
    else:
        raise IOError, 'Unexpected value for GO_block_size in acqp'

    if acqp['ACQ_dim'] == 2:
        # this is 2D
        encoding = [1, 1, 0, 0] # dimensions that require FFT
    else:
        encoding = [1, 1, 1, 0] # dimensions that require FFT

    #find the sequence of k-space lines
    # this information should be stored in the PVM_EncSteps1
    # parameter of the methods file. It appears to be missing from
    # EPI data
    try:
        phase_range = max(acqp['ACQ_spatial_phase_1']) - \
                        min(acqp['ACQ_spatial_phase_1'])
        Enc1Steps = (numpy.array(acqp['ACQ_spatial_phase_1']) -
                    min(acqp['ACQ_spatial_phase_1']))
        Enc1Steps *= (acqp['ACQ_spatial_size_1']-1)/phase_range
        # rounding of float to nearest int is a tricky business...
        Enc1Steps = numpy.floor(Enc1Steps+.5).astype('int')
    except KeyError:
        logger.info('ACQ_spatial_phase_1 not found in acqp, '+
                'trying PVM_EncSteps1 from method')
        try:
            PVM_EncSteps1 = method['PVM_EncSteps1']
            # ensure that it runs from 0 to max
            Enc1Steps = numpy.array(PVM_EncSteps1)-min(PVM_EncSteps1)
        except KeyError:
            logger.warning('PVM_EncSteps1 missing from method parameter file')
            Enc1Steps = numpy.arange(ACQ_size[1])
    if acqp['ACQ_dim'] == 3:
        try:
            phase_range = max(acqp['ACQ_spatial_phase_2']) - \
                            min(acqp['ACQ_spatial_phase_2'])
            Enc2Steps = (numpy.array(acqp['ACQ_spatial_phase_2']) -
                        min(acqp['ACQ_spatial_phase_2']))
            Enc2Steps *= (acqp['ACQ_spatial_size_2']-1)/phase_range
            # rounding of float to nearest int is a tricky business...
            Enc2Steps = numpy.floor(Enc2Steps+.5).astype('int')
        except KeyError:
            logger.info('ACQ_spatial_phase_2 not found in acqp, '+
                    'trying PVM_EncSteps2 from method')
            try:
                PVM_EncSteps1 = method['PVM_EncSteps2']
                print('='*80+'\nsave by PVM_EncSteps2')
                # ensure that it runs from 0 to max
                Enc2Steps = numpy.array(PVM_EncSteps1)-min(PVM_EncSteps1)
            except KeyError:
                logger.warning('PVM_EncSteps2 missing from method parameter file')
                Enc2Steps = numpy.arange(ACQ_size[2])
    else:
        Enc2Steps = [0]

    n_stored_datapoints = numpy.array(ACQ_size).prod() * (
                acqp['NSLICES']*acqp['ACQ_n_echo_images']*
                acqp['ACQ_n_movie_frames'] * NR )
    n_datapoints = numpy.array(acqp['ACQ_size']).prod() * (
                acqp['NSLICES']*acqp['ACQ_n_echo_images']*
                acqp['ACQ_n_movie_frames'] * NR )

    fid_size =  n_stored_datapoints * int(datatype[1]) * 2 
    file_size = os.stat(fidname).st_size
    # load data
    logger.info('loading %s' % fidname)
    if (fid_size < file_size):
        logger.warning(
            ' %s: filesize (%i) > expected, calculated size (%i) ' % (
                 fidname, file_size, fid_size) +
            '\nwill truncate fid ...'
             )
        data = numpy.fromfile(fptr, dtype = dtype)[0:(2*n_stored_datapoints)]
    elif (fid_size>file_size):
        if resetNR: 
            NR = int(file_size//(int(datatype[1])*2*numpy.array(ACQ_size).prod() *
                               acqp['NSLICES']*acqp['ACQ_n_echo_images']*
                               acqp['ACQ_n_movie_frames']))
            n_stored_datapoints = numpy.array(ACQ_size).prod() * (
                acqp['NSLICES']*acqp['ACQ_n_echo_images']*
                acqp['ACQ_n_movie_frames'] * NR )
            n_datapoints = numpy.array(acqp['ACQ_size']).prod() * (
                acqp['NSLICES']*acqp['ACQ_n_echo_images']*
                acqp['ACQ_n_movie_frames'] * NR )

            data = numpy.fromfile(fptr, dtype = dtype)[0:(2*n_stored_datapoints)]
        else:
            raise IOError('filesize (%i) < expected size of fid (%i)' % 
                          (file_size, fid_size))
    else:
        data = numpy.fromfile(fptr, dtype = dtype)

    # byteorder: 'little' don't swap, 'big' do swap
    logger.debug('BYTORDA={0}'.format(acqp['BYTORDA']))
    if acqp['BYTORDA'] =='big':
        data.byteswap(True)  # swap ENDIANness in place
        
    # convert to complex data, usually we would do:
    #    >>>> fid = data[::2]+1j*data[1::2]
    # the following is faster by a factor of 6 to 7!!!
    fid = data.astype(numpy.float32).view(numpy.complex64)

    logger.info('ACQ_size = {0}, NR={1}, ACQ_obj_order={2}, Enc1Steps={3}'.
                   format(ACQ_size, NR, acqp['ACQ_obj_order'], Enc1Steps))

    if untouched:
        return {'data':fid,
                'isImage':False,
                'header':{'acqp': acqp, 'method': method}}
    else:
# ==============================================================================
#       below code has been discussed on stackoverflow
#       http://stackoverflow.com/questions/5422184/numpy-efficient-execution-of-a-complex-reshape-of-an-array
#       it turned out that vectorizing the index building (numpy.vectorize)
#       combined  with array-based array access was actually giving some
#       speed improvements. The slowest part is the last assignment in the
#       for loop but this can be improved by reshaping to use array indices
#       directly.  At 1.5s per (128, 64, 1, 1200) we are stuck with what we
#       have and would have to resort to cython or worse to improve
# ==============================================================================

        # reshape into a large 2D array with dimensions

        fid = fid.reshape(n_stored_datapoints/ACQ_size[0], ACQ_size[0])
        # Now is the time to lop off any spurious zero-filling added on disk
        # due to  parameter acqp['GO_block_size'] == 'Standard_KBlock_Format'
        # If this is set, blocks of 1k bytes make up the read lines (and might
        # be zerofilled if not long enough...)
        tempfid = fid[:,0:(acqp['ACQ_size'][0]/2)]

        # ACQ_size - might 1, 2, or thre elements depnding on ACQ_dim
        # append 1 to make it 3D
        ACQ_size = numpy.hstack([acqp['ACQ_size'], 
                                 numpy.tile(1,3-acqp['ACQ_dim'])])
        
        fid_reorder = numpy.ndarray((ACQ_size[0]/2, 
                        ACQ_size[1], 
                        ACQ_size[2], 
                        acqp['NSLICES'], # caveat: in 3D this is usually =1
                        acqp['ACQ_n_echo_images'], # for multi-echos, diffusion?
                        acqp['ACQ_n_movie_frames'], # could be LL step, diffusion?
                        NR # repetition of an experiment
                        ),dtype = 'complex')

        # this idx enumerates all the kspace lines (scans in BRUKER pulse 
        # programming lingo)
        idx = numpy.arange(n_datapoints/ACQ_size[0], dtype=int)
                             
        # the following assignments create lists that map any of index values
        # in idx to some value in the dimension corresponding to that specialized
        # index. E.g. PE1_idx takes the Enc1Steps (usually a linear int function
        # from 1..nr of enc steps) and copies them out a sufficient number of times
        # until all idx elements have a corresponding entry in PE1_idx.
        PE1_idx = numpy.tile(numpy.repeat(Enc1Steps,
                                          acqp['NSLICES']*
                                          acqp['ACQ_n_echo_images']*
                                          acqp['ACQ_n_movie_frames']),
                             n_datapoints/(ACQ_size[0]*ACQ_size[1]*
                                           acqp['NSLICES']*
                                           acqp['ACQ_n_echo_images']*
                                           acqp['ACQ_n_movie_frames']))
        echo_sequence = numpy.array(acqp['ACQ_obj_order'][
                                                0:acqp['ACQ_n_echo_images']])
        echo_idx = numpy.tile(echo_sequence,
                              n_datapoints/(ACQ_size[0]*acqp['ACQ_n_echo_images']))

        PE2_idx = numpy.tile(numpy.repeat(Enc2Steps,
                                          n_datapoints/(NR*ACQ_size[0]*ACQ_size[2])),
                             NR)
        #ACQ_obj_order lists echoes and slices and movie frames!
        slice_sequence = numpy.array(acqp['ACQ_obj_order'][
                                     :(acqp['NSLICES']*acqp['ACQ_n_echo_images'])
                                     :acqp['ACQ_n_echo_images']],
                                    dtype=int)//(acqp['ACQ_n_echo_images'])
        slice_idx = numpy.tile(numpy.tile(slice_sequence,
                                          n_datapoints/(NR*ACQ_size[0]*
                                                        acqp['NSLICES'])),
                               NR)
        mov_idx = numpy.tile(numpy.repeat(numpy.arange(acqp['ACQ_n_movie_frames']), 
                                           acqp['NSLICES']),
                              n_datapoints/(ACQ_size[0]*
                                            acqp['NSLICES']*
                                            acqp['ACQ_n_movie_frames']))
                               
                             
        NR_idx = numpy.repeat(numpy.arange(NR), n_datapoints/(ACQ_size[0]*NR))

        # The next assignment causes the inherent requiremet for all these
        # indices PE1_idx, PE2_idx etc to be of the same length as idx.
        fid_reorder[:, PE1_idx, PE2_idx, slice_idx, 
                    echo_idx, mov_idx, NR_idx] = tempfid[idx, :].T

        if squeezed:
            fid_return = numpy.squeeze(fid_reorder)
        else:
            fid_return = fid_reorder

        return {'data':fid_return,
                'isImage':False,
                'header':{'acqp': acqp,
                          'method': method,
                          'encoding': encoding  # indicates the dims that require FFT
                          }
                }

def readfidspectro(fptr=None, 
                   acqp=None,
                   method=None,
                   untouched=False):
    """
    Returns BRUKER's fid file as a properly dimensioned & rearranged array

    :param fptr: filename of fid file or filehandle to open fid file
    :type fptr: string or FileType
    :param untouched: Do not rearrange lines into slices and echos and such
                      in the fid in its form as found on disk
    :type untouched: boolean
    :return: array of kspace data
    :rtype: numpy.array
    :raises: AssertError for various inconsistencies in data size or simply
             misunderstandings of how the data is to be interpreted.
    """

    if isinstance(fptr, FileType):
        fidname = fptr.name
    if isinstance(fptr, StringType):
        fidname = fptr
    dirname = os.path.abspath(os.path.dirname(fidname))
    acqp = acqp or readJCAMP(os.path.join(dirname,'acqp'))
    method = method or readJCAMP(os.path.join(dirname,'method'))

    assert "Spectroscopic" in acqp['ACQ_dim_desc'] ,(
            "Problem: Could this be an imaging instead of a spectro scan?")

    # determine data type
    if acqp['GO_raw_data_format'] == 'GO_32BIT_SGN_INT':
        datatype = 'i4'
    elif acqp['GO_raw_data_format'] == 'GO_16BIT_SGN_INT':
        datatype = 'i2'
    elif acqp['GO_raw_data_format'] == 'GO_32BIT_FLOAT':
        datatype = 'f4'
    else:
        raise IOError('Unknown ##$GO_raw_data_format = '\
                          + acqp['GO_raw_data_format'])
    dtype = numpy.dtype(datatype)
    # byteorder: 'little' don't swap, 'big' do swap
    BYTORDA = acqp['BYTORDA']

    # determine array dimensions
    ACQ_size = acqp['ACQ_size'][:]
    # we acquire complex data which requires numbers in the read direction
    ACQ_size[0] /= 2
    # number of objects (increments? e.g., for inversion recovery?)
    NI = acqp['NI']
    # find BRUKER object order
    ACQ_obj_order = acqp['ACQ_obj_order']
#    ACQ_phase_factor = acqp['ACQ_phase_factor']
#    ACQ_phase_encoding_mode = acqp['ACQ_phase_encoding_mode']
#    ACQ_phase_enc_start = acqp['ACQ_phase_enc_start']
#    ACQ_rare_factor = acqp['ACQ_rare_factor']
    #see ho many repetitions
    NR = acqp['NR']

    dtype = numpy.dtype(datatype)

    # load data
    data = numpy.fromfile(fptr, dtype=dtype)
    if BYTORDA =='big':
        data.byteswap(True)  # swap ENDIANness in place

    # convert to complex data

    #fid = data[::2]+1j*data[1::2]
    # the following is faster by a factor of 6 to 7!!!
    fid = data.astype(numpy.float32).view(numpy.complex64)

    logger.debug('ACQ_size = {0}, NR={1}, ACQ_obj_order={2}'.
                  format(ACQ_size, NR, ACQ_obj_order))

    if untouched:
        return {'data':fid,
                'isImage':False,
                'header':{'acqp': acqp, 'method': method}}
    else:
        # reshape into a large 2D array with dimensions [readsize, nr(objorder)*phase*NR]
        tempfid = fid.reshape(NR, len(ACQ_obj_order), ACQ_size[0])

        assert NI == len(ACQ_obj_order), (
                "I don't understand how the various acqp parameters interact\n\
                workaround: use readfidspectro with option untouched")

        fid_reorder = numpy.empty((ACQ_size[0], len(ACQ_obj_order), NR),
                            dtype = 'complex')

        idx = numpy.arange(len(ACQ_obj_order)*NR)

        # using array-based index access
        ObjNr = idx % len(ACQ_obj_order)

        fid_reorder[:, ObjNr] = tempfid[:, idx].T

        return {'data':fid_reorder.squeeze(), # squeezing might be problematic
                'isImage':False,
                'header':{'acqp': acqp,
                          'method': method,
                          'encoding': None  # indicates the dims that require FFT
                          }
                }

def read2dseq(scandirname,
              visu_pars=None):
    """
    Returns BRUKER's 2dseq file as a properly dimensioned array

    :param string scandirname: filename of the scan directory
    :param dict reco, d3proc, visu_pars:
        parameter files (as dict) that can be provided by caller
        default: None which means they will be loaded by this function
    :return: dictionary with data, and header information. The data is
             an array of BRUKER-reconstructed image data in the respecive proc
             directory.
    :rtype: dict with 'data':numpy array 'header':dict{'recpo':..., 'd3proc':...}
    :raises: IOERROR if directory non-existent

    This relies on numpy's array functionality
    """

    try:
        reco = readJCAMP(os.path.join(scandirname,'reco'))
        RECO_transposition = reco['RECO_transposition']
    except IOError:
        RECO_transposition = 0
    # get relevant information from the visu_pars files
    visu_pars = visu_pars or readJCAMP(os.path.join(scandirname,'visu_pars'))

    # determine ENDIANness and storage type

    if visu_pars['VisuCoreWordType'] =='_16BIT_SGN_INT':
        datatype = 'i2'
    elif visu_pars['VisuCoreWordType'] =='_32BIT_SGN_INT':
        datatype = 'i4'
    elif visu_pars['VisuCoreWordType'] =='_32BIT_FLOAT':
        datatype = 'f'
    else:
        raise IOError('unknown ##$VisuCoreWordType = '+
                        visu_pars['VisuCoreWordType'])

    if visu_pars['VisuCoreByteOrder'] == 'littleEndian':
        datatype = '<'+datatype
    else:
        datatype = '>'+datatype
    dtype = numpy.dtype(datatype)

    # ANDREW YUNG:
    # load data based on the visu_pars file.  first two dimensions are the in-
    # plane pixels. The higher dimensions refer to things like slice, echo,
    # diffusion weighting, which are described by VisuGroupDepVals and
    # VisuFGOrderDesc.  If these parameters do not exist, assume that there are
    # only 3 dimensions and label the 3rd dimension as "frame".  Output a
    # dictionary which contains a) the image data with appropriate data slopes
    # and offsets applied, b) string list describing each dimension, and c)
    # header data containing parameters from the visu_pars and reco files.
    # files.

    # make a copy of this parameter so as not to change the original!
    matrix_size = visu_pars['VisuCoreSize'][:]

    if RECO_transposition == 0:
        dimdesc=['readout','PE1']
    else:
        dimdesc=['PE1','readout']
    dimcomment=['','']

    # if VisuCoreSize has 3 elements, this is a 3D acquisition
    if len(matrix_size)==3:
        dimdesc.append('PE2')
        dimcomment.append('')

    # Determine size and descriptors of frame groups (RECO_size, dimdesc and
    # dimcomment).  VisuFGOrderDesc is a  struct which describes the number of
    # images in the frame group, the type of frame (e.g. FG_SLICE, FG_ECHO),
    # and index ranges for the parameter array VisuGroupDepVals, which denotes
    # the names of parameter arrays which depend on that particular frame
    # group.  If there is a dependent parameter array called VisuFGElemComment
    # for the frame group, it should be stored in dimcomment (e.g. DTI
    # generated procnos have FA, Tensor trace, etc.)

    if 'VisuFGOrderDesc' in visu_pars:
        for v in visu_pars['VisuFGOrderDesc']:
            FGdim = v[0]
            matrix_size.append(FGdim)
            dimdesc.append(v[1])
            depvalstart = v[3]
            depvalend = depvalstart + v[4]
            more_dimcomment = ''
            for depval in range(depvalstart,depvalend):
                if visu_pars['VisuGroupDepVals'][depval][0] == '<VisuFGElemComment>':
                    FGcommentstart=visu_pars['VisuGroupDepVals'][depval][1]
                    fullFGcomments = visu_pars['VisuFGElemComment']
                    more_dimcomment=fullFGcomments[FGcommentstart:
                                                   FGcommentstart+FGdim]
            dimcomment.append(more_dimcomment)

    # extract binary data from 2dseq. For now, format the data shape so all
    # the image frames are lumped together in the 3rd dimension.
    reco_offset = numpy.asarray(visu_pars['VisuCoreDataOffs'])
    assert len(set(reco_offset)) == 1, 'Cannot deal with multiple VisuCoreDataOffs'
    reco_slope = numpy.asarray(visu_pars['VisuCoreDataSlope'])
    assert len(set(reco_slope)) == 1, 'Cannot deal with multiple VisuCoreDataSlope'


    n_frames = numpy.asarray(matrix_size)[2:].prod()
    logger.info('Guessing we have {0} when VisuCoreFrameCount = {1}'.format(
                    n_frames, visu_pars['VisuCoreFrameCount']))
    matrix_size.reverse()

    filename = os.path.join(scandirname,'2dseq')
    logger.info('read2dseq: loading %s' % filename)
    data = numpy.fromfile(file=filename, dtype=dtype)

    logger.info('Dat shape = {0} while matrix_size={1}'.format(
                data.shape, matrix_size))
    data=data.reshape(n_frames, matrix_size[-2],matrix_size[-1]).astype('float64')

    # now apply the data slopes and offsets to transform the stored binary
    # number into a absolute number
    data = reco_offset[0] + reco_slope[0]*data

    # Finally, shape the data so all frames are put into separate dimensions.
    data = data.reshape(matrix_size)
    # there are two kinds of transposition needed:
    #  (a) transpose so that time, z, y, x -> x, y, z, time
    #  (b) account for the row-major preference in python and for column
    #      major in paravision -> this equates to an in-plane transpose
    #      (swapping x for y)
    swp_axis = range(len(matrix_size)-1)
    swp_axis.reverse()
    swp_axis.insert(1,len(matrix_size)-1)
    data = data.transpose(swp_axis)
    # and now we have to apply thesame logic to te axis descriptors
    # which have already been built from revere accomplishing the
    # transposition (a) above. watch the beautiful pythn index swap in
    # action. a,b = b,a ... genius!
    dimcomment[0], dimcomment[1] = dimcomment[1], dimcomment[0]
    dimdesc[0], dimdesc[1] = dimdesc[1], dimdesc[0]

    return {'data':data,
            'dimcomment':dimcomment,
            'dimdesc':dimdesc,
            'header':{'visu_pars': visu_pars}}


def dict2string(d):
    '''
    convert dictionary to nicely looking multi-line string

    :param dict d: input dictionary
    :return: list of strings
    :rtype: list

    This might be useful when turning the JCAMP-style dictionaries
    into something that goes into a text display. Example use would be::

        >>> d={'TE':2.3, 'TR':5, 'NAME':'random name'}
        >>> print(dict2string(d))
                          TE : 2.3
                          TR : 5
                        NAME : 'random name'
    '''
    strlist = []
    for k, v in d.iteritems():
        strlist.append('{0!s:>20} : {1!s}'.format(k, repr(v)))
    return '\n'.join(strlist)

def fftbruker(array, encoding=None, DCoffset=False):
    '''
    wrapper to fft bruker FIDs

    returns the fft of a multi-dimensional BRUKER FID. It uses the
    parameter 'encoding' to decide over which dimensions to Fourier transform.
    Typically, a 2D only needs FFT over 1st and 2nd dimension
    (encoding = [1, 1, 0, 0]) whereas 3D files get a FT over three dimensions
    (encoding = [1, 1, 1, 0]). The 4th dimension (repetitions) doesn't
    usually get FTed.
    '''

    encoding = encoding or [1, 1, 0, 0]

    #find all axes that should be FTed
    FTaxes = numpy.where(numpy.array(encoding) != 0)[0]

    img = numpy.fft.fftn(numpy.fft.fftshift(array, axes = FTaxes), axes=FTaxes)
    if DCoffset:
        # remove DC offset
        if encoding[2]: # do this only for central point in 3D
            for i in range(img.shape[3]):
                img[0, 0, 0, i] = numpy.median(img[0:2, 0:2, 0:2, i])
        else: # do this for every slice separately
            for i in range(img.shape[3]):
                for j in range(img.shape[2]):
                    img[0, 0, j, i] = numpy.median(img[0:2, 0:2, j, i])

    return numpy.fft.fftshift(img, axes = FTaxes)

def fftfid(fptr=None,
           readfidresult=None,
           AxisFlip=True, 
           PEshift=True,
           **kwargs):
    ''' Take filename, retrieve fid and FT as best as you can
    
    This relies on the FT-able axis to be in he 1st three positions.
    
    :param fptr: input fid file (handed through to readfid)
    :param AxisFlip: Flip axis to achieve display in agreement with radiological
                    convention (default: True)
    :param PEshift: find out shifts of slice packs and apply in the 
                    PE directions (default: True)
    :param squeezed: squeeze FT array before returning (default: true)
    :param **kwargs: handed through to readfid (only squeezed is handled 
                    internally)
    :return: complex-valued array from FT of fid
    :rtype: ndarray
    
    Remarks: It does not account for the possibility of several slice packs
    with different orientations etc. As a result, tripilot scans are not 
    treated properly.'''
    if kwargs.pop('squeezed',True):
        # squeezed was set!
        squeezed=True
    else:
        squeezed=False
    if readfidresult is None:
        readfidresult = readfid(fptr, squeezed=False, **kwargs)
    fid = readfidresult['data']
    if PEshift:
        method = readfidresult['header']['method']
        # Shifting along the phase encodes is not implemented through any modification
        # of the acquisition. Instead, one has to roll the arrays along those
        # dimensions during reconstruction. The amount of rolling (in pixels) 
        # depends on the amount of shift and the FoV in that direction.
        xyz_shift_mm = numpy.array([method['PVM_SPackArrReadOffset'][0]]+
                                [method['PVM_SPackArrPhase1Offset'][0]]+
                                [method['PVM_SPackArrPhase2Offset'][0]], 
                                dtype=float).flatten()
        
        xyz_pix = fid.shape[0:3] # always 3D
        # make 3D even for 2D acq
        xyz_FoV_mm = numpy.hstack([method['PVM_Fov'],
                                   numpy.tile(1,3-len(method['PVM_Fov']))])
        xyz_shift_fraction = xyz_shift_mm/xyz_FoV_mm    
        # create meshed arrays that have linearly increasing values along the 1st
        # and 2nd dimension
        PE1 = numpy.meshgrid(numpy.zeros(xyz_pix[0]),
                             numpy.arange(xyz_pix[1]),
                             numpy.zeros(xyz_pix[2]),
                             numpy.zeros(fid.shape[3]),
                             numpy.zeros(fid.shape[4]),
                             numpy.zeros(fid.shape[5]),
                             numpy.zeros(fid.shape[6]), indexing='ij')[1]
        pe1shift = numpy.exp(complex(0,-2*numpy.pi)*PE1*xyz_shift_fraction[1])
        fid *= pe1shift
        if xyz_pix[2]>1:
            PE2 = numpy.meshgrid(numpy.zeros(xyz_pix[0]),
                                 numpy.zeros(xyz_pix[1]),
                                 numpy.arange(xyz_pix[2]),
                                 numpy.zeros(fid.shape[3]),
                                 numpy.zeros(fid.shape[4]),
                                 numpy.zeros(fid.shape[5]),
                                 numpy.zeros(fid.shape[6]), indexing='ij')[2]
            # echo position 
            echopos = method.get('PVM_EchoPosition',None)
            if echopos is None:
                echopos = 0.5
            else:
                echopos /=100.
            pe2shift = numpy.exp(complex(0,-2*numpy.pi)*PE2*
                                     (xyz_shift_fraction[2]-
                                      (echopos*xyz_pix[2]-1)/xyz_pix[2]))
            fid *= pe2shift
    # shifting along axes that require FT by half the size removes the 
    # horrible phase ripples
    fidshifted = numpy.fft.fftshift(fid, axes=[0,1,2])
    # Performing the FT can always be done along 3 dimensions! This is because
    # the readfid routine slots read, phase and slice encode into array 
    # dimensions 0, 1, 2. Multi-slice fids have the slices in the 4th
    # dimension...
    ii = numpy.fft.fftshift(numpy.fft.fftn(fidshifted,
                                           axes=[0,1,2]),
                            axes=[0,1])
    # The 1st, 2nd, 3rd direction always contains the Read, phase, and slice.
    # Typically, for display purposes we like AP, LR, HF direction along
    # certain directions depending on the image orientation:
    # axial updown[0] = AP, leftright[1] = LR
    if AxisFlip:
        method = readfidresult['header']['method']
        if method['PVM_SPackArrSliceOrient'][0]=='axial':
            if method['PVM_SPackArrReadOrient'][0]=='L_R':
                ii = numpy.transpose(ii, [1,0,2,3,4,5,6])
        # sag: updown[0]= HF, leftright[1] = AP
        elif method['PVM_SPackArrSliceOrient'][0]=='sagittal':
            if method['PVM_SPackArrReadOrient'][0]=='A_P':
                ii = numpy.transpose(ii, [1,0,2,3,4,5,6])
        # cor: updown[0]= HF, leftright[1] = LR                           
        elif method['PVM_SPackArrSliceOrient'][0]=='coronal':
            if method['PVM_SPackArrReadOrient'][0]=='L_R':
                ii = numpy.transpose(ii, [1,0,2,3,4,5,6])
        # flip along all PE & FE dimension
        retarr = ii[::-1,::-1,::-1,...]
    else:
        retarr = ii
    if squeezed:
        return numpy.squeeze(retarr)
    else:
        return retarr

def readRFshape(filename):
    '''reads BRUKERS RF shape files spnam0 etc that are stored with experiments

    returns a dictionary with ready-to-use attributes for amplitude,
    phase, bandwidth etc.'''
    RFstringdict = readJCAMP(filename)

##TITLE= /camelot/PV.2.0.Devel/exp/stan/nmr/lists/wave/imaghermite
##JCAMP-DX= 5.00 BRUKER JCAMP library
##DATA TYPE= Shape Data
##ORIGIN= BRUKER MEDICAL GMBH
##OWNER= <dwe>
##DATE= 98/04/29
##TIME= 09:18:57
##MINX= 0.000000e+00
##MAXX= 1.000000e+02
##MINY= 0.000000e+00
##MAXY= 1.800000e+02
##$SHAPE_EXMODE= Excitation
##$SHAPE_TOTROT= 90.000000e+00
##$SHAPE_BWFAC= 5.4000000e+00
##$SHAPE_INTEGFAC= 1.794e-01
##$SHAPE_REPHFAC = 50
##$SHAPE_TYPE = conventional
##$SHAPE_MODE= 0
##NPOINTS= 1024
##XYPOINTS= (XY..XY)

    RFshape = {}
    for tag in RFstringdict.keys():
        try:
            RFshape[tag] = float(RFstringdict[tag])
        except ValueError:
            RFshape[tag] = RFstringdict[tag]

    # make a consecutive list of float values stored in the XYPOINTS field
    XYlist = ' '.join(RFshape['XYPOINTS'].split(', ')).split()
    RFshape['amp'] = [float(dummy) for dummy in XYlist[1::2]]
    RFshape['phase'] = [float(dummy) for dummy in XYlist[2::2]]

    return RFshape

# main part - run test cases if called as a module
if __name__ == "__main__":
    import doctest
    doctest.testmod()
#    doctest.run_docstring_examples(readfid,globals())
