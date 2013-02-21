#!/usr/bin/env python

# Copyright (C) 2012-2013 Stefan A Reinsberg and SARlab members
# full license details see LICENSE.txt
"""Collection of BRUKER input routines

Handy functions to read BRUKER data and header files.
"""

import numpy
import os.path
import re
from types import StringType, FileType

DEBUG = 0

def readJCAMP(filename, removebrackets=True, typecast=False):
    """
    Parse text file in JCAMP format

    :param string filename: filename of fid file
    :param boolean removebrackets: format strings without extra brackets?
    :param boolean typecast:
        attempt to cast values of records to int/floats/strings/
        a list of the above (Default: False) - this feature is a tad
        experimental
    :return: Dictionary of labelled data records (LDR) with LDR-names as keys and their content as values.
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
    structures) that can be present. A currently experimental feature is the
    typecasting of the records into all these different datatypes.
    """
    import sys

    try:
        if DEBUG >= 1:
            print "opening {0}".format(filename)
        JCAMPfile = open(filename, "r")
        JCAMPdata = JCAMPfile.read().splitlines() # read and lose the "\n"
        JCAMPfile.close()

    except IOError as (errno, strerror):
        print "There was an I/O error({0}): {1}".format(errno, strerror)
        raise IOError, strerror
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise

    else:

        # let's loop through the file, remove comments and put each
        # LDR on its own line
        LDRlist = [] # start with empty list
        for line in JCAMPdata:

            # match location comment
            if re.match("\\$\\$ \\/.*", line):
                line = [re.sub("\\$\\$ ", "##$FILE_LOCATION=", line)]

            # match date comment (assumes there is nothing else in the comments)
            elif re.match("\\$\\$ [^@].*", line):
                line = re.sub("\\$\\$ ", "##$DATE=", line)
                uname = "##$USERNAME="+(line.split(' '))[-1]
                line = re.sub("[^ ]+$", "", line) # remove username
                line = [line, uname]

            # match all other comments
            elif re.match("\\$\\$ @.*", line):
                line = []

            # match lists that are normal LDRs
            elif re.match("##.*", line):
                # we should remove parenthesis
                if removebrackets:
                    line = re.sub("\\(.*\\)", "", line)
                line = [line]

            # this must be a line that belongs to the preceeding LDR
            # attach this to the line that was previously appended to
            # the LDRlist
            else:
                LDRlist[-1] = LDRlist[-1] + " " + line
                line = []

            #add this to the list of LDRs
            LDRlist.extend(line)

        # use python list comprehension to strip every labelled
        # data record of its preceding ## or ##$.
        LDRlist = [re.sub('##[\\$]*', '', LDR) for LDR in LDRlist]

        # use python list comprehension to split every LDR at the
        # " = " sign and turn it into a dictionary entry
        LDRdict = dict([LDR.split("=") for LDR in LDRlist])

        if typecast:
            return typecast_dict_elements(LDRdict)
        else:
            return LDRdict

def typecast_dict_elements(hetero_dict):
    '''
    Return a dictionary of previously only string values in a typecast form

    :param dict hetero_dict:
        input dictionary of heteregenous type (int, float, etc.)
    :return: dictionary with values cast to int, float etc
    :rtype: dict

    Here are some example calls:
        
        >>> typecast_dict_elements({'key':'1'})
        {'key': 1}
        >>> typecast_dict_elements({'key':'1.31'})
        {'key': 1.31}
        >>> typecast_dict_elements({'key':'Spam'})
        {'key': 'Spam'}
        >>> typecast_dict_elements({'key':'1.31 123'})
        {'key': [1.31, 123.0]}
        >>> typecast_dict_elements({'key':'1.31 123 string'})
        {'key': '1.31 123 string'}
        >>> typecast_dict_elements({'key':'(1, 2, <>) (3, 4, <>)'})
        {'key': [['1', ' 2', ' <>'], ['3', ' 4', ' <>']]}

    '''
    for dict_item in hetero_dict:
        try:
            # is it just an int?
            hetero_dict[dict_item] = int(hetero_dict[dict_item])
        except ValueError:
            try:
                #maybe it is just a float then?
                hetero_dict[dict_item] = float(hetero_dict[dict_item])
            except ValueError:
                # no, then this is currently a string in one of three flavours
                # Case 1: array of ints. 'ACQ_phase_enc_start': ' -1 -1',
                # Case 2: array of floats. 'ACQ_spatial_phase_1': ' -1 -0.9583 -0.9166 -0.875 -0.8333 ...
                # Case 3: weird combo. 'TPQQ': ' (<hermite.exc>, 16.4645986123031, 0) (<fermi.exc>, 115.8030276379, 0)
                split_string = [s for s in re.split(' ',
                                                    hetero_dict[dict_item]) if s]
                try:
                    #is this a homogeneous list of ints?
                    hetero_dict[dict_item] = [int(s) for
                                          s in split_string]
                except ValueError:
                    try:
                        #no? maybe it is a homogeneous list of floats?
                        hetero_dict[dict_item] = [float(s) for
                                              s in split_string]
                        #this was a homogeneous list of floats
                    except ValueError:
                        # Case 3: Array of floats, ints, and strings with funny brackets and such
                        # Example: 'TPQQ': ' (<hermite.exc>, 16.4645986123031, 0) (<fermi.exc>, 115.8030276379, 0)
                        list_level_one = [s.strip(' ()')
                                for s in re.split('\) *\(', hetero_dict[dict_item])]
                        if len(list_level_one)>1:
                            hetero_dict[dict_item] = [list_element.split(',') for
                                    list_element in list_level_one]
                        else:
                            try:
                                # does the RH unpack to  single-element list?
                                (hetero_dict[dict_item],) = list_level_one[0].split(',')
                            except ValueError:
                                # no? OK, let's keep the whole list
                                hetero_dict[dict_item] = list_level_one[0].split(',')
                    

    return hetero_dict

def readfid(fptr=None, untouched=False):
    """
    Returns BRUKER's fid file as a properly dimensioned & rearranged array.

    :param fptr: filename of fid file or filehandle to open fid file
    :type fptr: string or FileType
    :param untouched: leave the fid as found without rearranging lines into slices and echos?
    :type untouched: boolean
    :return: Flat (untouched = True) or Rearranged and assembled array of the acquire k-space
    :rtype: numpy array
    :raises: IOERROR if filesize and matrix description appear to be inconsistent

    >>> import os
    >>> fname = os.path.expanduser('~/data/Moosvi.ii1/3/fid')
    >>> fid = readfid(fname)
    /Users/stefan/data/Moosvi.ii1/3
    >>> fid['data'].shape
    (128, 128, 3, 1)
    
    This is what the manual D.14 File Formats (Raw data files)

The reconstruction assumes that the raw data to be reconstructed will be 
located in a file named either fid or ser in the EXPNO directory. Any averaging 
of the raw data must be done before the data is written to the file, no 
post-acquisition averag- ing is done by the reconstruction.

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
   * ACQ_obj_order - Permutation of the order of the acquired objects. 
   * ACQ_phase_factor - Number of subsequent scans in the raw dataset belonging
     to the same Object. Typically, this parameter is set to 1 which will keep 
     the cod- ing constant within a Multiplex step. If ACQ_phase_factor > 1, 
     this parameter will give the number of consecutively acquired phase encoding 
     steps that belong to a single Object. Concerning phase increment: 
     see ACQ_rare_factor.
   * ACQ_phase_enc_mode[ ] - Ordering scheme for scans within the k-space. 
   * ACQ_phase_enc_start[ ] - For positioning of first scan in phase encoding 
     scheme. 
   * ACQ_rare_factor - For positioning of the ACQ_phase_factor scans within the
     k-space. In the case of ACQ_phase_factor > 1 the phase encoding increment 
     is determined by ACQ_size[1] / ACQ_rare_factor. 
     ACQ_rare_factor = ACQ_phase_factor is used for RARE experiments. 
   * NR - Number of repeated experiments within the dataset.


    """
    if isinstance(fptr, FileType):
        fidname = fptr.name
    if isinstance(fptr, StringType):
        fidname = fptr
    dirname = os.path.abspath(os.path.dirname(fidname))
    print(dirname)
    acqp = readJCAMP(dirname + "/acqp", typecast=True)
    
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
    ACQ_dim = acqp['ACQ_dim'] #Spatial/spectroscopic dimension
    ACQ_size = acqp['ACQ_size'] #ACQ_size[0] should be even (real valued counts)
    assert acqp['ACQ_experiment_mode'] == 'SingleExperiment',(
            'I am not clever enough to read Parallel acquired data, yet')
    assert ACQ_dim == len(ACQ_size),(
            'Not enough/too much information about the size of length in\n'+ 
            'each dim imae dimension')

    
    #There is the possibility of 'zero filling' since objects (aka kspace 
    # lines) are written n 1Kb blocks if GO_block_size != continuous
    if acqp['GO_block_size'] == 'continuous':
        # we acquire complex data which requires numbers in the read direction
        ACQ_size[0] /= 2
    elif acqp['GO_block_size'] == 'Standard_KBlock_Format':
        true_obj_blocksize = (int(datatype[1])*ACQ_size[0])
        obj_blocksize = 1024 * ((true_obj_blocksize-1) // 1024 + 1)
        ACQ_size[0] = obj_blocksize/2/int(datatype[1])
    else:
        raise IOError, 'Unexpected value for GO_block_size in acqp'
        
    #see ho many repetitions
    NR = acqp['NR']
    # find BRUKER object order
    ACQ_obj_order = acqp['ACQ_obj_order']
    ACQ_phase_factor = acqp['ACQ_phase_factor']
    ACQ_phase_enc_mode = acqp['ACQ_phase_enc_mode']
    ACQ_phase_enc_start = acqp['ACQ_phase_enc_start']
    ACQ_rare_factor = acqp['ACQ_rare_factor']

    # 2D vs 3D issues about slices etc.
    if ACQ_dim == 2:
        # this is 2D
        assert len(ACQ_obj_order) == acqp['NSLICES'],(
                'NSLICES not equal to number of ACQ_obj_order')

        ACQ_size.append(len(ACQ_obj_order))
        encoding = [1, 1, 0, 0] # dimensions that require FFT
    else:
        encoding = [1, 1, 1, 0] # dimensions that require FFT


    # RARE factor sort of indicates how many PE shots are spend on one slice
    # before moving on to next object (slice?) as defined in AQ_obj_order
    ACQ_rare_factor = acqp['ACQ_rare_factor']
    #find the sequence of k-space lines
    method = readJCAMP(dirname + '/method', typecast=True)
    try:
        PVM_EncSteps1 = method['PVM_EncSteps1']
        #ensure that it runs from 0 to max
        PVM_EncSteps1 = numpy.array(PVM_EncSteps1)-min(PVM_EncSteps1)
    except KeyError:
        print('Warning: PVM_EncSteps1 missing from method parameter file')
        PVM_EncSteps1 = numpy.arange(ACQ_size[1])

    # load data
    data = numpy.fromfile(fptr, dtype = dtype)
    if BYTORDA =='big':
        data.byteswap(True)  # swap ENDIANness in place

    # convert to complex data
    #fid = data[::2]+1j*data[1::2]
    # the following is faster by a factor of 6 to 7!!!
    fid = data.astype(numpy.float32).view(numpy.complex64)

    if DEBUG >= 1:
        print('ACQ_size = {0}, NR={1}, ACQ_obj_order={2}, EncSteps={3}'.
            format(ACQ_size, NR, ACQ_obj_order, PVM_EncSteps1))

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
        # [readsize, nr(objorder)*phase*NR]
        if DEBUG >= 1:
            print('size(fid) = {0} ?=? ACQ_size*NR={1}'.
                  format(numpy.shape(fid), ACQ_size[0] *
                           ACQ_size[1] * ACQ_size[2] * NR))
        if numpy.size(fid) < (ACQ_size[0] * ACQ_size[1] * ACQ_size[2] * NR):
            print ('''WARNING:    parameter file describes a file that is
            larger than what is actually found in */fid''')
            NR = numpy.size(fid) / (ACQ_size[0]*ACQ_size[1]*ACQ_size[2])
            print ('            NR reset to {0}'.format(NR))
            fid = fid[0:ACQ_size[0]*ACQ_size[1]*ACQ_size[2]*NR]

        tempfid = fid.reshape(ACQ_size[1]*ACQ_size[2]*NR, ACQ_size[0])
        fid_reorder = numpy.empty((ACQ_size[0], ACQ_size[1],
                                   len(ACQ_obj_order), NR),
                            dtype = 'complex')

        idx = numpy.arange(ACQ_size[1] * len(ACQ_obj_order) * NR)

# = =============================================================================
#         for i in idx:
#             # work out where the rearranged matrix gets its data from
#             slicenr = ACQ_obj_order[(i/ACQ_rare_factor) % len(ACQ_obj_order)]
#             # linear encoding: PEnr = (i/len(ACQ_obj_order)) % ACQ_size[1]
#             # a more general scheme taking into account RARE PE is trickier:
#             PEidx = ((i % (ACQ_size[1]*len(ACQ_obj_order))) /
#                     (ACQ_rare_factor*len(ACQ_obj_order))*ACQ_rare_factor
#                  + (i % ACQ_rare_factor))
#             PEnr = PVM_EncSteps1[PEidx]
#             REPnr = i / (ACQ_size[1]*len(ACQ_obj_order))
#             print(i, ':', PEnr, slicenr, REPnr)
#             fid_reorder[:, PEnr, slicenr, REPnr] = tempfid[i, :]
# = =============================================================================

        # alternative using array-based index access
        # We are vectorizing the functions so they can eb c
        PEnr = numpy.array(PVM_EncSteps1)[(idx % (ACQ_size[1]*len(ACQ_obj_order))) /
                    (ACQ_rare_factor*len(ACQ_obj_order))*ACQ_rare_factor
                 + (idx % ACQ_rare_factor)]
        slicenr = numpy.array(ACQ_obj_order)[(idx/ACQ_rare_factor) % len(ACQ_obj_order)]

        REPnr = idx / (ACQ_size[1] * len(ACQ_obj_order))

        fid_reorder[:, PEnr, slicenr, REPnr] = tempfid[idx, :].T

        return {'data':fid_reorder,
                'isImage':False,
                'header':{'acqp': acqp,
                          'method': method,
                          'encoding': encoding  # indicates the dims that require FFT
                          }
                }

def readfidspectro(fptr=None, untouched=False):
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
    print(dirname)
    acqp = readJCAMP(dirname+"/acqp", typecast=True)

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
    ACQ_size = acqp['ACQ_size']
    # we acquire complex data which requires numbers in the read direction
    ACQ_size[0] /= 2
    # number of objects (increments? e.g., for inversion recovery?)
    NI = acqp['NI']
    # find BRUKER object order
    ACQ_obj_order = acqp['ACQ_obj_order']
    ACQ_phase_factor = acqp['ACQ_phase_factor']
    ACQ_phase_enc_mode = acqp['ACQ_phase_enc_mode']
    ACQ_phase_enc_start = acqp['ACQ_phase_enc_start']
    ACQ_rare_factor = acqp['ACQ_rare_factor']
    #see ho many repetitions
    NR = acqp['NR']

    dtype = numpy.dtype(datatype)

    #load the method file
    method = readJCAMP(dirname+'/method')
    PVM_EncSteps1 = method['PVM_EncSteps1'].split()
    #ensure that it runs from 0 to max
    PVM_EncSteps1 = numpy.array(PVM_EncSteps1)-min(PVM_EncSteps1)

    # load data
    data = numpy.fromfile(fptr, dtype=dtype)
    if BYTORDA =='big':
        data.byteswap(True)  # swap ENDIANness in place

    # convert to complex data

    #fid = data[::2]+1j*data[1::2]
    # the following is faster by a factor of 6 to 7!!!
    fid = data.astype(numpy.float32).view(numpy.complex64)

    if DEBUG >= 1:
        print('ACQ_size = {0}, NR={1}, ACQ_obj_order={2}, EncSteps={3}'.
            format(ACQ_size, NR, ACQ_obj_order, PVM_EncSteps1))

    if untouched:
        return {'data':fid,
                'isImage':False,
                'header':{'acqp': acqp, 'method': method}}
    else:

        if DEBUG >= 1:
            print(ACQ_size, ACQ_obj_order, NR)
            print(fid.shape)

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

        print idx
        print ObjNr

        fid_reorder[:, ObjNr] = tempfid[:, idx].T

        return {'data':fid_reorder.squeeze(), # squeezing might be problematic
                'isImage':False,
                'header':{'acqp': acqp,
                          'method': method,
                          'encoding': None  # indicates the dims that require FFT
                          }
                }

def read2dseq(procdirname, arrayTranspose = False):
    """
    Returns BRUKER's 2dseq file as a properly dimensioned array

    :param procdirname: filename of directory that contains the processed data
    :type procdirname: string
    :return: dictionary with data, and headerinformation. The data is
             an array of BRUKER-reconstructed image data in the respecive proc
             directory.
    :rtype: dict with 'data':numpy array 'header':dict{'recpo':..., 'd3proc':...}
    :raises: IOERROR if directory non-existent

    This relies on numpy's array functionality
    """

    # get relevant information from the  reco file

    reco = readJCAMP(procdirname+'/reco')
    d3proc = readJCAMP(procdirname+'/d3proc')

    # determine array dimensions
    RECO_size = [int(dummy) for dummy in (d3proc['IM_SIZ'], \
                                          d3proc['IM_SIY'], \
                                          d3proc['IM_SIX'])]

    # determine ENDIANness and storage type
    if reco['RECO_wordtype'] =='_16BIT_SGN_INT':
        datatype = 'i2'
    elif reco['RECO_wordtype'] =='_32BIT_SGN_INT':
        datatype = 'i4'
    elif reco['RECO_wordtype'] =='_32BIT_FLOAT':
        datatype = 'f'
    else:
        raise IOError('unknown ##$RECO_wordtype = '+reco['RECO_wordtype'])

    if reco['RECO_byte_order'] == 'littleEndian':
        datatype = '<'+datatype
    else:
        datatype = '>'+datatype

    dtype = numpy.dtype(datatype)

    # load data
    data = numpy.fromfile(file=procdirname+'/2dseq',
                          dtype=dtype).reshape(RECO_size)


    if arrayTranspose:
        print 'Changing data to result in XYZ intead of ZYX'
        data = data.transpose((2,1,0))


    return {'data':data,
            'isImage':True,
            'header':{'reco': reco, 'd3proc':d3proc}}

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

def fftbruker(array, encoding=[1, 1, 0, 0], DCoffset=False):
    '''
    wrapper to fft bruker FIDs

    returns the fft of a multi-dimensional BRUKER FID. It uses the
    parameter 'encoding' to decide over which dimensions to Fourier transform.
    Typically, a 2D only needs FFT over 1st and 2nd dimension
    (encoding = [1, 1, 0, 0]) whereas 3D files get a FT over three dimensions
    (encoding = [1, 1, 1, 0]). The 4th dimension (repetitions) doesn't
    usually get FTed.
    '''

    #find all axes that should be FTed
    FTaxes = []
    _encoding = encoding[:]
    for i in range(len(_encoding)):
        if _encoding.pop():
            FTaxes.append(len(_encoding))

    img = numpy.fft.fftn(array, axes=FTaxes)
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



def readRFshape(filename):
    '''reads BRUKERS RF shape files spnam0 etc that are stored with experiments

    returns a dictionary with ready-to-use attributes for amplitude,
    phase, bandwidth etc.'''
    try:
        RFstringdict = readJCAMP(filename)
    except IOError:
        print('problem reading RF file shape')
        return None

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
    RFshape['amp'] = [float(dummy) for dummy in XYlist[::2]]
    RFshape['phase'] = [float(dummy) for dummy in XYlist[1::2]]

    return RFshape

# main part - run test cases if called as a module
if __name__ == "__main__":
    import os
    fname = os.path.expanduser('~/data/Moosvi.ii1/3/fid')
    fid = readfid(fname)
    import doctest
    doctest.testmod()
