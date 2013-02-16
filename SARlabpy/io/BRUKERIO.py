#!/usr/bin/env python

# Copyright (C) 2012-2013 Stefan A Reinsberg and SARlab members
# full license details see LICENSE.txt
# Random comment added by FM
"""Collection of BRUKER input routines

Handy functions to read BRUKER data and header files.
x`"""

import numpy
import os.path
from types import StringType, FileType

DEBUG = 1

def readJCAMP(filename, removebrackets=True):
    """
    Parse text file in JCAMP format

    :param filename: filename of fid file
    :type filename: string
    :param removebrackets: format strings without extra brackets?
    :type removebrackets: boolean
    :return: Dictionary of labelled data records (LDR) with LDR-names as keys and their content as values.
    :rtype: dict
    :raises: IOERROR if opening file fails or passes on any other error

    The *JCAMP format* is a self-documenting, ASCII text file format
    that is maintained by IUPAC (`see a report
    here <http://www.iupac.org/objID/Article/pac7108x1549>`_).
    It consists of labelled data records (LDR) that start with a ##
    and end when the next record begins.
    They can span several lines. The data label is
    enclosed between '##' and '='. If it starts with a $, we are
    dealing with a private LDR.

    The issue of reading these is complicated due to the various
    types of data (integers, floats, strings, arrays and nested
    structures) that can be present. Currently no attempt is made to
    perform type conversion before returning a dictionary of the JCAMP
    file.

    """
    import re
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

        # with this dictionary, find all the values that contain a
        # array index at the start of the value, e.g. "( 16 )"
        #for k, v in LDRdict.iteritems():
            #matchdim = re.match(r"\([^\)]+\)", v)
            #if matchdim: # we have a match for dimension
                #remainder = v[matchdim.end():]

        return LDRdict

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

    """
    if isinstance(fptr, FileType):
        fidname = fptr.name
    if isinstance(fptr, StringType):
        fidname = fptr
    dirname = os.path.abspath(os.path.dirname(fidname))
    print(dirname)
    acqp = readJCAMP(dirname + "/acqp")

    # determine array dimensions
    ACQ_size = acqp['ACQ_size'].split() # matrix size
    ACQ_size = [int(dummy) for dummy in ACQ_size]
    # we acquire complex data which requires numbers in the read direction
    ACQ_size[0] /= 2
    #see ho many repetitions
    NR = int(acqp['NR'])
    # find BRUKER object order
    ACQ_obj_order = acqp['ACQ_obj_order'].split()
    ACQ_obj_order = [int(dummy) for dummy in ACQ_obj_order]

    # 2D vs 3D issues about slices etc.
    if len(ACQ_size) == 2:
        # this is 2D
        if len(ACQ_obj_order) != int(acqp['NSLICES']):
            raise IOError, 'NSLICES not equal to number of ACQ_obj_order'
        else:
            ACQ_size.append(len(ACQ_obj_order))
        encoding = [1, 1, 0, 0] # dimensions that require FFT
    else:
        encoding = [1, 1, 1, 0] # dimensions that require FFT


    # RARE factor sort of indicates how many PE shots are spend on one slice
    # before moving on to next object (slice?) as defined in AQ_obj_order
    ACQ_rare_factor = int(acqp['ACQ_rare_factor'])
    # determine data type
    if acqp['GO_raw_data_format'] == 'GO_32BIT_SGN_INT':
        datatype = 'i4'
    else:
        raise IOError('Unknown ##$GO_raw_data_format = '\
                          + acqp['GO_raw_data_format'])
    # not sure about byteorder !?
    dtype = numpy.dtype(datatype)

    #find the sequence of k-space lines
    method = readJCAMP(dirname + '/method')
    PVM_EncSteps1 = method['PVM_EncSteps1'].split()
    PVM_EncSteps1 = [int(dummy) for dummy in PVM_EncSteps1]
    #ensure that it runs from 0 to max
    PVM_EncSteps1 = numpy.array(PVM_EncSteps1)-min(PVM_EncSteps1)

    # load data
    data = numpy.fromfile(fptr, dtype = dtype)

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
    acqp = readJCAMP(dirname+"/acqp")

    assert "Spectroscopic" in acqp['ACQ_dim_desc'] , "Problem: Could this be an imaging isntead of a spectro scan?"

    # determine array dimensions
    ACQ_size = acqp['ACQ_size'].split() # matrix size
    ACQ_size = [int(dummy) for dummy in ACQ_size]
    # we acquire complex data which requires numbers in the read direction
    ACQ_size[0] /= 2
    #see ho many repetitions
    NR = int(acqp['NR'])

    # number of increments are used, e.g., for inversion recovery
    NI = int(acqp['NI'])

    # find BRUKER object order
    ACQ_obj_order = acqp['ACQ_obj_order'].split()
    ACQ_obj_order = [int(dummy) for dummy in ACQ_obj_order]

    # determine data type
    if acqp['GO_raw_data_format'] == 'GO_32BIT_SGN_INT':
        datatype = 'i4'
    else:
        raise IOError('Unknown ##$GO_raw_data_format = '\
                          +acqp['GO_raw_data_format'])
    # not sure about byteorder !?
    dtype = numpy.dtype(datatype)

    #load the method file
    method = readJCAMP(dirname+'/method')
    PVM_EncSteps1 = method['PVM_EncSteps1'].split()
    PVM_EncSteps1 = [int(dummy) for dummy in PVM_EncSteps1]
    #ensure that it runs from 0 to max
    PVM_EncSteps1 = numpy.array(PVM_EncSteps1)-min(PVM_EncSteps1)

    # load data
    data = numpy.fromfile(fptr, dtype=dtype)

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

        assert NI == len(ACQ_obj_order), "I don't understand how the various acqp parameters interact\n\
        Need more coding for correct reco\n\
        workaround: use readfidspectro with option untouched"

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


def read2dseq(procdirname):
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
    RECO_size = [int(dummy) for dummy in (d3proc['IM_SIX'], \
                                          d3proc['IM_SIY'], \
                                          d3proc['IM_SIZ'])]

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
    return {'data':data,
            'isImage':True,
            'header':{'reco': reco, 'd3proc':d3proc}}

def dict2string(d):
    '''
    convert dictionary to nicely looking multi-line string

    :param d: input dictionary
    :type: dict
    :return: list of strings
    :rtype: list
    
    this might be useful when turning the JCAMP-style dictionaries
    into something that goes into a text display
    '''
    strlist = []
    for k, v in d.iteritems():
        strlist.append('%-20s:\t%s' % (k, repr(v)))
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
    print "test case for readJCAMP"
    readJCAMP('/home/stefan/data/HPGS3/HPGS3Tb1.4s1/3/acqp')
    print "test case for read2dseq"
    x = read2dseq('/home/stefan/data/HPGS3/HPGS3Tb1.4s1/3')

