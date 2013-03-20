"""Collection of BRUKER input routines

Handy functions to read BRUKER data and header files.
Stefan A Reinsberg - Spring 2011, Vancouver
Modified by Andrew Yung - Spring 2013
"""
DEBUG=0

import numpy as np
import time
import os.path
from types import StringType, FileType, ListType, FloatType, IntType


# XXXX readBrukerParx
def readBrukerParx(filename):
    """parse Bruker parameter file (JCAMP format) and return as dictionary
    
    The *JCAMP format* is a self-documenting, ASCII text file format that is 
    maintained by IUPAC (http://www.iupac.org/objID/Article/pac7108x1549).
    It consists of labelled data records (LDR) that start with a ## and end
    when the next record begins. They can span several lines. The data label is
    enclosed between '##' and '='. If it starts with a $, we are dealing with a
    private LDR.
    
    This is a version of Stefan's readJCAMP function.  The only difference is
    that this function parses the parameter values into their data type (it
    doesn't just leave the parameter values as strings).  The only exception is
    if a structure array is detected; the function only extracts the elements
    of this structure array as a string.
    
    """
    import re
    import sys
    
    try:
        print('opening ',filename)
        JCAMPfile = open(filename, "r")
        JCAMPdata = JCAMPfile.read().splitlines() # read and lose the "\n"
        JCAMPfile.close()

    except IOError as (errno, strerror):
#        print "There was an I/O error({0}): {1}".format(errno, strerror)
        raise IOError, strerror
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise
    
    else:

       # let's loop through the file, remove comments and put each
        # LDR on its own line 
        LDRlist = [] # start with empty list
        for line in JCAMPdata:

            # match location comment: it begins with $$ and the beginning of path definition
            if re.match("\\$\\$ /.*", line):  #match pattern = "$$ /"
                line = [re.sub("\\$\\$ ", "##$FILE_LOCATION=", line)]

            # match date comment (assumes there is nothing else in the comments)
            elif re.match("\\$\\$ [^@].*", line):  #match pattern = "$$ " but without @ character 
                line = re.sub("\\$\\$ ", "##$DATE=", line)
                uname="##$USERNAME="+(line.split(' '))[-1]
                line = re.sub("[^ ]+$", "", line) # remove username
                line=[line, uname]

            # match all other comments and discard
            elif re.match("\\$\\$ @.*", line):
                line=[]

            # match lists that are normal LDRs
            elif re.match("##.*",line):
                # we should remove parenthesis
                 line=[line]

            # this must be a line that belongs to the preceeding LDR
            # (like for strings or arrays)
            # attach this to the line that was previously appended to 
            # the LDRlist
            else:
                LDRlist[-1]=LDRlist[-1]+" "+line
                line=[]

            #add this to the list of LDRs
            LDRlist.extend(line)

        # strip every labelled data record of its preceding ## or ##$. 

        LDRlist = [re.sub('##[\\$]*','',LDR) for LDR in LDRlist]
 
        # split every LDR at the "=" sign and turn it into a dictionary entry
        LDRdict = dict([LDR.split("=") for LDR in LDRlist])
        
        # enclose the DATE string in <> (there are parentheses that mess up
        # the struct comprehension)
        LDRdict['DATE'] = '< ' + LDRdict['DATE'] + '>'
        
        #######################################################################
        # DEFINE REGULAR EXPRESSIONS TO DENOTE TEXT REPRESENTATION OF LDR VALUES
        #######################################################################

        # Bruker sometimes puts  a string in <> brackets 
        # especially for data labels.  Spaces within an element are allowed       
        enclosedstringpattern = r'<.*>'
    
        # Bruker puts a struct (combination of numbers and strings) into (), 
        # with at least one comma delimiter
        structpattern = r'\(.*,.*\)'        
        
        # integer:  optional + or - followed by one or more digits
        integerpattern = r'[+\-]?\d+'
    
        # real number with scientific notation
        real_scinot_pattern = r'[+\-]?\d(\.\d+|)[Ee][+\-]\d\d?\d?'
        # real number with decimal notation
        real_decimal_pattern = r'[+\-]?(\d+\.\d*|\d*\.\d+)'
        # real number with either notation    
        real_pattern = r'(' + real_scinot_pattern + '|' + real_decimal_pattern + ')'

        # for each parameter entry, do the conversion from string to final data type        
        for k,v in LDRdict.iteritems():
            # First, check if the parameter value is composed of structs '(struct1) (struct2) ...'
            if re.search(structpattern, v):
                # strip off the () and convert string into a list of strings
                # split function adds empty strings so delete them
                v = re.split(r'(?:\) \(|^\s\(|\)$)',v)
                v = [element for element in v if element]
            # if not a struct array, first check if the parameter value is a <>-enclosed string
            elif re.search(enclosedstringpattern,v):
                
                # strip off the <> brackets
                v = re.split(r'(?:\s*<\s*|> <|\s*>\s*)',v)
                v = [element for element in v if element]
            # if not struct array or enclosed string, it must be a list of one or more numbers 
            else:
                # it is now safe to use whitespace as the delimeter between numbers
                v = re.split(r' ',v)
                v = [element for element in v if element]
                # convert to into or float, according to what kind of number was detected
                for index, element in enumerate(v):
                    if re.match(real_pattern, element):
                        v[index] = float(element)
                    elif re.match(integerpattern, element):
                        v[index] = int(element)
                        
            # now that the value string was converted to the appropriate data
            # type, store it in the dictionary.  If there is only one element,
            # convert it to a single number
            if len(v) == 1:
                LDRdict[k] = v[0]
            else:
                LDRdict[k] = v
                    
        return LDRdict  






#XXXXX read2dseq XXXXXXXXXXXXXXXXXXXXXXXXX
def read2dseq_visu(procdirname):
    """returns BRUKER's 2dseq file as a properly dimensioned array 

    This relies on numpy's array functionality
    """
    import re

    # get relevant information from Bruker parameter files
    recopath = os.path.join(procdirname,"reco")
    visupath = os.path.join(procdirname,"visu_pars")
    
    # some procnos do not have reco files (e.g. DTI calculated maps)
    reco_found = True
    try:
        reco=readBrukerParx(recopath)
    except:
        reco_found = False
        reco = []
        

    # determine the readout direction   
    if reco_found:
       RECO_transposition = reco['RECO_transposition']
    else:
       RECO_transposition = 0

    visu=readBrukerParx(visupath)

    # determine ENDIANness and storage type
    if visu['VisuCoreWordType']=='_16BIT_SGN_INT':
        datatype='i2'
    elif visu['VisuCoreWordType']=='_32BIT_SGN_INT':
        datatype='i4'
    elif visu['VisuCoreWordType']=='_32BIT_FLOAT':
        datatype='f'
    else:
        raise IOError('unknown ##VisuCoreWordType='+visu['VisuCoreWordType'])

    if visu['VisuCoreByteOrder'] == 'littleEndian':
        datatype='<'+datatype
    else:
        datatype='>'+datatype

    dtype=np.dtype(datatype)

    # load data based on the visu_pars file.  first two dimensions are the in-               
    # plane pixels. The higher dimensions refer to things like slice, echo, 
    # diffusion weighting, which are described by VisuGroupDepVals and
    # VisuFGOrderDesc.  If these parameters do not exist, assume that there are     
    # only 3 dimensions and label the 3rd dimension as "frame".  Output a
    # dictionary which contains a) the image data with appropriate data slopes
    # and offsets applied, b) string list describing each dimension, and c) 
    # header data containing parameters from the visu_pars and reco files.
    # files.

    RECO_size=visu['VisuCoreSize']
    if RECO_transposition == 0:
        dimdesc=['readout','PE1']
    else:
        dimdesc=['PE1','readout']
    dimcomment=[[''],['']]    
 
    # if VisuCoreSize has 3 elements, this is a 3D acquisition
    if len(RECO_size)==3:
        dimdesc.append('PE2')
        dimcomment.append([''])

    # Determine size and descriptors of frame groups (RECO_size, dimdesc and 
    # dimcomment).  VisuFGOrderDesc is a  struct which describes the number of
    # images in the frame group, the type of frame (e.g. FG_SLICE, FG_ECHO),
    # and index ranges for the parameter array VisuGroupDepVals, which denotes
    # the names of parameter arrays which depend on that particular frame 
    # group.  If there is a dependent parameter array called VisuFGElemComment
    # for the frame group, it should be stored in dimcomment (e.g. DTI 
    # generated procnos have FA, Tensor trace, etc.)
    
    if 'VisuFGOrderDesc' in visu:
        structpattern = r'^(\d+),\s*<(.*)>,\s*<(.*)>,\s*(\d+),\s*(\d+)$'
        for v in visu['VisuFGOrderDesc']:
            match=re.search(structpattern,v)
            if match:
               FGdim = int(match.group(1))
               RECO_size.append(FGdim)
               dimdesc.append(match.group(2))
               depvalstart = int(match.group(4))
               depvalend = depvalstart + int(match.group(5))
               GroupDepVals = visu['VisuGroupDepVals']
               for depval in range(depvalstart,depvalend):
                   commentmatch=re.search('<VisuFGElemComment>, (\d+)', GroupDepVals[depval])
                   if commentmatch:
                       FGcommentstart=int(commentmatch.group(1))
                       fullFGcomments = visu['VisuFGElemComment'];
                       dimcomment.append(fullFGcomments[FGcommentstart:FGcommentstart+FGdim])
               if not commentmatch:
                   dimcomment.append([''])

    # extract binary data from 2dseq. For now, format the data shape so all
    # the image frames are lumped together in the 3rd dimension.
    reco_offset = np.asarray(visu['VisuCoreDataOffs'])
    reco_slope = np.asarray(visu['VisuCoreDataSlope'])    
    n_frames = np.asarray(RECO_size)[2:].prod()
    datapath=os.path.join(procdirname,'2dseq')
    data=np.fromfile(file=datapath, dtype=dtype)
    data=data.reshape(RECO_size[0],RECO_size[1],n_frames)

    # now apply the data slopes and offsets to transform the stored binary 
    # number into a absolute number
    if reco_slope.size == 1:  #single slope value which applies to all frames
        data = reco_offset + reco_slope*data
    else:  # separate slope + offset for each frame
        for i in range(0,n_frames):
            data[:,:,i] = reco_offset[i] + reco_slope[i]*data[:,:,i]
    
    # Finally, shape the data so all frames are put into separate dimensions.
    data = data.reshape(RECO_size)    
    
    return {'data':data, 
            'isImage':True,
            'dimcomment':dimcomment,
            'dimdesc':dimdesc,
            'header':{'visu': visu, 'reco': reco}}




    
# XXXX main XXXX
# main part - run test cases if called as a module
if __name__ == "__main__":
   import pprint 
    print "test case for readBrukerParx"
#    acqp_dict=readBrukerParx('I:\\Documents\\Dropbox\\T1VFAphant.3U1\\34\\acqp')
    visu_dict=readBrukerParx('I:\\Documents\\Dropbox\\KDVPM101.dP1\\10\\pdata\\2\\visu_pars')
    pprint.PrettyPrinter().pprint(visu_dict)

    print "test case for read2dseq_visu"
    x=read2dseq_visu('I:\\Documents\\Dropbox\\T1VFAphant.3U1\\34\\pdata\\1')
#    x=read2dseq_visu('I:\\Documents\\Dropbox\\KDVPM101.dP1\\10\\pdata\\2')

    
