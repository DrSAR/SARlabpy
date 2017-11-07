# -*- coding: utf-8 -*-
"""
Class Definitions for BRUKER data
"""
import os
import re
import glob
from datetime import datetime

import nibabel
import sys
import collections

import configobj
import pandas

from ..helpers import natural_sort
from .lowlevel import (readJCAMP, readfid, readfidspectro, read2dseq,
                       fftfid,
                       dataroot)
from .AData_classes import AData, load_AData

from .lazy_property import lazy_property

from . import JCAMP_comparison

## From http://docs.python.org/2/howto/logging.html#logging-basic-tutorial
import logging
logger=logging.getLogger(__name__)

masterlist_root = os.path.expanduser('~/sdata/masterlists')

class AttrDict(object):
    '''a class that behaves like a dict but gives you access to its attributes
    directly'''
   
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
    def __getitem__(self, key):
        return self.__dict__[key]
    def __setitem__(self, key, value):
        self.__dict__[key] = value
    def __delitem__(self, key):
        del self.__dict__[key]
    def __contains__(self, key):
        return key in self.__dict__
    def __len__(self):
        return len(self.__dict__)
    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        name_str = self.__class__.__name__+'(\n{}\n)'
        return name_str.format(',\n'.join(args))
    def keys(self):
        return self.__dict__.keys()
    def iteritems(self):
        return self.__dict__.iteritems()
    def values(self):
        return self.__dict__.values()


def strip_all_but_classname(obj, class_str):
    '''
    Nicely format the class_str that is used in some of the __str__ calls
    '''
    m = re.search("'[^']+'", str(obj.__class__))
    if m:
        class_name = m.group(0).strip("'")
    else:
        class_name = class_str
    return class_name

def last_path_components(absdirname, depth=1):
    '''
    Return the last two parts of a path. Should be platform independent and
    fairly immune to trailing '/' and other irregularities. Used to determine
    the shortdirname in Scan and Study.

    >>> last_path_components('~/bdata/stefan/nmr/readfidTest.ix1/9')
    '9'
    >>> last_path_components('~/bdata/stefan/nmr/readfidTest.ix1/9/')
    '9'
    >>> last_path_components('~/bdata/stefan/nmr/readfidTest.ix1/9/', depth=2)
    'readfidTest.ix1/9'
    '''
    head = absdirname.rstrip(os.sep)
    rval = []
    for i in range(depth):
        head, tail = os.path.split(head)
        rval.insert(0,tail)
    return os.sep.join(rval)

masterlist_lookup = None
def find_all_patients_in_masterlists():
    global masterlist_lookup
    if masterlist_lookup is None:
        masterlist_lookup={}
        config_file_names = os.listdir(masterlist_root)
        for config_file_name in config_file_names:
            fname = os.path.join(masterlist_root, config_file_name)
            conf = configobj.ConfigObj(fname)
            for k in list(conf.keys()):
                if k != 'General':
                    assert masterlist_lookup.get(k) is None, \
                        'Non-unique section name {0} in masterlist {1} (previously in {2})'.format(k, fname, masterlist_lookup[k])
                    masterlist_lookup[k]=fname
    else:
        raise RuntimeError("Config Files have already been perused.")

# ===========================================================

class JCAMP_file(AttrDict):
    '''
    Represents a JCAMP encoded parameter file.

    Parameters become attributes in this class
    '''
    def __init__(self, filename):
        '''read subject (JCAMP) file and fill dictionary using the 
        __init__ routine of parent class'''
        if not os.path.isfile(filename):
            raise IOError('File "%s" not found' % filename)
        self.filename = filename
        super(JCAMP_file, self).__init__(**readJCAMP(self.filename))

class PDATA_file(object):
    '''
    Initialize a processed data set that usually sits in `*/pdata/[1-9]`.

    constructor expects a directory name (filename). It will look
    for the 2dseq file, the d3proc and th reco file.
    '''
    def __init__(self, filename):
        self.filename = filename

    @lazy_property
    def reco(self):
        try:
            reco_obj = JCAMP_file(os.path.join(self.filename,'reco'))
        except IOError:
            logger.warning('reco file for %s not found \n' % self.filename +
                    '(this is normal for some BRUKER sets: e.g. DTI)' )
            return None
        else:
            return reco_obj

    @lazy_property
    def visu_pars(self):
        try:
            return JCAMP_file(os.path.join(self.filename,'visu_pars')) 
        except IOError: 
            logger.warning('visu_pars file %s not found\n' % self.filename) 
            return None

    @lazy_property
    def d3proc(self):
        return JCAMP_file(os.path.join(self.filename,'d3proc'))

    @lazy_property
    def data(self):
        dta = read2dseq(os.path.join(self.filename),
                                 visu_pars=self.visu_pars.__dict__)
        setattr(self, 'dimdesc', dta['dimdesc'])
        setattr(self, 'dimcomment', dta['dimcomment'])
        return dta['data']

    @lazy_property
    def dimdesc(self):
        # need to load the data which we'll force by some parameter access
        self.data.shape
        return self.dimdesc

    @lazy_property
    def dimcomment(self):
        # need to load the data which we'll force by some parameter access
        self.data.shape
        return self.dimcomment


    def uid(self):
        try:
            return self.visu_pars.VisuUid
        except AttributeError:
            return ''

    def store_adata(self, *args, **kwargs):
        '''
        Store some secondary data for this PData scan.

        See call signature of AData_classes.AData.fromdata

        This triggers an update in the adata dictionary of the
        parent object (Scan).
        '''
        raise NotImplementedError('Deprecated - should not be confused with\n'+
                                  'Scan.store_adata()')
        return AData.fromdata(self, *args, **kwargs)

    def export2nii(self,filename,rescale = None, std_mod=None):
        '''
        Export the data originating from a 2dseq (BRUKER) reconstruction to
        a Nifti file format.

        :param string filename: where to write the file

        A lot of the magic happens in visu_pars_2_Nift1Header in order to
        figure out geometry.

        Examples:
            >>> import tempfile
            >>> scn = Scan(os.path.join('readfidTest.ix1','3'))
            >>> fname = os.path.join(tempfile.gettempdir(),'readfid.nii')
            >>> print('writing tempfile %s' % fname) #doctest:+ELLIPSIS
            writing tempfile ...readfid.nii
            >>> scn.pdata[0].export2nii(fname)

        '''
        from .visu_pars_2_Nifti1Header import visu_pars_2_Nifti1Header
        import sarpy.analysis.getters
        header = visu_pars_2_Nifti1Header(self.visu_pars)
        aff = header.get_qform()

        if rescale is None:
            img_pair = nibabel.nifti1.Nifti1Image(self.data, aff, header=header)
            img_pair.to_filename(filename)
        else:
   
            img_pair = nibabel.nifti1.Nifti1Image(self.data, aff, header=header)         

            h = img_pair.get_header()
            
            minVal, maxVal = sarpy.analysis.getters.get_image_clims(self.data, std_mod)
            h['cal_min'] = minVal
            h['cal_max'] = maxVal
            
            img_pair = nibabel.nifti1.Nifti1Image(self.data, aff, header=h)
            img_pair.to_filename(filename)            


class Scan(object):
    '''
    Object to represent a BRUKER scan consisting, typically, of FID_file
    and PDATA_file(s).

    The __init__() method expects the filename to be a directory. It will
    try around a bit (filenames inside the directory) before throwing an
    exception. It will also give up if it can't find any of the
    acqp, method, fid or one 2dseq file.

    Simple Example:

        >>> import os
        >>> exp3 = Scan('readfidTest.ix1/3')
        >>> dir(exp3)    # doctest: +NORMALIZE_WHITESPACE
        ['__class__', '__delattr__', '__dict__', '__doc__', '__format__',
         '__getattribute__', '__hash__', '__init__', '__module__', '__new__',
         '__reduce__', '__reduce_ex__', '__repr__', '__setattr__',
         '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'acqp',
         'adata', 'dirname', 'fid', 'method', 'pdata', 'pdata_uids',
         'shortdirname', 'store_adata']

        >>> exp3.acqp.ACQ_protocol_name
        'FLASH_bas(modified)'
        >>> exp3.acqp.ORIGIN
        'Bruker BioSpin MRI GmbH'
        >>> exp3.acqp.DATE
        'Thu Feb 21 19:01:26 2013 PST (UT-8h)'

    Also, Scan has a pretty way of falling over when no fid file is detected:

        >>> Scan('readfidTest.ix1/9') # doctest:+ELLIPSIS
        Traceback (most recent call last):
            ...
        IOError:...
        
    Scan should be able to read spectro scans:
    
        >>> Scan('VFAT1phant1.mo1/8').fid.shape
        (2048, 30)
        
        >>> Scan('VFAT1phant1.mo1/9').fid.shape
        (2048, 30)
    '''
    @lazy_property
    def acqp(self):
        return JCAMP_file(os.path.join(self.dirname,'acqp'))

    @lazy_property
    def method(self):
        return JCAMP_file(os.path.join(self.dirname,'method'))

    @lazy_property
    def fid(self):
        try:
            kspace = readfid(os.path.join(self.dirname,'fid'),
                                      squeezed=False,
                                      acqp=self.acqp.__dict__,
                                      method=self.method.__dict__)['data']
        except TypeError:
            try:
                kspace = readfidspectro(os.path.join(self.dirname,'fid'),
                                      acqp=self.acqp.__dict__,
                                      method=self.method.__dict__)['data']
            except IOError:
                # is there a ser instead o fid file?
                kspace = readfidspectro(os.path.join(self.dirname,'ser'),
                                          acqp=self.acqp.__dict__,
                                          method=self.method.__dict__)['data']
        return kspace

    @lazy_property
    def fftfid(self):
        readfidresult = {'data': self.fid,
                         'header':{'method':self.method.__dict__}}
        return fftfid(os.path.join(self.dirname,'fid'),
                               readfidresult=readfidresult)

    @lazy_property
    def adata(self):
        adata_dict = load_AData(self.pdata, self.dirname)
        #check whether dependent adata is outdated
        for ad in adata_dict.items():
            if  ad[1].meta['depends_on'] != 'UNKNOWN':
                for dep in ad[1].meta['depends_on']:
                    if dep != '':
                        scn_name, adata_lbl = dep.split(';')
                        t1 = datetime.strptime(ad[1].meta['created_datetime'],'%c')
                        try:
                            t2str = Scan(scn_name).adata[adata_lbl].meta['created_datetime']
                        except KeyError:
                            msg = ('adate["{2}"] depends on "{0}" which '+
                                   'is missing adata "{1}"'
                                   ).format(scn_name, adata_lbl,ad[1].key)
                            logger.warning(msg)
                        except OSError:
                            msg = ('adate["{1}"] depends on "{0}" which '+
                                   'is missing?'
                                   ).format(scn_name, ad[1].key)
                            logger.warning(msg)
                        else:
                            t2 = datetime.strptime(t2str,'%c')
                            if t2>t1:
                                msg = ('dependency "{0};{1}" is younger than '+
                                      'it should be').format(scn_name, adata_lbl)
                                logger.warning(msg)
                                  
        if len(adata_dict) == 0:
            logger.info('No analysis data in directory "{0}"'.
                         format(self.dirname))
        return adata_dict

    @lazy_property
    def pdata(self):
        pdata_list = []
        for f in natural_sort(glob.glob(os.path.join(self.dirname,'pdata','*'))):
            try:
                pdata = PDATA_file(f)
                pdata_list.append(pdata)
            except IOError:
                pass
        return pdata_list

    @lazy_property
    def pdata_uids(self):
        return [pdata.uid() for pdata in self.pdata]

    @lazy_property
    def _masterlist(self):
        global masterlist_lookup
        PatName = re.match('[^.]+', self.shortdirname).group()
        # look for PatName amongst the section defined in any of the masterlists
        # masterlist_lookup gets populated once at the module level
        if masterlist_lookup is None:
            find_all_patients_in_masterlists()
        configfile = masterlist_lookup.get(PatName, None)
        
        #if this is successfull then the dictionary will be accessible as scn.__masterlist
        if configfile is not None:
            conf = configobj.ConfigObj(configfile)
            return conf
        else:
            return None
        
    def __init__(self, root, absolute_root=False):
        '''
        Is the filename a direcotry and can we at least find
        one of acqp, fid, or 2dseq

        :param string filename: directory name for the BRUKER data set
        '''
        global masterlist_lookup
        if masterlist_lookup is None:
            find_all_patients_in_masterlists()

        if absolute_root:
            filename = root
        else:
            extended_search = glob.glob(os.path.join(dataroot,'*','nmr',root))
            if len(extended_search) > 1:
                raise IOError(
                    'Data root (%s) is not unique ' % root +
                    '\nDataset present for several users!')
            elif len(extended_search) == 0:
                raise IOError(
                    'Data root (%s) does not exist' % root +
                    '\nwrong filename?')
            else:
                filename = extended_search[0]
        if os.path.isdir(filename):
            # good, that could be it
            self.dirname = filename
        elif os.path.isdir(os.path.dirname(filename)):
            logger.warning(('Scan initialized with a filename (%s) not a' +
                            ' directory.') % filename)
            # did the punter point the filename to some file inside
            # the directory?
            self.dirname = os.path.dirname(filename)
        else:
            # no good
            raise IOError(
                'Filename "%s" is neither a directory nor a file inside one' %
                filename)

        self.shortdirname = last_path_components(self.dirname, depth=2)
        
        mtch = re.match('[^.]+', self.shortdirname)
        # patientname should match SUBJECT_id as given in file subject for 
        # every Study. However, at the Scan level we have no notion
        # of Studies. Ideally one would assert consistency when Scans
        # are assembled into Studies
        if mtch is not None:
            self.patientname = mtch.group()
        else:
            self.patientname = None
            
        mtch = re.search('\....', self.shortdirname)
        if mtch is not None:
            self.studyname = mtch.group()[1:]
        else:
            self.studyname = None
            
        mtch = re.search('\/[0-9]+', self.shortdirname)
        if mtch is not None:
            self.scannumber = mtch.group()[1:]
            
        self.masterlist_filename = masterlist_lookup.get(self.patientname)

        # see whether we can find an fid file
        # in all likelihood this means that an acqp and method file
        # is also present - this was true for ca 9000 scans we tested thi in
        if not(os.path.isfile(os.path.join(self.dirname,'fid')) or
               os.path.isfile(os.path.join(self.dirname,'ser'))):
            raise IOError(
                ('Directory "{0}" did not contain a BRUKER fid or ser '+
                'file').format(self.dirname))

        # if there are no 2dseq files (processed data) this is potentially
        # solvable of the data is retro-actively reconstructed. In the
        # meantime -> warning
        if len(glob.glob(os.path.join(self.dirname,'pdata','*','2dseq'))) == 0:
            logger.warning(('Directory "{0}" did not contain processed data '+
                         '(2dseq)').format(self.dirname))


    def __repr__(self):
        '''
        Simple print representation of the Scan object:

            >>> a=Scan('readfidTest.ix1/5')
            >>> print(a)
            __main__.Scan("readfidTest.ix1/5")
        '''
        return '{0}("{1}")'.format(
                    strip_all_but_classname(self, 'Scan'),
                    self.shortdirname)

    def __str__(self):
        '''
        More elaborate string representation of the Scan object:

            >>> a=Scan('readfidTest.ix1/5')
            >>> a
            Scan object: Scan("readfidTest.ix1/5")
               protocol: MSME_bas
               TE = [14]
               TR = [50]
               FA = 180
               FoV = [3, 4, 2.5]
               matrix = [266, 105, 25]
        '''
        try:

            try:
                return ('Scan object: Scan("{0}")\n'+
                                            '   protocol: {1}\n'+
                                            '   TE = {2}\n'+
                                            '   TR = {3}\n'+
                                            '   FA = {4}\n'+
                                            '   FoV = {5}\n'+
                                            '   matrix = {6}\n'
                                            '   -----BS-----\n'+
                                            '   BSFreq = {7}\n'+
                                            '   BSPower = {8}').format(self.shortdirname,
                                                    self.acqp.ACQ_protocol_name,
                                                    self.acqp.ACQ_echo_time,
                                                    self.acqp.ACQ_repetition_time,
                                                    self.acqp.ACQ_flip_angle,
                                                    self.acqp.ACQ_fov,
                                                    self.acqp.ACQ_size,
                                                    self.method.BSFreqOffset,
                                                    self.method.BSPulse[3])
            except:
                return ('Scan object: Scan("{0}")\n'+
                        '   protocol: {1}\n'+
                        '   TE = {2}\n'+
                        '   TR = {3}\n'+
                        '   FA = {4}\n'+
                        '   FoV = {5}\n'+
                        '   matrix = {6}').format(self.shortdirname,
                                self.acqp.ACQ_protocol_name,
                                self.acqp.ACQ_echo_time,
                                self.acqp.ACQ_repetition_time,
                                self.acqp.ACQ_flip_angle,
                                self.acqp.ACQ_fov,
                                self.acqp.ACQ_size)
        except AttributeError:
            return self.__str__()

    def store_adata(self, pdata_idx=0, force=False, **kwargs):
        '''
        Store an AData set for one of the processed children of this scan.
        Typically the first one (pdata_idx=0)

        :param boolean force:
            overwrite pre-existing AData sets, (default  False)
        :param boolean compressed:
            compress adata files on write to disk, (default  True)
        '''
        import getpass
        # only get rid of adata when forced and when from our user! 
        if self.adata.get(kwargs['key']) is not None:
            owner = self.adata.get(kwargs['key']).meta.get('username', None)
            if owner != getpass.getuser():
                raise AttributeError(
                            'existing ADATA "%s" belongs to other user (%s)' % 
                            (kwargs['key'], owner))
            elif force:
                self.adata.pop(kwargs['key'], None)
        self.adata[kwargs['key']] = AData.fromdata(
                    self.pdata[pdata_idx], **kwargs)

    def rm_adata(self, key):
        '''
        delete adata set
        '''
        deleted = False
        for k in list(self.adata.keys()):
            if re.search(key+'$', k) is not None:
                if self.adata.pop(k, None) is not None:
                    logger.info('adata %s deleted for %s' %(key, self.shortdirname))
                    deleted = True
        if not deleted:
            logger.warning('adata %s NOT found in %s' %(key, self.shortdirname))
                    
    def masterlist_attr(self, attr, level=None):
        '''look up attr in the hierarchy of dictionaries of the masterlist
        config file. As a rule, the attribute stored at the more specific
        (deeper) level takes precedence.
        
        level = 'Experiment', 'Patient', 'Study', 'Scan', None
        '''
        
        (attr_experiment, attr_patient, attr_study, attr_scan) = (None,
                                                                  None,
                                                                  None,
                                                                  None)
        attr_experiment = self._masterlist.get(attr)
        patient_dic = self._masterlist.get(self.patientname)
        if patient_dic is not None:
            attr_patient = patient_dic.get(attr)
            study_dic = patient_dic.get('study '+self.studyname)
            if study_dic is not None:
                attr_study = study_dic.get(attr)
                scanlabels_dic = study_dic.get('scanlabels')
                if scanlabels_dic is not None:
                    for k,v in scanlabels_dic.items():
                        if v == self.scannumber:
                            scan_dic = study_dic.get(k)
                    if scan_dic is not None:
                        attr_scan = scan_dic.get(attr)
            
        # attributes at all the levels are assembled
        attribs = (attr_experiment, attr_patient, attr_study, attr_scan)
        logger.debug('Attributes found: {0}'.format(attribs))
        # find out which one takes precedence
        level_number = {'Experiment':0,
                        'Patient':1,
                        'Study':2,
                        'Scan':3}.get(level, 3)
        
        # climb up the levels in search for the attribute and see whether
        # any of them (i) is not None; stop when reaching level=0=Experiment
        while (attribs[level_number] is None) and (level_number > 0):
            level_number = level_number-1
        return attribs[level_number]

class Study(object):
    '''
    A study in BRUKER parlance is a collection of scans performed on
    on subject (typically within a day, without interruption). A study
    directory is expected to contain a JCAMP file 'subject' and have one
    or more subdirectories which are scans.
    '''
    def __init__(self, root, absolute_root=False):
        '''
        Initialize the Study through loading metadata from the subject file.

        :param string filename: directory name for the BRUKER study
        '''
        if absolute_root:
            filename = root
        else:
            extended_search = glob.glob(os.path.join(dataroot,'*','nmr',root))
            if len(extended_search) > 1 :
                raise IOError(
                    'Data root (%s) is not unique ' % root +
                    '\nStudy present for several users!')
            elif len(extended_search) < 1:
                raise IOError(
                    'Study root (%s) does not exist ' % root)
            else:
                filename = extended_search[0]

        if os.path.isdir(filename):
            # good, that could be it
            self.dirname = filename
        elif os.path.isdir(os.path.dirname(filename)):
            logger.warning(('Scan initialized with a filename (%s) not a' +
                            ' directory.') % filename)
            # did the punter point the filename to some file inside
            # the directory?
            self.dirname = os.path.dirname(filename)
        else:
            # no good
            raise IOError(
                'Filename "%s" is neither a directory nor a file inside one' %
                filename)

        self.shortdirname = last_path_components(self.dirname)
        self.subject = JCAMP_file(os.path.join(self.dirname,'subject'))

    @lazy_property
    def masterlist_filename(self):
        if self.scans:
            default_name = self.scans[0].masterlist_filename
            for scn in self.scans:
                assert scn.masterlist_filename == default_name, \
                    'Scans within Study described by multiple masterlist files! ({0} != {1})'.format(
                    default_name, scn.masterlist_filename)
            return default_name
        else:
            return None

    @lazy_property
    def _masterlist(self):
        if self.masterlist_filename is None:
            return None
        else:
            return configobj.ConfigObj(self.masterlist_filename)
    
    @lazy_property
    def scans(self):
        scans = []
        eligible_dirs = natural_sort(os.listdir(self.dirname))
        for fname in eligible_dirs:
            filename = os.path.join(self.dirname, fname)
            if (os.path.isdir(filename) and
                re.match('[0-9]+', fname)):
                try:
                    scans.append(Scan(filename))
                except IOError:
                    logger.warning('could not open scan: %s (ignoring)' % filename)
                    pass
        return scans

    def __repr__(self):
        '''
        Simple print representation of the Study object

            >>> a=Study('readfidTest.ix1')
            >>> print(a)
            __main__.Study("readfidTest.ix1")
        '''
        return '{0}("{1}")'.format(
                    strip_all_but_classname(self, 'Study'),
                    self.shortdirname)

    def __str__(self):
        '''
        More elaborate string representation of the Study object:

        >>> a=Study('readfidTest.ix1')
        >>> a
        Study object: Study("readfidTest.ix1")
           subject: Moosvi, readfidTest
           scans: 1-TriPilot-multi
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  MSME_bas
                  MSME_bas
                  MSME_turbo_bas
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  MSME_bas
                  EPI_bas
                  EPI_nav_bas
                  DtiStandard_bas
                  UTE2D
                  UTE3D
                  ZTE
                  FLASH_bas(modified)
        '''
        try:
            scan_list_repr = [scan.acqp.ACQ_protocol_name for
                                scan in self.scans]
            return ('Study object: Study("{0}")\n'+
                    '   subject: {1}\n'+
                    '   scans: {2}').format(
                                self.shortdirname,
                                self.subject.SUBJECT_name_string,
                                '\n          '.join(scan_list_repr))
        except AttributeError:
            return self.__str__()

    def find_scan_by_protocol(self, protocol_name):
        found_scans = []

        for s in self.scans:
            try:
                if re.match(protocol_name, s.acqp.ACQ_protocol_name):
                    found_scans.append(s)
            except AttributeError as exc:
                logger.warning('Scan in dir {0} has no acqp attribute'.format(
                                s.shortdirname), exc_info=True)
        return(found_scans)

    def scan_finder(self, **kwargs):
        '''
        This is the non-generator version and is provided for convenience.
        See difference between range and xrange
        '''
        return list(self.xscan_finder(**kwargs))

    def xscan_finder(self, **kwargs):
        '''
        Generator of all scans in urrent study that fit the criteria as
        given by kwargs.

        All possible keys are listed and associated to a comparison function
        in the submodule JCAMP_comparison.
        '''
        chosen_comparison = {}
        # find comparison type (regex, array comparison or plain vanilla '==')
        for key in list(kwargs.keys()):
                # remember, the keys here are frozensets
            possible_comp = [k for k in
                             JCAMP_comparison.dictionary
                             if key in k] or [frozenset(['default'])]
            chosen_comparison[key] = JCAMP_comparison.dictionary[
                                                            possible_comp[0]]
        for scn in self.scans:
            if (all(k in scn.acqp.__dict__ and
                    chosen_comparison[k](v, scn.acqp.__dict__[k])
                                         for k, v in list(kwargs.items())) or
                all(k in scn.method.__dict__ and
                    chosen_comparison[k](v, scn.method.__dict__[k])
                                         for k, v in list(kwargs.items()))):
                yield scn

    def _find_adata(self):
        '''
        All keys of adata sets attached to scans in this study
        '''
        adatas = set()
        for scn in self.scans:
            for k in list(scn.adata.keys()):
                adatas.add(k)
        return adatas

    def find_adata_scans(self):
        '''
        All keys *AND SCANS* of adata sets attached to scans in this study
        '''
        
        klist = self._find_adata()
        ad_dict = collections.OrderedDict()
        
        for k in klist: # Populate the ad_dict with the existing keys
            ad_dict[k] = []
            
        for scn in self.scans:
            for ke in list(scn.adata.keys()):
                ad_dict[ke].append(scn.shortdirname)

        return ad_dict
        
    def rm_adata(self, key):
        '''
        Remove adata with given *key* by iterating over the lot
        '''
        for scn in self.scans:
            scn.rm_adata(key)

class StudyCollection(object):
    '''
    A StudyCollection can have multiple studies (in BRUKER speak) which can in
    turn have multiple scans. It can be initialised by pointing it
    to a scan or study. This should be a superclass too (?), e.g. the Patient and
    the Experiment class.
    '''

    def __init__(self):
        '''
        Initialize without loading any data.
        '''
        self.study_instance_uids = []
        self.studies = []

    def __repr__(self):
        '''
        Simple print representation of the Study object

            >>> a=StudyCollection()
            >>> print(a)
            __main__.StudyCollection()
        '''
        return '{0}()'.format(
                    strip_all_but_classname(self, 'StudyCollection'))

    def __str__(self):
        '''
        More elaborate string representation of the StudyCollection object:

            >>> a=StudyCollection()
            >>> a
            __main__.StudyCollection()
        '''
        lines_to_print = [self.__str__()]
        lines_to_print.extend([str(s) for s in self.studies])
        return  '\n'.join(lines_to_print)

    def add_study(self, study=None, by_name=None, by_uid=None):
        '''search through database to find anything that matches root

        :param Study study:
            add the indicated study
        :param string by_name:
            search through database with the indicated string
        :param string by_uid:
            search through database by UID. This is potentially very slow
        '''

        done = False
        if study:
            assert isinstance(study, Study), (
                'Object added to StudyCollection is not a Study')
            done = True

        if by_name:
            if done:
                raise ValueError(
                    'Object added to StudyCollection identified by more than '+
                    'one criterion')
            raise NotImplementedError
            done = True

        if by_uid:
            if done:
                raise ValueError(
                    'Object added to StudyCollection identified by more than '+
                    'one criterion')
            raise NotImplementedError
            done = True

        if done:
            if (study.subject.SUBJECT_patient_instance_uid in
                self.study_instance_uids):
                logger.warning('study previously added')
            else:
                self.studies.append(study)
                try:
                    self.study_instance_uids.append(
                            study.subject.SUBJECT_patient_instance_uid)
                except AttributeError:
                    logger.warning(
                        'SUBJECT_patient_instance_uid not found in study: %s' %
                        study.dirname)
                try:
                    logger.info('study "%s" added to StudyCollection' %
                                study.subject.SUBJECT_study_name)
                except AttributeError:
                    logger.warning(
                        'SUBJECT_study_name not found in study: %s' %
                        study.dirname)
        else:
            logger.warning('no study added to StudyCollection')

    def get_SUBJECT_id(self):
        return [x.subject.SUBJECT_id for x in self.studies]

    def find_scan_by_protocol(self, protocol_name):
        found_scans = []
        for study in self.studies:
            found_scans.extend(study.find_scan_by_protocol(protocol_name))
        return(found_scans)

    def xscan_finder(self, **kwargs):
        '''
        Generator to yield all scans that match te criteria given in kwargs.

        The meat of the code is in Study.xscan_finder()
        '''
        for stdy in self.studies:
            for scn in stdy.xscan_finder(**kwargs):
                yield scn
    def scan_finder(self, **kwargs):
        '''
        This is the non-generator version and is provided for convenience.
        See difference between range and xrange
        '''
        return list(self.xscan_finder(**kwargs))

    def _find_adata(self):
        adatas=set()
        for stdy in self.studies:
            ret = stdy._find_adata()
            adatas = adatas.union(ret)
        return adatas
        
    def find_adata_scans(self, flatten = False):
        
        klist = self._find_adata()
        
        ad_dict = {}
        for k in klist:
            ad_dict[k] = []
            
        for stdy in self.studies:
            c_dict = stdy.find_adata_scans()
            for k in stdy._find_adata():
                ad_dict[k].append(c_dict[k])

        if flatten is False:
            return ad_dict

        else:
            ad_dict_flatten = {}

            for k,v in ad_dict.items():
                ad_dict_flatten[k] = [item for sublist in ad_dict[k] for item in sublist]

            return ad_dict_flatten
       
    def rm_adata(self, key):
        '''
        Remove adata with given *key* by iterating over all studies
        '''
        for stdy in self.studies:
            stdy.rm_adata(key)

class Patient(StudyCollection):
    '''
    A Patient is a special Collection of studies in that the subject_id
    has to agree from one study to the next
    '''
    def __init__(self, patient_name):
        super(Patient, self).__init__()
        self.patient_id = None

        searchdir = os.path.join(dataroot, '*', 'nmr', patient_name) + '*'
        directories = natural_sort(glob.glob(searchdir))
        for dirname in directories:
            study = Study(dirname)
            if not self.patient_id:
                self.patient_id = study.subject.SUBJECT_id
            if self.patient_id != study.subject.SUBJECT_id:
                logger.warning(
                    'loading studies of different SUBJECT_id(s) not permitted')
            else:
                self.add_study(study)

    def __repr__(self):
        '''
        Simple print representation of the Patient object

            >>> a=Patient('readfidTest')
            >>> print(a)
            __main__.Patient("readfidTest")
        '''
        return '{0}("{1}")'.format(
                strip_all_but_classname(self, 'Patient'),
                self.patient_id)

    def __str__(self):
        '''
        More elaborate string representation of the Patient object:

        >>> a=Patient('readfidTest')
        >>> a
        Patient object: Patient("readfidTest")
           studies: Study object: Study("readfidTest.ix1")
           subject: Moosvi, readfidTest
           scans: 1-TriPilot-multi
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  MSME_bas
                  MSME_bas
                  MSME_turbo_bas
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  FLASH_bas(modified)
                  MSME_bas
                  EPI_bas
                  EPI_nav_bas
                  DtiStandard_bas
                  UTE2D
                  UTE3D
                  ZTE
                  FLASH_bas(modified)
        '''
        try:
            study_list_repr = [study.__repr__() for study in self.studies]
            return ('Patient object: Patient("{0}")\n'+
                    '   studies: {1}').format(
                                self.patient_id,
                                '\n          '.join(study_list_repr))
        except AttributeError:
            return self.__str__()

class Experiment(StudyCollection):
    '''
    An Experiment is not BRUKER terminology. It is a grouping of several
    studies (that have scans themselves) by various criteria.
    Typically, the patient name contains a root pointing to a common goal
    of the performed scans. It is a bit like a patient with the relaxation
    of the requirement of having identical SUBJECT_ids.
    '''
    def __init__(self, root=None, absolute_root=False):
        super(Experiment, self).__init__()
        if root:
            self.find_studies(root=root, absolute_root=absolute_root)
            
    @classmethod
    def from_filter(cls, filterstring):
        return cls(filterstring)
    @classmethod
    def from_masterlist(cls, masterlistname):
        bare_experiment = cls(root=None, absolute_root=None)
        bare_experiment.patients = AttrDict()
        bare_experiment.labels = AttrDict()
        if masterlistname is not None:
            conf = configobj.ConfigObj(os.path.join(masterlist_root, masterlistname))
            for k in list(conf.keys()):
                if k != 'General':
                    bare_experiment.patients[k] = collections.OrderedDict()
                    for stdy_str in (s for s in conf[k] if s.startswith('study ')):
                        short_stdy_str = stdy_str.split()[1]
                        long_stdy_str = k+'.'+short_stdy_str
                        study = Study(long_stdy_str)
                        bare_experiment.add_study(study)
                        # populate the dict attributes for ease of access
                        sclbs = conf[k][stdy_str]['scanlabels']
                        for kk in sclbs:
                            sclbs[kk] = os.path.join(long_stdy_str, sclbs[kk]) 
                            if kk not in bare_experiment.labels:
                                bare_experiment.labels[kk]=list()
                            bare_experiment.labels[kk].append(sclbs[kk])
                        bare_experiment.patients[k].update(sclbs)
        bare_experiment.root=masterlistname
        return bare_experiment # not so bare by now since we have added studies
    
    def __repr__(self):
        '''
        Simple print representation of Experiment object

        >>> d=Experiment('NecS3Hs10')
        >>> print(d)
        __main__.Experiment("NecS3Hs10")
        '''
        return '{0}("{1}")'.format(
                strip_all_but_classname(self, 'Experiment'),
                self.root)

    def __str__(self):
        '''
        More elaborate string representation of the Experiment object:

        >>> a=Experiment('NecS3')
        >>> a       # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
        Experiment object: Experiment("NecS3")
           studies: --Total (...)--
                    NecS3Hs01a.iJ1
                    ...
        '''
        try:

            study_list = [study.shortdirname for study in self.studies]

            return ('Experiment object: Experiment("{0}")\n'+
                    '   masterlist: {1}\n' +
                    '   studies: --Total ({2})--\n'+
                    '            {3}').format(self.root,
                                              os.path.basename(self.masterlist_filename),
                                              len(study_list),
                                '\n            '.join(study_list))
        except AttributeError:
            return self.__str__()


    @lazy_property
    def masterlist_filename(self):
        if self.studies:
            new_study_list = []
            for stdy in self.studies:
                if stdy.masterlist_filename is not None:
                    new_study_list.append(stdy)
                else:
                    logger.warning('removing study {0} from Experiment'.format(
                                   stdy.shortdirname))
            self.studies = new_study_list 
        if self.studies:
            default_masterlistname = self.studies[0].masterlist_filename
            for stdy in self.studies:
		
                assert default_masterlistname == stdy.masterlist_filename, \
                    'Studies ({0}, {1}) in Experiment described by multiple masterlists \n({2}, {3})'.format(
                    self.studies[0].shortdirname, stdy.shortdirname,
                    default_masterlistname, stdy.masterlist_filename)
            return default_masterlistname
        else:
            return None
    
    @lazy_property
    def _masterlist(self):
        return configobj.ConfigObj(self.masterlist_filename)

    def find_studies(self, root=None, absolute_root=False):
        if absolute_root:
            searchdir = root
        else:
            searchdir = os.path.join(dataroot, '*', 'nmr', root) + '*'

        self.root = root
        directories = natural_sort(glob.glob(searchdir))
        for dirname in directories:
            study = Study(dirname)
            self.add_study(study)
            
    @lazy_property
    def masterlist_df(self):
        '''iterates over all sections (=Patient labels) in config file
        with the exception of section 'General' and turns it into a
        pandas.DataFrame
        
        This will ignore all rows whose label starts with study. A design assumption is that the
        masterlist only contains a General section and all other sections actually describe
        Patients. The naming of those patients has to follow *exactly* the naming of BRUKER
        patients. Inside those patient sections there may be subsections title "study xxx"
        where xxx is the BRUKER three-letter study abbreviation.'''
        pat_iterator = ((k,v) for k,v in self._masterlist.items() if k!='General')
        df = pandas.DataFrame.from_items(pat_iterator)
        #identify all indices that don't start with "study "
        indexarr = [re.match('study ',x) is None for x in df.index]
        return df[indexarr]
    
    @lazy_property
    def study_masterlist_df(self):
        '''iterates over all sections referring to studies in config file
        and turns it into a pandas.DataFrame
        
        The string used in the index is reconstructed from patient name and
        three-letter study abbreviation. See comments on the assumption of section naming in
        the doc string to masterlist_df.
        '''
        temp_dict = {}
        for k,v in self._masterlist.items():
            for ki,vi in v.items():
                if re.match('study ', ki):                    
                    temp_dict[re.sub('study ', k+'.', ki)]=vi
                    
        df = pandas.DataFrame.from_items(iter(temp_dict.items()))
        return df


if __name__ == '__main__':
    import doctest
    doctest.testmod()
#    import numpy
#    stdy = Study('PhantomOrientation.iY1')
#    print stdy.__repr__()
#    gm=numpy.array([[1,0,0],[0,1,0],[0,0,1]])
#    sr = stdy.xscan_finder(ACQ_grad_matrix = gm)
#    for x in sr:
#        print x
#    sr = stdy.xscan_finder(ACQ_size = [512, 256])
#    for x in sr:
#        print x
#    print '-'*40+'\n initializing ...'
#    NecS3 =Experiment('NecS3')
#    print 'initialized the experiment.done.'
#    import time
#
#    # doing this the firt timeround is slow because all the loading is being
#    # triggered with the JCAMP_file lookups
#    start = time.time()
#    srA = NecS3.find_scan_by_protocol('05')
#    print time.time()-start
#
#    # this run is still slow but that s due to addition loading of the method
#    # files.
#    start = time.time()
#    srB = list(NecS3.xscan_finder(ACQ_protocol_name='05'))
#    print time.time()-start
#    print 'same result from both methods? '+{True:'Yes!',
#                                             False:'No!'}[srA == srB]
#
#    # these ones are positively blazing through once the data is in memory
#    start = time.time()
#    sr = NecS3.find_scan_by_protocol('05')
#    print time.time()-start
#
#    # these ones are positively blazing through once the data is in memory
#    start = time.time()
#    sr = list(NecS3.xscan_finder(ACQ_protocol_name='05'))
#    print time.time()-start
