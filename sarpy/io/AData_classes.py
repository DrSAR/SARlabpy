# -*- coding: utf-8 -*-
"""
Class definitions for the analysed data structures
"""

import os, errno
import glob
import cPickle
import json
import re
import BRUKER_classes
from lazy_property import lazy_property
from visu_pars_2_Nifti1Header import visu_pars_2_Nifti1Header
import nibabel
import numpy
from datetime import datetime

import logging
logger=logging.getLogger('sarpy.io.Adata_classes')

#this is where datagoes that is secondary to acquired data.
# for this to make any sense, there needs to be som mechanism by which the
# source data can be found
adataroot = os.path.expanduser(os.path.join('~','analysed-data'))
# here is an extra change

def mkdir_p(path):
    '''
    feature that allows the creation of the intermediate directories.
    it mimics the system call of mkdir -p
    '''
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

import collections
class ADataDict(collections.MutableMapping):
    """A dictionary which perfoms some magic on value-by-key assignment

    Pre-existing keys are bein checked against the new key.
    This is inspired by an
    `answer on stackoverflow <http://stackoverflow.com/a/3387975/607562>`_
    to the question of how to perfectly override  a dict.
    """

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs)) # use the free update to set keys

    def __getitem__(self, key):
        return self.store[key]

    def __setitem__(self, key, value):
        '''attribute setter for elements of this souped-up dictionary.
        CHecks for pre-existin keys and only overwrites of forced. This
        has to happen throuh wrapper script. Example:

        >>> import sarpy
        >>> scn = sarpy.Scan('PhantomOrientation.iY1/2')
        >>> scn.store_adata(key='times2',data=scn.pdata[0].data*2, force=True)

        # If the adata exists, do not overwrite!
        >>> scn.store_adata(key='times2',data=scn.pdata[0].data*2) # doctests: ELLIPSIS
        Traceback (most recent call last):
        ...
        AttributeError: AData (key="times2") exists
        To force overwrite, use store_adata with option force=True
        >>> scn.adata['times2']=sarpy.AData.fromdata(parent=scn.pdata[0],
        ...                     data=scn.pdata[0].data*2) # doctests: ELLIPSIS
        Traceback (most recent call last):
        ...
        AttributeError: AData (key="times2") exists
        To force overwrite, use store_adata with option force=True
        '''
        if key in self.store:
            raise AttributeError(('AData (key="%s") exists\nTo force overwrite'+
                                ', use store_adata with option force=True')
                                % key)

        # here we receive the key from the dictionary access
        # hence, we should store it in the meta-data
        value.meta['key'] = key

        # store data on disk

        # the handle is a fairly safe and human-readable way to identify
        # the originating scan and pdata
        handle = re.sub(os.path.sep,'-',
                        re.sub(BRUKER_classes.dataroot+os.path.sep, '',
                               value.meta['parent_filename']))
        folder = os.path.join(
                       value.meta['parent_uid'],
                        '-#-'.join([value.meta['key'],handle]))
        value.meta['dirname'] = folder
        value.meta['fileroot'] = '-#-'.join([value.meta['key'],
                                            value.meta['parent_uid']])
        absfolder = os.path.join(adataroot, folder)
        fileroot = os.path.join(absfolder, value.meta['fileroot'])
        # react if there is a pre-existing key of the same name
        if os.path.isdir(absfolder):
            try:
                os.remove(fileroot+'.json')
                os.remove(fileroot+'.pickle')
            except OSError:
                pass
            else:
                logger.warning('overwriting analysed data ({0}) for {1}'.
                       format(value.meta['key'], value.meta['parent_filename']))
        else:
            mkdir_p(absfolder)

        # note how no directional information is stored here.
        # AData.export2nii is better that way if analyzed-data should be made
        # available externally: in that case, visu_pars geometry information
        # will be used to populate the nifti header.
        with open(fileroot+'.pickle', 'w') as f:
            cPickle.dump(value.data, f)

        with open(fileroot+'.json','w') as paramfile:
            json.dump(value.meta, paramfile, indent=4)
        logger.info('Saving to {0}.(pickle, json)'.format(fileroot))

        logger.info('assigning {0} \nobject {1}'.format(key, value))
        self.store[key] = value


    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __str__(self):
        return str(self.store)

    def __repr__(self):
        return repr(self.store)

def load_AData(pdatas, dirname):
    logger.info('loading of adata (list) now from %s' % dirname)
    special_dict = ADataDict()
    for pdata in pdatas:
        adata_potential = os.path.join(adataroot, pdata.uid())
        if os.path.isdir(adata_potential):
            for adata_sets in os.listdir(adata_potential):
                dirname = os.path.join(adata_potential, adata_sets)
                logger.info('loading adata from %s' % dirname)
                try:
                    adata_candidate = AData.fromfile(dirname, parent=pdata)
                except IOError:
                    logger.warning('Could not load AData')
                else:
                    special_dict.store[adata_candidate.key]=adata_candidate
    return special_dict


class AData(object):
    def __init__(self, **kwargs):
        if 'key' in kwargs:
            self.key=kwargs['key']
        else:
            self.key=None
#            raise AttributeError('To create an adata set, '+
#                                 'an identifying key is required')
        if 'data' in kwargs:
            self.data = kwargs['data']
            logger.info('stored adata in memory: {0}'.format(self.data.shape))

        if 'meta' in kwargs:
            self.meta = kwargs['meta']
        else:
            self.meta = {'key': self.key,
                         'created_datetime': datetime.now().strftime('%c')}
        #some data is handy to have close by
        self.parent_uid = self.meta['parent_uid']
        self.parent = kwargs['parent']

    @lazy_property
    def data(self):
        #find in file when asked to load
        datafilename = os.path.join(adataroot,
                                    self.meta['dirname'],
                                    self.meta['fileroot']+'.pickle')
        logger.info('loading Nifti file %s' % datafilename)
        with open(datafilename) as f:
            data = cPickle.load(f)
        self.__yet_loaded = True
        return data

    def __str__(self):
        '''
        Simple representation (for print() and str() uses) of AData class
        '''
        return 'AData.fromfile("%s")' % self.meta['dirname']
    def __repr__(self):
        '''
        More elaborate representation
        '''
        return ("Analysed data ('{3}') based on PDATA (uid={0})\n"+
                '  created on {1}\n'+
                '  parent: {2}'+
                '').format(self.parent_uid,
                          self.meta['created_datetime'],
                          self.meta['parent_filename'],
                          self.key)

    @classmethod
    def fromfile(cls, filename, parent=None, lazy=True):
        '''
        classmethod that can be used to initiatlize the AData object
        by reading the content from file. To use, issue:

            >>> import glob, os
            >>> adata_dir = glob.glob(adataroot+'/*/*')[0] # take the first dir
            >>> adata_dir_short = re.sub(adataroot+os.path.sep, '', adata_dir)
            >>> AData.fromfile(adata_dir_short)  # doctest:+ELLIPSIS
            Analysed data ('...') based on PDATA (uid=...)
              created on ...
              parent: ...

        Admittedly this will rarely be required and most likely be done from
        some other constructor.

        :param string filename:
            points to the adataroot/uid/uid.*/ directory
        '''
        datadir = os.path.join(adataroot,filename,'*')
        paramfilename = glob.glob(datadir+'json')
        datafilename = glob.glob(datadir+'pickle')


        if len(paramfilename) != 1:
            logger.warn('Need precisely one parameter file *json, found: {0}'.
                            format(paramfilename))
            raise IOError
        if len(datafilename) != 1:
            logger.warn('Need precisely one data file *pickle, found: {0}'.
                            format(datafilename))
        with open(paramfilename[0], 'r') as paramfile:
            meta = json.load(paramfile)

        return cls(meta=meta,
                   key=meta['key'],
                   parent=parent)


    @classmethod
    def fromdata(cls, parent=None, data=None, **kwargs):
        '''
        This creates an AData object from an array (maybe somehing else later)
        It does not store it on disk. This would typically be done by the
        ADataDict on assignment to a key.

        :param BRUKER_classes.PDATA_file parent:
            A parent has to be an object that provides a UID and has a
            visu_pars attribute. Typically a PDATA_file
        :param array data:
            sourcedata
        :param dict(optional) meta:
            meta information in the form of a dictionary. The following keys
            will be added to it: paret_uid, created_datetime, parent_filename.
        '''
        if 'meta' in kwargs:
            meta = kwargs['meta']
        else:
            meta = {}
            meta['parent_uid'] = parent.uid()
            meta['created_datetime'] = datetime.now().strftime('%c')

        if parent is None:
            raise ValueError('analysed data has to be based on a parent')
        else:
            meta['parent_filename'] = parent.filename
            kwargs['parent'] = parent

        if data is None:
            raise ValueError('analysed data has to be passed: \n' +
                             ' ... fromdata(..., data=numpy.array(), ...)')
        else:
            kwargs['data'] = data

        kwargs['meta'] = meta
        return cls(**kwargs)

    def export2nii(self, filename):
        '''
        Export AData content to a named Nifti1 file using the visu_pars-defined
        geometry of the associated parent processed data (PData)

        :param string filename:
            where to write the file, what did you think?

        Example:
            >>> import sarpy
            >>> scn = sarpy.Scan('PhantomOrientation.iY1/2')
            >>> scn.store_adata(key='times2',data=scn.pdata[0].data*2, force=True)
            >>> scn.adata['times2'].export2nii('/tmp/PhantomOrientation-times2.nii.gz')
        '''
        aff, header = visu_pars_2_Nifti1Header(self.parent.visu_pars)
        img_nii = nibabel.nifti1.Nifti1Image(self.data,
                                              aff,
                                              header=header)
        img_nii.to_filename(filename)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
#    scn = BRUKER_classes.Scan('readfidTest.ix1/3')
#    print(scn.adata)
#    scn.adata['new']=AData.fromdata(parent=scn.pdata[0],
#                                data=numpy.empty([10,10,10]))
