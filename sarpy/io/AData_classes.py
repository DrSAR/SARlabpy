# -*- coding: utf-8 -*-
"""
Class definitions for the analysed data structures
"""

import os, errno
import glob
import json
import re
import BRUKER_classes
from lazy_property import lazy_property
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
        if key in self.store:
            logger.warn('ADate set is being overwriten for key %s' % key)
        
        # here we receve the key from the dictionary access
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
        # react if there is a pre-existingnkey of the same name
        if os.path.isdir(absfolder):
           os.remove(fileroot+'.json')
           os.remove(fileroot+'.nii')
           logger.warning('overwriting analysed data ({0}) for {1}'.
                       format(value.meta['key'], value.meta['parent_filename']))
        else:
            mkdir_p(absfolder)

        nibabel.Nifti1Image(value.data,numpy.eye(4)).to_filename(
                                                            fileroot+'.nii')
        with open(fileroot+'.json','w') as paramfile:
            json.dump(value.meta, paramfile, indent=4)
        logger.info('Saving to {0}.(nii, json)'.format(fileroot))
 
        print('assigning {0} \nobject {1}'.format(key, 
                          value))
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
                adata_candidate = AData.fromfile(dirname)
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
            
        if 'meta' in kwargs:
            self.meta = kwargs['meta']
        else:
            self.meta = {'key': self.key, 
                         'created_datetime': datetime.now().strftime('%c')}
        #some data is handy to have close by
        self.parent_uid = self.meta['parent_uid']
        
    @lazy_property
    def data(self):
        #find in file when asked to load
        datafilename = os.path.join(adataroot,
                                    self.meta['dirname'],
                                    self.meta['fileroot']+'.nii')
        logger.info('loading Nifti file %s' % datafilename)
        data = nibabel.load(datafilename)
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
        return ('Analysed data based on PDATA (uid={0})\n'+
                '  created on {1}\n'+
                '  parent: {2}'+
                '').format(self.parent_uid,
                          self.meta['created_datetime'],
                          self.meta['parent_filename'])

    @classmethod
    def fromfile(cls, filename, lazy=True):
        '''
        classmethod that can be used to initiatlize the AData object 
        by reading the content from file. To use, issue:
            
            >>> AData.fromfile('2.16.756.5.5.100.1384712661.15242.1362627392.1/12--NecS3Hs10.iK1-5-pdata-1')
        
        Admittedly this will rarely be required and most likely be done from
        some other constructor.
        
        :param string filename:
            points to the adataroot/uid/uid.*/ directory
        '''
        datadir = os.path.join(adataroot,filename,'*')
        paramfilename = glob.glob(datadir+'json')
        datafilename = glob.glob(datadir+'nii')
        
        
        assert len(paramfilename) == 1, (
            'Need precisely one parameter file *json, found: {0}'.format(
            paramfilename))
        if len(datafilename) != 1:
            logger.warn('Need precisely one nifti file *nii, found: {0}'.
                            format(datafilename))            
        with open(paramfilename[0], 'r') as paramfile:
            meta = json.load(paramfile)

        return cls(meta=meta, 
                   key=meta['key'])

            
    @classmethod
    def fromdata(cls, parent=None, **kwargs):
        '''
        This creates an AData object from an array (maybe somehing else later)
        It does not store it on disk. This would typically be done by the
        ADataDict on assignment to a key.
        
        :param BRUKER_classes.PDATA_file parent:
            A parent has to be a Scan or some object that provides a UID
        :param array data:
            sourcedata
        :param dict meta:
            meta information in the form of a dictionary. It will have the
            Parent uid added.
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
            
        kwargs['meta'] = meta
        return cls(**kwargs)

if __name__ == '__main__':
#    import doctest
#    doctest.testmod()
    scn = BRUKER_classes.Scan('readfidTest.ix1/3')
    print(scn.adata)
    scn.adata['new']=AData.fromdata(parent=scn.pdata[0], 
                                data=numpy.empty([10,10,10]))