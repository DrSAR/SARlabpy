# -*- coding: utf-8 -*-
"""
Class definitions for the analysed data structures
"""

import os, errno
import glob
import json
import re
import BRUKER_classes
import nibabel
import numpy
from datetime import datetime

import logging
logger=logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#this is where datagoes that is secondary to acquired data.
# for this to make any sense, there needs to be som mechanism by which the
# source data can be found
adataroot = os.path.expanduser('~/analysed-data')

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

class AData(object):
    def __init__(self, yet_loaded=True, meta=None, **kwargs):
        try:
            self.key=kwargs['key']
        except AttributeError:
            raise AttributeError('To create an adata set, '+
                                 'an identifying key is required')
        if yet_loaded:
            self.data = kwargs['data']
            self.__yet_loaded=True
        else:
            self.__yet_loaded = yet_loaded
        self.meta = meta
        #some data is handy to have close by
        self.parent_uid = meta['parent_uid']
        
    def __getattr__(self, attr):
        if attr == 'data':
            if not self.__yet_loaded:
                self.data = nibabel.load(os.path.join(adataroot, 
                                          self.meta['dirname'],
                                          self.meta['fileroot'])+'.nii')
                self.__yet_loaded = True
                logger.info('Loading of Adata.data now forced for (%s)' %
                        self.meta['dirname'].split(os.path.sep)[1])

            return self.data
        else:
            raise AttributeError('%s is not an attribute' % attr)
        
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
        assert len(datafilename) == 1, (
            'Need precisely one nifti file *nii, found: {0}'.format(
            datafilename))            
        with open(paramfilename[0], 'r') as paramfile:
            meta = json.load(paramfile)
        if not lazy:
            logger.info('loading Nifti file %s' % datafilename[0])
            data = nibabel.load(datafilename[0])
            yet_loaded = True
            return cls(data=data, 
                       meta=meta, 
                       key=meta['key'], 
                       yet_loaded=yet_loaded)
        else:
            yet_loaded = False
            return cls(meta=meta, 
                       key=meta['key'], 
                       yet_loaded=yet_loaded)

            
    @classmethod
    def fromdata(cls, parent=None, data=None, 
                 key=None, meta=None,
                 force=False):
        '''
        :param BRUKER_classes.PDATA_file parent:
            A parent has to be a Scan or some object that provides a UID
        :param array data:
            sourcedata
        :param dict meta:
            meta information in the form of a dictionary. It will have the
            Parent uid added.
        '''
        if parent is None:
            raise ValueError('analysed data has to be based on a parent')
        else:
            if not meta: 
                meta = {}
            # preference is given to key defined in meta parameter
            # but key through the key parameter is admissable as well
            try:
                key_joined = meta['key'] or key
            except KeyError:
                key_joined = key
                
            if key_joined is None:
                raise ValueError('To create an adata set, '+
                                 'an identifying key is required')
            meta['key']=key_joined
            meta['parent_uid'] = parent.uid()
            meta['created_datetime'] = datetime.now().strftime('%c')
            
            # the handle is a fairly safe and human-readable way to identify
            # the originating scan and pdata
            try:
                handle = re.sub(BRUKER_classes.dataroot+os.path.sep, '',
                                parent.filename)
                meta['parent_filename'] = parent.filename
            except AttributeError:
                handle = ''
                meta['parent_filename'] = ''
            handle = re.sub(os.path.sep,'-',handle)
            # folder = meta['dirname'] is where the AData will be stored
            # relative to adataroot
            folder = os.path.join(
                            meta['parent_uid'],
                            '-#-'.join([meta['key'],handle]))    
            meta['dirname'] = folder
            meta['fileroot'] = '-#-'.join([meta['key'], meta['parent_uid']])
            absfolder = os.path.join(adataroot, folder)
            fileroot = os.path.join(absfolder, meta['fileroot'])
            # react if there is a pre-existingnkey of the same name
            if os.path.isdir(absfolder):
                if force:
                   # overwrite
                   os.remove(fileroot+'.json')
                   os.remove(fileroot+'.nii')
                   logger.warning('overwriting analysed data ({0}) for {1}'.
                               format(meta['key'], meta['parent_filename']))
                else:
                   raise IOError(errno.EEXIST, 'Key exists')
            else:
                mkdir_p(absfolder)
                
            nibabel.Nifti1Image(data,numpy.eye(4)).to_filename(fileroot+'.nii')
            with open(fileroot+'.json','w') as paramfile:
                json.dump(meta, paramfile, indent=4)
            logger.info('Saving to {0}.(nii, json)'.format(fileroot))

            
        return cls(data=data, meta=meta, key=meta['key'])

class ADataCollection(object):
    '''
    A dictionary of AData objects. This is  decorator class that is used
    in the BRUKER_classes.Scan object
    '''
    def __init__(self):
        self.__yet_loaded = False
        self.adata = {}
                
    def __set__(self, *args): raise AttributeError('read only!')
        
    def __get__(self, instance, value):
        if not self.__yet_loaded:
            self.__forced_load(instance.pdata_uids) 
        return self.adata
        
    def __forced_load(self, pdata_uids):
        logger.info('delayed loading of ADataCollection now forced ...')
        for pdata_uid in pdata_uids:
            adata_potential = os.path.join(adataroot, pdata_uid)
            if os.path.isdir(adata_potential):
                for adata_sets in os.listdir(adata_potential):
                    adata_candidate = AData.fromfile(
                            os.path.join(adata_potential, adata_sets))
                    self.adata[adata_candidate.key]=adata_candidate
            else:
                logger.info('No analysis data in directory "{0}"'.
                             format(self.dirname))

        self.__yet_loaded = True

if __name__ == '__main__':
    import doctest
    doctest.testmod()