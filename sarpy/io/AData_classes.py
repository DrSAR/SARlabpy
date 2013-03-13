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
    def __init__(self, data=None, meta=None):
        self.data = data
        self.meta = meta
        #some data is handy to have close by
        self.parent_uid = meta['parent_uid']
        
    def __str__(self):
        '''
        Simple representation (for print() and str() uses) of AData class
        '''
        return 'AData.fromfile("%s")' % self.meta['filename']
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
    def fromfile(cls, filename):
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
        data = nibabel.load(datafilename[0])
        return cls(data=data, meta=meta)

            
    @classmethod
    def fromdata(cls, parent=None, data=None, meta=None):
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
            meta['parent_uid'] = parent.uid()
            meta['created_datetime'] = datetime.now().strftime('%c')
            
            # try to find some recognizable information on parent
            try:
                handle = re.sub(BRUKER_classes.dataroot, '', parent.filename)
                meta['parent_filename'] = parent.filename
            except AttributeError:
                handle = ''
                meta['parent_filename'] = ''
            handle = re.sub(os.path.sep,'-',handle)
            
            #create new uid based on already existing uids
            try:
                existing_adatas = os.listdir(os.path.join(adataroot, 
                                                          meta['parent_uid']))
            except OSError:
                existing_adatas = ''
            used_ids = []
            for folder in existing_adatas:
                used_ids.append(int(re.sub('-.*','',folder)))
            if used_ids:
                new_id = str(max(used_ids) + 1)
            else:
                new_id = '1'
            folder = os.path.join(
                            meta['parent_uid'],
                            '-'.join([new_id, handle]))    
            meta['filename'] = folder
            target = os.path.join(adataroot, folder)
            assert not(os.path.exists(target)), ('target directory '+
                    '(%s) exists - this should not have happened' % target)
            mkdir_p(target)
            fileroot = os.path.join(target, 
                                    '.'.join([meta['parent_uid'], new_id]))
            nibabel.Nifti1Image(data,numpy.eye(4)).to_filename(fileroot+'.nii')
            with open(fileroot+'.json','w') as paramfile:
                json.dump(meta, paramfile, indent=4)

            
        return cls(data=data, meta=meta)

if __name__ == '__main__':
    import doctest
    doctest.testmod()