# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 21:10:05 2013

@author: stefan
"""
import re
import os
import subprocess
import collections
from collections import defaultdict
import json
import numpy
import pandas
import sarpy


def flatten_list(l):

    if len(l) < 2:
        print('List is already flat')
        return l

    else:
        # From: http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
        return [item for sublist in l for item in sublist]

def natural_sort(l):
    '''
    Sort a list by a natural sort order (number ascending) and even if
    bracketed by blocks of alpha-characters.

    This is based on code a blog post by `Jeff Atwood <http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html>`_.
    It is also discussed on `stackoverflow <http://stackoverflow.com/questions/4836710/does-python-have-a-built-in-function-for-string-natural-sort>`_.
    '''
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

def shellquote(s):
    return "'" + s.replace("'", "'\\''") + "'"

def git_repo_state():
    init_path = os.path.abspath(__file__)
    worktree = os.path.dirname(os.path.dirname(init_path))
    gitdir = os.path.join(worktree, '.git')

    # this would actually fail if there were not some tags up the tree
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                                ' describe --dirty',
                               stdout=subprocess.PIPE, shell=True)
    (describe, err) = proc.communicate()
    dirty = re.search('-dirty$', describe.decode("utf-8")) is not None 
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                            ' status --porcelain', 
                            stdout=subprocess.PIPE, shell=True)
    (status, err) = proc.communicate()
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                            ' rev-parse HEAD',
                            stdout=subprocess.PIPE, shell=True)
    (sha1, err) = proc.communicate()
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                            ' log -1 --format="%ci" '+sha1.decode('utf-8'),
                            stdout=subprocess.PIPE, shell=True)
    (date, err) = proc.communicate()
    date = date.decode('utf-8').strip()
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                            ' rev-parse --abbrev-ref HEAD',
                            stdout=subprocess.PIPE, shell=True)
    (branch, err) = proc.communicate()

    if dirty:
        date = date + ' +++'
    return {'describe': describe.decode('utf-8').strip(),
            'dirty':dirty,
            'date':date,
            'status':status.decode('utf-8').strip(),
            'sha1':sha1.decode('utf-8').strip(),
            'branch':branch.decode('utf-8').strip()}

def replace_nan_with(arr, value=0):
    '''replace the NaN values (isnan) with value (default=0) for better ImageJ
    compatibility'''
    arr[numpy.where(numpy.isnan(arr))]=value
    return arr

def smooth_SG(y, window_size, order, deriv=0, rate=1):

    # Implementation of a Savitzky-Golay filter
    # Acquired from: 
    # http://stackoverflow.com/questions/22988882/how-to-smooth-a-curve-in-python

    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int (%s)" % msg)
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = list(range(order+1))
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def adata2df(scn_name,adatadict,infodict,roi_adata=None):

    # Function to turn adata into dataframes
    # See http://localhost:8889/user/fmoosvi/notebooks/OxygenMRI/q-dOE-MRI/adata2dataframe.ipynb

    '''s='OEP8Cs12.Q21/7'
        Usage example:
        sarpy.helpers.adata2df(scn_name=s,adatadict={'dOE1':'OEdraft2',
                                     'dOE2':'OEdraft2',
                                     'dOE3':'OEdraft2'}, infodict={'AnimalID':s.split('.')[0],
                                                                   'Treatment':'Ctrl',
                                                                   'Loaction':'sc'}) '''

    def array2df(arr, columname=None):
        idx = numpy.where(numpy.isfinite(arr))
        # if arr is a 3D array then we can recreate the tuples by 
        # using zip(idx[0],idx[1],idx[2]) but that breakes down for
        # arrays that are less than 3D
        df = pandas.DataFrame(arr[idx], index=[tuple(x) for x in numpy.array(idx).transpose()],
                             columns=[columname])
        # the 'name' of the returned column will be '0' - not too descriptive, I know.
        return df
    
    scn = sarpy.Scan(scn_name)
    if roi_adata: #in case we need to mask out a region (e.g muscle vs tumour)
        roi = scn.adata[roi_adata].data
    else:
        roi=1
   
    assert type(adatadict) is dict
    assert type(infodict) is dict
    
    # append all the adatas into a list
    adata_arrayList = []
    df_list = []

    for colName,adataName in adatadict.items():
        adata_arrayList.append(roi*scn.adata[adataName].data)
        df_list.append(array2df(adata_arrayList[-1],columname=colName))

    # now iterate through the dfs and then merge them
    if len(df_list)==1:
        mergeddf = df_list[0]
    else:
        mergeddf = df_list[0]
        for idf in range(1,len(df_list)):
            mergeddf=pandas.merge(mergeddf,df_list[idf],left_index=True, right_index=True,how='outer')     

    # Now add all the other information like treatments, animal name etc...        
    finaldf = pandas.concat([mergeddf,pandas.DataFrame(infodict,index=mergeddf.index)], axis=1) 
    
    return finaldf