# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 21:10:05 2013

@author: stefan
"""
import re

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

def git_repo_state():
    import os
    import re
    import subprocess
    import sarpy
    init_path = os.path.abspath(sarpy.__file__)
    worktree = os.path.dirname(os.path.dirname(init_path))
    gitdir = os.path.join(worktree, '.git')

    # this would actually fail if there were not some tags up the tree
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                                ' describe --dirty',
                               stdout=subprocess.PIPE, shell=True)
    (describe, err) = proc.communicate()
    dirty = re.search('-dirty$', describe) is not None 
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                            ' status --porcelain', 
                            stdout=subprocess.PIPE, shell=True)
    (status, err) = proc.communicate()
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                            ' rev-parse HEAD',
                            stdout=subprocess.PIPE, shell=True)
    (sha1, err) = proc.communicate()
    proc = subprocess.Popen('git --git-dir='+gitdir+' --work-tree='+worktree+
                            ' log -1 --format="%ci" '+sha1,
                            stdout=subprocess.PIPE, shell=True)
    (date, err) = proc.communicate()
    if dirty:
        date = date.strip() + ' +++'
    return {'describe': describe.strip(),
            'dirty':dirty,
            'date':date.strip(),
            'status':status.strip(),
            'sha1':sha1.strip()}
    
