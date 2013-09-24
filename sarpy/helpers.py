# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 21:10:05 2013

@author: stefan
"""
import re
import sarpy
import os
import collections
import json
import numpy

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

def generate_summary(masterlist_name):

    root = os.path.join(os.path.expanduser('~/sdata'),
                        masterlist_name,
                        masterlist_name)
    if os.path.exists(os.path.join(root+'_updated.json')):
        fname_to_open = root+'_updated.json'
    else: 
        fname_to_open = root+'.json'
    with open(os.path.join(fname_to_open),'r') as master_file:
        json_str = master_file.read()
        master_list = json.JSONDecoder(
                           object_pairs_hook=collections.OrderedDict
                           ).decode(json_str)
        
    # Generate the filename of the output file based on the masterlist
    filename = root + '_summary.pdf'
    
    setvals = set()
    pats = []
    for k,v in master_list.iteritems():
        pats.append(k)
        for val in v:
            setvals.add(val)
    
    # The order of this setvals is a bit weird, so reverse the set, and turn it into a list:
    headings = natural_sort(list(setvals))
    headings.insert(0,'')
    
    checks =[]
    checks.append(headings)
    for p in pats:
        curr=[]
        curr.append(p)
    
        #Recall the first one is blank to match up with the patients column
        for k in headings[1:]:
            
            if master_list[p][k][0] == '':
                curr.append('No')
            else:
                curr.append('Yes')
        checks.append(curr)
    
        
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter, inch, landscape
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle
     
    doc = SimpleDocTemplate(filename, pagesize=landscape(letter))
    # container for the 'Flowable' objects
    elements = []
     
    checks_size = numpy.array(checks).shape


    colwidths = checks_size[1]*[0.7*inch]
    colwidths.insert(0,0.95*inch)
    
    rowHeights = checks_size[0]*[0.3*inch]
    
    t=Table(checks,colwidths, rowHeights)
    styled = []
    styled = [ ('TEXTCOLOR',(0,0),(0,-1),colors.blue),
               ('TEXTCOLOR',(0,0),(-1,0),colors.blue),
               ('ALIGN',(0,0),(-1,-1),'CENTER'),
               ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
               ('INNERGRID', (1,1), (-1,-1), 0.25, colors.black),
               ('BOX', (1,1), (-1,-1), 0.25, colors.black),
               ('BACKGROUND',(1,1),(-1,-1),colors.lightgreen)
            ]
    
    for cols in xrange(checks_size[0]):
    
        indexes = [i for i,x in enumerate(checks[cols]) if x == 'No']
        
        for idx in indexes:
            styled.append(('BACKGROUND',(idx,cols),(idx,cols),colors.orange))
            
    
    t.setStyle(TableStyle(styled))
    elements.append(t)
    # write the document to disk
    doc.build(elements)             
    
