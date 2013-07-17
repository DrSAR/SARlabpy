#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Create tumour boards from histology data.

A config file is used to describe the construction of histo boards. 

Copyright: SARlab members, UBC, Vancouver, 2013
"""
import argparse
import ConfigParser
import os
import glob
import re
from collections import defaultdict
import tempfile

def shellquote(s):
    return "'" + s.replace("'", "'\\''") + "'"
    
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


conf_parser = argparse.ArgumentParser(
    # Turn off help, so we print all options in response to -h
        add_help=False
        )

conf_parser.add_argument("-c", "--conf_file",
                         help="Specify config file", metavar="FILE")
#parse_known_args does not produce an error for unknown args
args, remaining_argv = conf_parser.parse_known_args()
if args.conf_file:
    config = ConfigParser.SafeConfigParser()
    config.read([args.conf_file])
    defaults = dict(config.items("Defaults"))
else: 
    defaults = {}

# Don't surpress add_help here so it will handle -h
parser = argparse.ArgumentParser(
    # Inherit options from config_parser
    parents=[conf_parser],
    # print script description with -h/--help
    description=__doc__,
    # Don't mess with format of description
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )
parser.set_defaults(**defaults)
parser.add_argument("--test","-t",default=False, action="store_true", 
                    help='dry run of tumour-board creation')
parser.add_argument("--force",default=False, action="store_true", 
                    help='ignore violations and attempt to run anyway')                    
parser.add_argument("--scale",default=False, action="store_true", 
                    help='create thumbnail of 5000 pixels width')                    
args = parser.parse_args(remaining_argv)

markers = dict(config.items("Markers"))

if args.test:
    print('test mode:\n')
    print(args)
    print(markers)

# slurp in directory names and check for inconsistencies.

root = os.path.join(os.path.expanduser('~/hdata'), args.experimentname)

patient_slice_ID = set()
for k,v in markers.iteritems():
    print('checking files in "%s"' % v)
    image_files = (glob.glob(os.path.join(root, v, '*jpg')) +
                   glob.glob(os.path.join(root, v, '*jpeg')) + 
                   glob.glob(os.path.join(root, v, '*tif')) + 
                   glob.glob(os.path.join(root, v, '*tiff')) +
                   glob.glob(os.path.join(root, v, '*png')))
    # strip out the patient+slice ID from filename
    patient_slice_ID.update([os.path.basename(x).split(' ')[0] 
                            for x in image_files])

patients=set()
violation = False
for img_ID in patient_slice_ID:
    #print img_ID
    img_ID_parsed = re.match('([0-9]+)([a-z])([0-9]+)([ab]*)', img_ID)
    if img_ID_parsed is None:
        print('convention violation (type I): %s' % img_ID)
        violation = True
    else:
        # since the img_ID is kosher we need to now check whether it exists
        # in all marker sub directories
        for k,v in markers.iteritems():
            files_for_img_ID = glob.glob(os.path.join(root, v, img_ID+'*'))
            if len(files_for_img_ID) == 0:
                print('convention violation (type II): %s not found in %s ' %
                        (img_ID, v))
            elif len(files_for_img_ID) > 1:
                print('convention violation (type III): found {0} in {1}'.format(img_ID, v))
            else: # no violation found -> proceed to assemble names                
                patnr = img_ID_parsed.group(1)
                patients.add(patnr)

if violation and not args.force:
    raise ValueError('violation(s) found\n'+
                     '- please get your filenames in order\n'+
                     '- or consider running with option --force')

# with all violations flagged, we can now assume that an img_ID occurs once
# and only once in each of the marker subfolders.
        
patient_ID_sorted = natural_sort(patient_slice_ID)

for patient in patients:
    # find all the images for the current patient
    imges = [x for x in patient_ID_sorted 
                        if re.match(patient+'[a-z].*', x) is not None]                
    tmpfiles = []
    for k,v in markers.iteritems():
        img_fnames = [glob.glob(os.path.join(root, v, x+'*jpg'))[0] for x in imges]
        print patient, v
        tmpfiles.append(tempfile.mktemp(prefix='histoboardtmp', suffix='.png'))
        cmd = 'convert +append {0} {1}'.format(
                        ' '.join([shellquote(x) for x in img_fnames]),
                        tmpfiles[-1])
        print(cmd)
        if not args.test: os.system(cmd)
        
    cmd = 'convert -append {0} {1}'.format(' '.join(tmpfiles), patient+'.png')
    print(cmd)        
    if not args.test: 
        os.system(cmd)
        if args.scale:
            cmd = 'convert -scale 5000 {0} {1}'.format(patient+'.png', patient+'-scaled.png')
            print(cmd)        
            os.system(cmd)        
        