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
from sarpy import helpers

def generate(**kwargs):
    '''create histoboards

    conf_file=None: name of config file as kw parameter
    '''

    args = argparse.Namespace()
    # since argsvars will now be pointing at the attribute dictionary of
    # object args you will find below how changes in argsvars modifies the
    # content of args and vice versa. Convenient but maybe a little non-obvious.
    argsvars = vars(args)
    argsvars.update(kwargs)
    args.test = kwargs.get('test', False) # default for dry run
    args.force = kwargs.get('force', False) # default: don't ignore syntax errors

    if args.conf_file:
        print("loading config file %s" % args.conf_file)
        config = ConfigParser.SafeConfigParser()
        base_fname = os.path.join(os.path.expanduser('~'),
                                  'hdata',
                                  args.conf_file)
        if config.read([base_fname]):
            argsvars.update(dict(config.items("Defaults")))
        else:
            raise IOError('Could not read config file %s' % base_fname)

    if args.test:
        print('test mode:\n')
        print(args)
            
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
        for x in [os.path.basename(x) for x in image_files]:
            regmatch = re.match('([0-9]+-)([0-9]+\.{0,1}[0-9]*)([ab]*)(.*)', x)
            if regmatch is not None:
                to_join = ''.join(regmatch.groups()[0:3])
                patient_slice_ID.add(to_join)
            else:
                print 'did not match'+x
    
    patients=set()
    print patient_slice_ID
    violation = False
    for img_ID in patient_slice_ID:
#        img_ID_parsed = re.match('([0-9]+)([a-z])([0-9]+)([ab]*)', img_ID)
        img_ID_parsed = re.match('([0-9]+-)([0-9]+\.{0,1}[0-9]*)([ab]*)(.*)', img_ID)
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
            
    patient_ID_sorted = helpers.natural_sort(patient_slice_ID)
    for patient in patients:
        # find all the images for the *current* patient
        imges = [x for x in patient_ID_sorted 
                 if re.match(patient+'.*', x) is not None]                
        histo_row_files = []
        for k,v in markers.iteritems():
            img_fnames = [glob.glob(os.path.join(root, v, x+'*'))[0] for x in imges]
            print patient, v
            
            histo_row_files.append(os.path.join(args.output, patient+k+'.jpg'))
            
            # check whether file exists and don't bother recreating if it present
            if os.path.exists(histo_row_files[-1]) and not args.overwrite:
                print('previously created file found: %s' % histo_row_files[-1])
            else:
                # annotate file
                labelled_img_fnames=[]
                for image_file in img_fnames:
                    #lbl = os.path.basename(image_file).split()[0]
                    regmatch = re.match('([0-9]+-)([0-9]+\.{0,1}[0-9]*)([ab]*)(.*)',
                                        os.path.basename(image_file))
                    lbl = ''.join(regmatch.groups()[0:3])
                    print image_file
                    #cmd = 'convert -size 100x14 xc:none -gravity center '+ \
                    #      '-stroke black -strokewidth 2 -annotate 0 "%s"' % lbl + \
                    #      '-background none -shadow 100x3+0+0 +repage ' +\
                    #      '-stroke none -fill white -annotate 0 "%s" ' % lbl + \
                    #      "'%s' " % image_file +'+swap -gravity south ' + \
                    #      ' -geometry +0-3 -composite '+\
                    #      "'%s.jpg'" % os.path.join(args.output, lbl+'.png')
                    labelled_outfile = os.path.join(args.output, lbl+'.jpg')
                    cmd = "convert %s -pointsize 155 " % helpers.shellquote(image_file) + \
                          "label:'%s' +swap -gravity Center " % lbl+ \
                          "-append %s" %  helpers.shellquote(labelled_outfile)
                    labelled_img_fnames.append(labelled_outfile)
                    print cmd
                    os.system(cmd)      
                cmd = 'convert +append {0} {1}'.format(
                            ' '.join([helpers.shellquote(x) for x in labelled_img_fnames]),
                            histo_row_files[-1])
                print('running: %s' % cmd)
                if not args.test: os.system(cmd)
                for f in labelled_img_fnames:
                    print('removing %s' % f)
                    os.remove(f)
                    
        largeoutfile = os.path.join(args.output, patient+'.jpg')
        cmd = 'convert -append {0} {1}'.format(
                          ' '.join([helpers.shellquote(x) for x in histo_row_files]), 
                          largeoutfile)
        if os.path.exists(largeoutfile) and not args.overwrite:
            print('not running: %s' % cmd)
        else:
            print('running: %s' % cmd)
            if not args.test: os.system(cmd)
            
        if args.scale:
            outfile = os.path.join(args.output, patient+'-width%s.jpg' % args.factor)
            cmd = 'convert -scale {0} {1} {2}'.format(
                                args.factor, largeoutfile, outfile)
            if os.path.exists(outfile) and not args.overwrite:
                print('not running: %s' % cmd)
            else:
                print('running: %s' % cmd)
                if not args.test: os.system(cmd)
                    
if __name__ == "__main__":
    # we are being run from the commandline
    conf_parser = argparse.ArgumentParser(
        # Turn off help, so we print all options in response to -h
            description=__doc__)
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify full path to config file", metavar="FILE")
    conf_parser.add_argument("-t", "--test",default=False, action="store_true",
                        help='dry run of mriBoard creation')
    conf_parser.add_argument("--force", default=False, action="store_true", 
                        help='ignore violations and attempt to run anyway')     
    conf_parser.add_argument('--output', '-o', 
                        help='output directory')
    conf_parser.add_argument('--overwrite',
                        help='overwrite existing files')                    
    conf_parser.add_argument("--scale", default=False, action="store_true", 
                        help='create thumbnail of factor pixels width')  
    conf_parser.add_argument('--factor','-f', type=int,
                        help='create thumbnail of factor pixels width')
    args = conf_parser.parse_args()
    # at this point args is a 'Namespace' which is defined in the argparse module
    # to get access to its attributes you can simply say, e.g., args.conf_file
    # to treat it like a dictionary you have to use d=vars(args). Then you can
    # read and write to it: d['newattribute']=42 at which point you can verify
    # that args.newatttribute exists and equals 42. The magic of shallow copies...
    generate(**vars(args))
                