#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Create a master list for an "EXPERIMENT", i.e. a collection of studies

A config file is used to describe the construction of this master list.

Copyright: SARlab members, UBC, Vancouver, 2013
"""
import os
import argparse
import collections
import ConfigParser
import re
import json
import pprint
import sarpy
import sarpy.fmoosvi.getters

def generate(**kwargs):
    '''Generate a masterlist from a config file.

    config_file=None: name of config file as kw parameter
    test=False: Only  perform dry-run '''
    
    args = argparse.Namespace()
    argsvars = vars(args)
    argsvars.update(kwargs)
    args.test = kwargs.get('test', False) # default for dry run
    
    if args.conf_file:
        config = ConfigParser.SafeConfigParser()
        base_fname = os.path.join(os.path.expanduser('~'),  
                                  'sdata',
                                  args.conf_file)
        config.read([base_fname])
        argsvars.update(dict(config.items("MasterList")))

    
    labels = [lbl for lbl in config.sections() if not((lbl=='MasterList') or
                                                re.match('EXCEPTION.',lbl))]
    exc_labels=collections.defaultdict(list)
    for lbl in  [lbl for lbl in config.sections() if re.match('EXCEPTION.',lbl)]:
        (exc, lbl, patname) = lbl.split('.')
        exc_labels[lbl].append(patname)
        
    try:
        patient_exclude = [x.strip() for x in 
                       config.get('MasterList','patient_exclude').split(',')]
    except ConfigParser.NoOptionError:
        patient_exclude = []
        
    # get unique name of patients
    expt = sarpy.Experiment(args.experimentname)
    patname_list=sarpy.natural_sort(list(set(expt.get_SUBJECT_id())))
    
    master_sheet = collections.OrderedDict()
    
    # first do the regular assignments
    # as an example, this is what the config file states:
    #    [anatomy]
    #    protocol_name=05_IR-RARE-anatomy_highres
    #    study=1
    #    scan=0
    
    for patname in (x for x in patname_list if x not in patient_exclude):
        print patname
        pat = sarpy.Patient(patname)
        master_sheet[patname]={}
    
        for lbl in labels:
            configuration = dict(config.items(lbl))
            master_sheet[patname][lbl] = ["",[]]          
    
    # now, see whether there are exceptions defined or this label
    # as an example, this is what the config file looks like
    #    [EXCEPTION.24h-LL.HerP2Bs03]
    #    protocol_name=04_ubcLL2
    #    study=3 # instead of 1
    #    scan=0
    
            if patname in exc_labels[lbl]:
                configuration = dict(config.items(
                                '.'.join(['EXCEPTION',lbl,patname])))
    
            bbox=configuration.get('bounding_box', '0 1 0 1').split()
            prot_name = configuration['protocol_name']
            study_nr = int(configuration['study'])
    
            try:
                scns = pat.studies[study_nr].find_scan_by_protocol(prot_name)
            except IndexError:
                print('warning (STUDY #%i not found) "%s" for %s' 
                    % (study_nr, lbl, patname))                
            else:
                scn_nr = int(config.get(lbl,'scan'))
                try:
                    scn = scns[scn_nr]
                except IndexError:
                    print('warning (SCAN #%i not found amongst the %s) "%s" for %s' 
                        % (scn_nr, prot_name, lbl, patname))                
                else:
                    master_sheet[patname][lbl] = [scn.shortdirname, bbox]
    
    print('ignoring {0}'.format(' '.join(patient_exclude)))
    
    if args.test:
        print('test mode:\n')
        print('%s\n' % args)
        pprint.pprint(json.loads(json.dumps(master_sheet)))
    else:
        outputfilename = os.path.join(os.path.expanduser('~/sdata'),args.output)
        with open(outputfilename,'w') as outfile:
            json_str = json.dumps(master_sheet, indent=4)
            #what is going on?
            y = re.sub(r'\s\s+(\d+)', lambda match: r' {}'.format(
                        match.group(1), json_str), json_str)
            outfile.write(y)
            print('wrote to %s' % outputfilename)
        
        ###### Reload masterlist file as a json file after checking if roi exists 
        ###### if not, skip this part        
        check_list = expt.find_adata()        
        if 'roi' in check_list:
            with open(outputfilename,'r') as master_file:
                json_str = master_file.read()
                master_list = json.JSONDecoder(
                                   object_pairs_hook=collections.OrderedDict
                                       ).decode(json_str)                   
            ###### Updating bboxes #####
            for patname,v in master_list.iteritems():  
                # First get all the labels and put it in a list
                lbl_list = master_list[patname].keys()        
                # Search for all the labels that have roi and create an roi list 
                roi_labels = [r for r in lbl_list if 'roi' in r]        
                # Iterate over the roi list, get the updated bbox, check for same day-ness
                for r_lbl in roi_labels:        
                    scn_name = master_list[patname][r_lbl][0]
                    new_bbox = sarpy.fmoosvi.getters.get_roi_bbox(scn_name,'roi',type='pct')
                    search_string = r_lbl.split('-',1)[-1]        
                    # Iterate over the label list, transfer the new bbox into the master list
                    for lbl in lbl_list:        
                        if search_string in lbl and len(master_list[patname][lbl][1])==4: 
                            # If 0h/24h/48h/ exists and if scan + bbox exists, update bbox
                                master_list[patname][lbl][1] = new_bbox                
            #Save updated file to disk                        
            outputfilename = os.path.join(os.path.expanduser('~/sdata'),args.output)        
            json.dump(master_list, open(outputfilename,'w'),indent=4)        
            print('wrote updated bbox file to %s' % outputfilename)        
        else:        
            print("No rois were found, bboxes are likely incorrect. Proceed with caution!")

if __name__ == "__main__":
    # we are being run from the commandline
    conf_parser = argparse.ArgumentParser(
        # Turn off help, so we print all options in response to -h
            add_help=False)
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify full path to config file", metavar="FILE")
    #parse_known_args does not produce an error for unknown args
    args, remaining_argv = conf_parser.parse_known_args()
    if args.conf_file:
        config = ConfigParser.SafeConfigParser()
        base_fname = os.path.join(os.path.expanduser('~'),  
                                  'sdata',
                                  args.conf_file)
        config.read([base_fname])
        defaults = dict(config.items("MasterList"))
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
                        help='dry run of masterlist creation')
    args = parser.parse_args(remaining_argv)
    
    # at this point args is a 'Namespace' which is defined in the argparse module
    # to get access to its attributes you can simply say, e.g., args.conf_file
    # to treat it like a dictionary you have to use d=vars(args). Then you can
    # read and write to it: d['newattribute']='42 at which point you can verify
    # that args.newatttribute exists and equals 42. The magic of shallow copies...
    if defaults == {}:
        raise IOError('no config file read')
    else:
        generate(**vars(args))