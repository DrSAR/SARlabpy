#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Create a master list for an "EXPERIMENT", i.e. a collection of studies

A config file is used to describe the construction of this master list.

Copyright: SARlab members, UBC, Vancouver, 2013
"""
import argparse
import collections
import ConfigParser
import re
import json
import sarpy

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
                    help='dry run of tumour-board creation')
args = parser.parse_args(remaining_argv)

if defaults == {}:
    raise IOError('no config file read')

labels = [lbl for lbl in config.sections() if not((lbl=='MasterList') or
                                            re.match('EXCEPTION.',lbl))]
exc_labels = [lbl for lbl in config.sections() if (not(lbl=='MasterList') and
                                            re.match('EXCEPTION.',lbl))] 
try:
    patient_exclude = [x.strip() for x in 
                   config.get('MasterList','patient_exclude').split(',')]
except ConfigParser.NoOptionError:
    patient_exclude = []
    
# get unique name of patients
expt = sarpy.Experiment(defaults['experimentname'])
patname_list=sorted(list(set(expt.get_SUBJECT_id())))

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
        study_nr = int(config.get(lbl,'study'))
        prot_name = config.get(lbl,'protocol_name')
        master_sheet[patname][lbl] = ["",[]]          
        try:
            bbox =  [float(x.strip()) for x in 
                 config.get(lbl,'bounding_box').split()]
        except ConfigParser.NoOptionError:
            bbox=[0,1,0,1]

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
                
# now, take care of the exceptions.
# as an example, this is what the cponfig file looks like
#    [EXCEPTION.24h-LL.HerP2Bs03]
#    protocol_name=04_ubcLL2
#    study=3 # instead of 1
#    scan=0
for exc_lbl in exc_labels:
    (exc, lbl, patname) = exc_lbl.split('.')
    study_nr = int(config.get(exc_lbl,'study'))
    prot_name = config.get(exc_lbl,'protocol_name')
    pat = sarpy.Patient(patname)
    scns = pat.studies[study_nr].find_scan_by_protocol(prot_name)
    try:
        scn = scns[int(config.get(exc_lbl,'scan'))]
    except IndexError:
        master_sheet[patname][lbl] = ""
        print('WARNING: exception for %s could not be satisfied ' % exc_lbl)
    else:
        master_sheet[patname][lbl] = scn.shortdirname                

print('ignoring {0}'.format(' '.join(patient_exclude)))
if args.test:
    print('test mode:\n')
    print('%s\n' % args)
    print(json.dumps(master_sheet))#, indent=4))
else:
    with open(args.output,'w') as outfile:
        json_str = json.dumps(master_sheet, outfile, indent=4)
        y = re.sub(r'\s\s+(\d+)', lambda match: r' {}'.format(
                    match.group(1), json_str), json_str)
        outfile.write(y)
        print('wrote to %s' % args.output)
