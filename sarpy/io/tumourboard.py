#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
Create tumour boards from BRUKER MRI data and secondary analysed data.

A config file is used to describe the construction of tumour boards. 
Firstly, it needs to reference the masterlist for the experiment so as to
gather all the relevant scans associated with the appropriate labels.

Copyright: SARlab members, UBC, Vancouver, 2013
"""
import argparse
import ConfigParser
import os
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.plot([1,2,3])
plt.savefig('myfig')
import numpy
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
args = parser.parse_args(remaining_argv)

if args.test:
    print('test mode:\n')
    print(args)

# load the master_list to have easy access to data
with open(args.master_list,'r') as master_file:
    master_list = json.load(master_file)

# determine the layout
rows = [row for row in config.sections() if not((row=='Defaults'))]
n_rows = len(rows)
ref_lbl = config.get(args.ref_row,'label')
ref_pat = master_list.keys()[0]
scn = sarpy.Scan(master_list[ref_pat][ref_lbl])
try:
    ref_data = scn.adata[config.get(row,'adata')]
except ConfigParser.NoOptionError:
    ref_data = scn.pdata[0].data
n_cols=ref_data.shape[2]


# for every patient we will create the same board
for k,v in master_list.iteritems():
    title=k
    # assume an inch x inch for each square with an addition half inch all around
    fig = plt.figure(figsize = (n_cols+1, n_rows+1))
    fig.suptitle(k)
    G = plt.matplotlib.gridspec.GridSpec(n_rows, n_cols)   
    print('\n'+title)
    row_idx = 0
    for row in rows:
        lbl = config.get(row,'label')
        if config.get(row,'type') == 'histo':
            print(lbl,'HISTO!')
        else:
            scn = sarpy.Scan(v[lbl])
            print(lbl, v[lbl])            
            try:
                data = scn.adata[config.get(row,'adata')]
            except ConfigParser.NoOptionError:
                data = scn.pdata[0].data
                #if not 3D, then squash
                if data.ndim>3:
                    print('squashing array')
                    data = data.sum(axis=3)
            print(data.shape)
        for col_idx in xrange(n_cols):
            fig.add_subplot(G[row_idx, col_idx])
            plt.imshow(data[:,:,col_idx])        
            plt.axis('off')
            if row_idx == 0:
                plt.title('Slice {0}'.format(col_idx+1), fontsize=12)
            if col_idx == 0:
                plt.text(0,0,lbl, horizontalalignment='right', 
                         fontsize=10, rotation='vertical')
        row_idx += 1
        
    # Figure spacing adjustments
    fig.subplots_adjust(right = 0.85, wspace = 0.0001, hspace=0.0001)
    G.tight_layout(fig, h_pad = 0.1, w_pad = 0.1)
    #G.update(right = 0.87)
    
    # Saving Figure    
    filename = os.path.join(args.output, k + '.png')                
    plt.savefig(filename, bbox_inches=0, dpi=300)
    plt.close('all')
