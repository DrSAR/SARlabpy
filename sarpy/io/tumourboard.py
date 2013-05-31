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
matplotlib.use('Agg')# where did I come from !?
import matplotlib.pyplot as plt
import numpy
import sarpy
import sarpy.ImageProcessing.resample_onto as sir
import tempfile
from matplotlib.backends.backend_pdf import PdfPages


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
ref_scn = sarpy.Scan(master_list[ref_pat][ref_lbl][0])
try:
    ref_data = ref_scn.adata[config.get(defaults['ref_row'],'adata')]
except ConfigParser.NoOptionError:
    ref_data = ref_scn.pdata[0]

ref_filename = tempfile.mktemp(suffix='png')
ref_data.export2nii(ref_filename)

n_cols=ref_data.data.shape[2]
for i in xrange(n_cols):
    ref_data.data[:,:,i]=i

testPDF = PdfPages('testPDF.pdf')


# for every patient we will create the same board
for k,v in master_list.iteritems():

    # Start a PDF file of all the animals
    
    title = k
    # assume an inch x inch for each square with an addition half inch all around
    fig = plt.figure(figsize = (n_cols+1, n_rows+1))
    fig.suptitle(k)
    G = plt.matplotlib.gridspec.GridSpec(n_rows, n_cols)   
    print('\n'+title)
    row_idx = 0
    for row in rows:
        lbl = config.get(row,'label')
        fname = v[lbl][0]
        bbox = numpy.array([float(x) for x in v[lbl][1]])
        fig.add_subplot(G[row_idx, 0])
        plt.axis('off')
        plt.text(-0.1,0.6,lbl, horizontalalignment='center', 
                     fontsize=5, rotation='vertical')

        if config.get(row,'type') == 'histo':
            print(lbl,'HISTO!')
        elif fname == "":
            print('no scan in masterlist')
            continue
        else:
            scn = sarpy.Scan(fname)
            print(lbl, fname,scn.acqp.ACQ_protocol_name)            
            try:
                data = scn.adata[config.get(row,'adata')]
            except ConfigParser.NoOptionError:
                data = scn.pdata[0]
            xdata = data.data
        try:
            resample = config.getboolean(row,'resample')
        except ConfigParser.NoOptionError:
            resample = False
        if resample:
            print('resampling {0}\n{1}\n onto {2}'.format(scn,data, ref_data))
            src_filename = tempfile.mktemp(suffix='png')
            data.export2nii(src_filename)            
            xdata, xdata_sitk_image = sir.resample_onto(src_filename, ref_filename)
            #find out where the 0 1 etc end up.
#            print(numpy.mean(numpy.mean(xdata, axis=0),axis=0))

        for col_idx in xrange(min(n_cols, xdata.shape[2])):
            fig.add_subplot(G[row_idx, col_idx])
            bbox_pxl = (bbox.reshape(2,2).T*xdata.shape[0:2]).T.flatten()
            plt.imshow(xdata[bbox_pxl[0]:bbox_pxl[1],
                             bbox_pxl[2]:bbox_pxl[3],col_idx])        
            plt.axis('off')
            plt.title('Slice {0}'.format(col_idx+1), fontsize=8)
        row_idx += 1
        
    # Figure spacing adjustments
    #fig.subplots_adjust(right = 0.85, wspace = 0.0001, hspace=0.0001)
    #G.tight_layout(fig, h_pad = 0.01, w_pad = 0.01)
    #G.update(right = 0.87)
        
    testPDF.savefig(fig)      
    
    # Saving Figure    
    filename = os.path.join(args.output, k + '.png')                
    plt.savefig(filename, bbox_inches=0, dpi=500)
    plt.close('all')

testPDF.close()