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
import pylab
import numpy
import sarpy
import sarpy.fmoosvi.getters
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
with open(os.path.join(os.path.expanduser('~/mdata'),args.master_list),'r') as master_file:
    master_list = json.load(master_file)

# determine the layout
rows = [row for row in config.sections() if not((row=='Defaults'))]
n_rows = len(rows)
ref_lbl = config.get(args.ref_row,'label')

# Start a PDF file of all the animals

#TODO: figure out a way to capture the root of the experiment fast. This willwork
# for XXXSY. But will fail for Y >9
 
rootName = sarpy.Experiment(master_list.keys()[0]).root[0:5]
testPDF = PdfPages(os.path.join(os.path.expanduser('~/mdata'),rootName,rootName+'.pdf'))

# for every patient we will create the same board
for k,v in master_list.iteritems():
    ref_scn = sarpy.Scan(v[ref_lbl][0])
    try:
        ref_data = ref_scn.adata[config.get(defaults['ref_row'],'adata')]
    except ConfigParser.NoOptionError:
        ref_data = ref_scn.pdata[0]
    ref_filename = tempfile.mktemp(suffix='.nii')
    ref_data.export2nii(ref_filename)
    n_cols=ref_data.data.shape[2]

    
    title = k
    fig = pylab.figure()
    fig.suptitle(k)
    G = pylab.matplotlib.gridspec.GridSpec(n_rows, n_cols, wspace=0.0, hspace=0.0)   
    print('\n'+'-'*80+'\n'+title)

    row_idx = 0
    for row in rows:
        row_conf = dict(config.items(row))
        lbl =row_conf.pop('label')
        subtitle = row_conf.pop('subtitle','')
        fname = v[lbl][0]
        bbox = numpy.array([float(x) for x in v[lbl][1]])
        fig.add_subplot(G[row_idx, 0])
        pylab.axis('off')
        pylab.text(-0.35,0.5,'\n'.join([lbl,subtitle,fname]), 
                 horizontalalignment='center', 
                 fontsize=3, rotation='vertical')

        if fname == '':
            pylab.text(0.5,0.5,'Data not available',
                       horizontalalignment='center',
                       fontsize=8)

            row_idx += 1
            continue

        # allowed types can be:
        # img, plot, vtc, histo

        if row_conf.get('type', None) == 'img':
            
            #TODO: this statement is hre because of **row_conf in imhow
            row_conf.pop('type', None)

            scn = sarpy.Scan(fname)
            print(lbl, fname,scn.acqp.ACQ_protocol_name)            
            adata_key = row_conf.pop('adata', None)
            if adata_key is not None:
                data = scn.adata[adata_key]
                # Set the image limits for adata
                
                cl = sarpy.fmoosvi.getters.get_image_clims(data.data)

            else:
                data = scn.pdata[0]
                cl = None
                
            xdata = data.data
            
            
            
            resample_flag = row_conf.pop('resample', False) 
            if resample_flag and config.getboolean(row,'resample'):
                raise NotImplementedError('please fix resampling')
            #print('resampling {0}\n{1}\n onto {2}'.format(scn,data,ref_data))
            #src_filename = tempfile.mktemp(suffix='.nii')
            #data.export2nii(src_filename)            
            #xdata, xdata_sitk_image = sir.resample_onto(src_filename, ref_filename)
            #xdata = sir.resample_onto_pdata(data, ref_data)
            #os.remove(src_filename)
            #find out where the 0 1 etc end up.
#            print(numpy.mean(numpy.mean(xdata, axis=0),axis=0))

            for col_idx in xrange(min(n_cols, xdata.shape[2])):
                fig.add_subplot(G[row_idx, col_idx])
                bbox_pxl = (bbox.reshape(2,2).T*xdata.shape[0:2]).T.flatten()
                t=pylab.imshow(xdata[bbox_pxl[0]:bbox_pxl[1],
                                 bbox_pxl[2]:bbox_pxl[3],col_idx],
                           **row_conf)
                           
                if cl is not None:
                    t.set_clim(cl)
                pylab.axis('off')
                
                if row_idx == 0:
                    pylab.title('Slice {0}'.format(col_idx+1), fontsize=8)
                    
            row_idx += 1
        
        elif row_conf.get('type', None) == 'vtc':
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            raise NotImplementedError('do not know how to draw VTCs')
        elif row_conf.get('type', None) == 'plot':

            scn = sarpy.Scan(fname)
            print(lbl, fname,scn.acqp.ACQ_protocol_name)            
            adata_key = row_conf.pop('adata', None)
            if adata_key is not None:
                data = scn.adata[adata_key]
            else:
                data = scn.pdata[0]
            xdata = data.data

            for col_idx in xrange(min(n_cols, xdata.shape[0])):
                fig.add_subplot(G[row_idx, col_idx])
                pylab.plot(xdata[col_idx,0,:],xdata[col_idx,1,:])  

                pylab.xlim(0,numpy.nanmax(xdata[:,0,:]))
                pylab.ylim(0,numpy.nanmax(xdata[:,1,:]*1.3))
                pylab.locator_params(axis='both',which='both',
                                     nbins=3) 
                pylab.xlabel('Time (ms)',fontsize=5)
                pylab.ylabel('Enhancement (au)',fontsize=3)
                
                ax = pylab.gca()
                
                if row_idx == 0:
                    pylab.title('Slice {0}'.format(col_idx+1), fontsize=8)
                if col_idx != 0:
                    ax.axes.get_yaxis().set_visible([])
                    ax.axes.get_xaxis().set_visible([])
                
#TODO: finda way to put the axis labels on th INSIDE of the plot
                
                for label in ax.get_xticklabels():
                    label.set_fontsize(2)

                for label in ax.get_yticklabels():
                    label.set_fontsize(2)

            row_idx += 1
            
            
        elif row_conf.get('type', None) == 'histo':
            raise NotImplementedError('do not know how to draw histos')

        
    testPDF.savefig(fig)      
    
    # Saving Figure    
    outputfilename = os.path.expanduser(args.output)
    filename = os.path.join(outputfilename, k + '.png')                
    pylab.savefig(filename, bbox_inches=0, dpi=500)
    pylab.close('all')

#os.remove(ref_filename)

testPDF.close()
