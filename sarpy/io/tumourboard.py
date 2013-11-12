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
import collections
import json
import matplotlib
matplotlib.use('Agg')# where did I come from !?
import pylab
import numpy
import sarpy
import sarpy.fmoosvi.getters
import tempfile
import scipy
import copy
import re
from matplotlib.backends.backend_pdf import PdfPages


def determine_figure_size(n_rows,n_cols):
    
    aspect = numpy.true_divide(n_rows,n_cols)
       
    if n_rows <= 4 and n_cols <= 6:
        #moderately tested
        figure_size = (6,5*aspect)
        font_modifier = 0
        
    elif (n_rows >= 4 and n_rows <=8) and (n_cols>=6 and n_cols <= 10) :
        #moderately tested               

        figure_size = (8,8*aspect)
        font_modifier = 3
    
    elif n_rows <=4 and n_cols >=10:
        #moderately tested               
        figure_size = (10,14*aspect)
        font_modifier = 1

    elif n_rows >=8 and n_cols <=10:
        #moderately tested        
        figure_size = (16,8*aspect)
        font_modifier = 6      

    elif n_rows >= 6 and n_cols > 10:

        figure_size = (18,19*aspect)
        font_modifier = 6
        
    else:
        print('row_size = {0} and col_size = {1}'.format(n_rows,n_cols))
        raise NotImplementedError('Please code in this situation, missed it, figure size')
        
    return figure_size,font_modifier

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
    conf_file_name = [os.path.join(os.path.expanduser('~/sdata'),args.conf_file)]
    config.read(conf_file_name)
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
with open(os.path.join(os.path.expanduser('~/sdata'),args.master_list),'r') as master_file:
    json_str = master_file.read()
    master_list = json.JSONDecoder(
                       object_pairs_hook=collections.OrderedDict
                       ).decode(json_str)    

# determine the layout
rows = [row for row in config.sections() if not((row=='Defaults'))]
n_rows = len(rows)
ref_lbl = config.get(args.ref_row,'label')

# Start a PDF file of all the animals

#TODO: figure out a way to capture the root of the experiment fast. This willwork
# for XXXSY. But will fail for Y >9

rootName = str(conf_file_name[0]).split('/')
pdfName = os.path.splitext(rootName[-1])[0]

testPDF = PdfPages(os.path.expanduser(os.path.join('~/sdata',rootName[-2],args.output,pdfName+'.pdf')))

# for every patient we will create the same board
for k,v in master_list.iteritems():

    try:
        ref_scn = sarpy.Scan(v[ref_lbl][0])
    except(IOError,KeyError):
        print('Ref Scan failed for {0},{1} \n \n'.format(k,ref_lbl))
        continue
    try:
        ref_data = ref_scn.adata[config.get(defaults['ref_row'],'adata')]
    except ConfigParser.NoOptionError:
        ref_data = ref_scn.pdata[0]
    ref_filename = tempfile.mktemp(suffix='.nii')
    ref_data.export2nii(ref_filename)
    n_cols=ref_data.data.shape[2]

    fig_size,mod = determine_figure_size(n_rows,n_cols)
    
    title = k
    fig = pylab.figure(figsize=fig_size)
    fig.suptitle(k,fontsize=10+mod)
    G = pylab.matplotlib.gridspec.GridSpec(n_rows, n_cols, wspace=0.0, hspace=0.0)   
    print('\n'+'-'*80+'\n'+title)

    row_idx = 0
    for row in rows:
        row_conf = dict(config.items(row))
        lbl =row_conf.pop('label')
        subtitle = row_conf.pop('subtitle','')
        fname = v[lbl][0]
        bbox_pct = numpy.array([float(x) for x in v[lbl][1]])
        #bbox = sarpy.fmoosvi.getters.get_bbox(v,lbl)
        fig.add_subplot(G[row_idx, 0])
        pylab.axis('off')
        pylab.text(-0.65,0.5,'\n'.join([lbl,subtitle]), 
                 horizontalalignment='center', 
                 verticalalignment='center',
                 fontsize=5+mod, rotation='vertical')

        pylab.text(-0.35,0.5,'{0}'.format(fname), 
                 horizontalalignment='center', 
                 verticalalignment='center',
                 fontsize=2+mod, rotation='vertical')


        if fname == '':
            pylab.text(0.5,0.5,'Data not available',
                       horizontalalignment='center',
                       fontsize=4+mod)

            row_idx += 1
            continue

        # allowed types can be:
        # img, plot, vtc, histo

        if row_conf.get('type', None) == 'img':

            # Needed to flip colormaps for rois...default is False
            colormapFlip = False
            
            #TODO: this statement is here because of **row_conf in imshow
            row_conf.pop('type', None)
            clim_min = row_conf.pop('clim_min',None)
            clim_max = row_conf.pop('clim_max',None)

            scn = sarpy.Scan(fname)
            print(lbl, fname,scn.acqp.ACQ_protocol_name)            
            adata_key = row_conf.pop('adata', None)
            if adata_key is not None:

                if re.search('roi',adata_key):
                    colormapFlip = True                
                try:
                    data = scn.adata[adata_key]
                except KeyError:
                    pylab.text(0.5,0.5,'Data not available',
                       horizontalalignment='center',
                       fontsize=4+mod)

                    row_idx += 1
                    print('Adata Scan failed for {0},{1} \n \n'.format(k,adata_key))
                    continue
                    
                # Set the image limits for adata
                if (clim_min is None) and (clim_max is None):
                    (clim_min, clim_max) = sarpy.fmoosvi.getters.get_image_clims(
                                                data.data)
                else:
                    (clim_min, clim_max) = (numpy.int(clim_min),
                                            numpy.int(clim_max))
            else:
                data = scn.pdata[0]
                                
            xdata = data.data


            # Used masked arrays to show nan values as black

            xdata_mask = numpy.ma.array(xdata,mask=numpy.isnan(xdata))

            resample_flag = row_conf.pop('resample', False) 
            if resample_flag and config.getboolean(row,'resample'):
                raise NotImplementedError('please fix resampling')
                print('resampling {0}\n{1}\n onto {2}'.format(scn,data,ref_data))
            #src_filename = tempfile.mktemp(suffix='.nii')
            #data.export2nii(src_filename)            
            #xdata, xdata_sitk_image = sarpy.ImageProcessing.resample_onto.resample_onto(src_filename, ref_filename)
            #xdata = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(data, ref_data)
            #os.remove(src_filename)
            #find out where the 0 1 etc end up.
#            print(numpy.mean(numpy.mean(xdata, axis=0),axis=0))

            for col_idx in xrange(min(n_cols, xdata.shape[2])):
                fig.add_subplot(G[row_idx, col_idx])
                bbox = (bbox_pct.reshape(2,2).T*xdata.shape[0:2]).T.flatten()
                t=pylab.imshow(xdata_mask[bbox[0]:bbox[1],
                                 bbox[2]:bbox[3],col_idx],
                                 **row_conf)

                # Use black to show nan values
                # Source: http://stackoverflow.com/questions/2578752/how-can-i-plot-nan-values-as-a-special-color-with-imshow-in-matplotlib
                newcmap = copy.copy(t.get_cmap())

                if colormapFlip:
                    newcmap.set_bad('w',1.)
                else:
                    newcmap.set_bad('k',1.)
                t.set_cmap(newcmap)         
                t.set_clim(clim_min, clim_max)

                pylab.axis('off')
                
                if row_idx == 0:
                    pylab.title('{0}'.format(col_idx+1), fontsize=6+mod)
                    
            row_idx += 1
        
        elif row_conf.get('type', None) == 'vtc':
            
            scn = sarpy.Scan(fname)
            print(lbl, fname,scn.acqp.ACQ_protocol_name)            
            adata_key = row_conf.pop('adata', None)
            
            assert(adata_key is not None), 'Please supply a valid label for VTC'
            data = scn.adata[adata_key]
            reps = scn.pdata[0].data.shape[-1]

            for col_idx in xrange(min(n_cols, data.data.shape[2])):
                
                dat=scn.pdata[0].data[:,:,col_idx,:]
                imgdata = numpy.mean(dat,axis=2)
                vtcdata = data.data[:,:,col_idx]
                bbox = sarpy.fmoosvi.getters.get_roi_bbox(scn.shortdirname,'auc60_roi')
                axs=fig.add_subplot(G[row_idx, col_idx])
                                
                axs.imshow(imgdata[bbox[0]:bbox[1],\
                                   bbox[2]:bbox[3]],\
                                   cmap='gray', interpolation='none')
                axs.set_axis_off()
                fig.canvas.draw()
                
                box = axs.get_position().bounds
                
                height = box[3] / (bbox[1]-bbox[0])
               
                for i in xrange(bbox[0], bbox[1]):
                    tmpax = fig.add_axes([box[0], 
                                         box[1]+box[3]-(i-bbox[0]+1)*height, 
                                         box[2], height])
                    tmpax.set_axis_off()
                    tmpax.plot(vtcdata[i,(bbox[2]*reps):(bbox[3]*reps)],
                                       color='r', linewidth=.2)
                                       
        elif row_conf.get('type', None) == 'plot':

            scn = sarpy.Scan(fname)
            print(lbl, fname,scn.acqp.ACQ_protocol_name)            
            adata_key = row_conf.pop('adata', None)
            if adata_key is not None:
                
                try:
                    data = scn.adata[adata_key]
                except KeyError:

                    pylab.text(0.5,0.5,'Data not available',
                       horizontalalignment='center',
                       fontsize=4+mod)

                    row_idx += 1                    
                    print('Adata Scan failed for {0},{1} \n \n'.format(k,adata_key))
                    continue                
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
                pylab.xlabel('Time (ms)',fontsize=0+mod)
                pylab.ylabel('Enhancement (au)',fontsize=0+mod)                
                
                ax = pylab.gca()
                
                if row_idx == 0:
                    pylab.title('{0}'.format(col_idx+1), fontsize=6+mod)
                if col_idx != 0:
                    ax.axes.get_yaxis().set_visible([])
                    ax.axes.get_xaxis().set_visible([])
                
#TODO: finda way to put the axis labels on th INSIDE of the plot
                
                for label in ax.get_xticklabels():
                    label.set_fontsize(0+mod)

                for label in ax.get_yticklabels():
                    label.set_fontsize(0+mod)

            row_idx += 1

        if row_conf.get('type', None) == 'text':

            scn = sarpy.Scan(fname)
            print(lbl, fname,scn.acqp.ACQ_protocol_name)            
            adata_key = row_conf.pop('adata', None)
            
            if adata_key is not None:                
                try:
                    data = scn.adata[adata_key]
                except KeyError:
                    print('Adata Scan failed for {0},{1} \n \n'.format(k,adata_key))
                    continue
            else:
                raise IOError('No text adata find')
                                
            xdata = data.data
            
            subtitle = numpy.round(scipy.stats.nanmean(xdata))

            pylab.text(-0.15,0.5,'Avg: {0}'.format(subtitle), 
                             horizontalalignment='center', 
                             verticalalignment='center',
                             fontsize=5+mod, rotation='vertical')            
            

            for col_idx in xrange(min(n_cols, xdata.shape[0])):
                fig.add_subplot(G[row_idx, col_idx])
                pylab.axis('off')
                pylab.text(0.3,0.5,'{0}'.format(numpy.round(xdata[col_idx])), fontsize=8+mod)
                    
            row_idx += 1
        

















            
            
        elif row_conf.get('type', None) == 'histo':
            raise NotImplementedError('do not know how to draw histos')

        
    testPDF.savefig(fig)      
    
    # Saving Figure    
    #filename = os.path.expanduser(os.path.join('~/sdata',rootName[-2],args.output,k + '.png'))
    #pylab.savefig(filename, bbox_inches=0, dpi=300)
    pylab.close('all')

#os.remove(ref_filename)

testPDF.close()


