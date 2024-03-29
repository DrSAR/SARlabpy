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
import configparser
import os
import matplotlib
matplotlib.use('Agg',warn=False)# where did I come from !?
from matplotlib_scalebar.scalebar import ScaleBar

import pylab
import numpy
#import sarpy
#import sarpy.analysis.getters
import tempfile
import scipy
import copy
import re
import scipy
import sarpy.analysis.getters
from scipy import ndimage
print(scipy.__version__)

#import sarpy.analysis.analysis
#import sarpy.analysis.colormaps as cmaps


from matplotlib.backends.backend_pdf import PdfPages

def generate(**kwargs):
    '''create tumourboard

    conf_file=None: name of config file as kw parameter
    test=False: Only  perform dry-run '''

    args = argparse.Namespace()
    # since argsvars will now be pointing at the attribute dictionary of
    # object args you will find below how changes in argsvars modifies the
    # content of args and vice versa. Convenient but maybe a little non-obvious.
    argsvars = vars(args)
    argsvars.update(kwargs)
    args.test = kwargs.get('test', False) # default for dry run

    if args.conf_file:
        print(("loading config file %s" % args.conf_file))
        config = configparser.ConfigParser()
        base_fname = os.path.join(os.path.expanduser('~'),
                                  'sdata',
                                  args.conf_file)
        if config.read([base_fname]):
            argsvars.update(dict(config.items("Defaults")))
        else:
            raise IOError('Could not read config file %s' % base_fname)

    if args.test:
        print('test mode:\n')
        print(args)   
    
    # determine the layout
    rows = [row for row in config.sections() if not((row=='Defaults'))]
    n_rows = len(rows)
    
    # Start a PDF file of all the animals
    
    exp_name = str(args.experiment_name)
    pdfName = os.path.splitext(str(args.conf_file).split('/')[-1])[0]

    pdfPath = os.path.expanduser(os.path.join('~/sdata',exp_name,args.output,pdfName+'.pdf'))
    testPDF = PdfPages(pdfPath)

    sepFiles = True
    import sarpy
 #   reload(sarpy) # this looks hacky. Is it?
    # if we really need this: importlib.reload for Python 3.4 and above
    # http://stackoverflow.com/questions/961162/reloading-module-giving-nameerror-name-reload-is-not-defined
    exp = sarpy.Experiment.from_masterlist(exp_name+'.config')        

    # for every patient we will create the same board
    for k in sorted(exp.patients.keys()):
        ref_lbl = config.get(args.ref_row,'label')

        try:
            ref_scan_name = exp.patients[k][ref_lbl]
            ref_scn = sarpy.Scan(ref_scan_name)
        except(AttributeError,IOError,KeyError):
            #print(('Ref Scan failed for {0},{1}. Trying alternate ref scan \n \n'.format(k,ref_lbl)))
            #$)#$I()
            ref_lbl = config.get(args.alternate_ref_row,'label')
            #raise NotImplementedError
            ref_scan_name = exp.patients[k][ref_lbl]
            ref_scn = sarpy.Scan(ref_scan_name)
        try:
            ref_data = ref_scn.adata[config.get(args.ref_row,'adata')].data
        except configparser.NoOptionError:
            ref_data = ref_scn.pdata[0].data
        except KeyError:
            ref_data = ref_scn.pdata[0].data

        #ref_filename = tempfile.mktemp(suffix='.nii')
        #ref_data.export2nii(ref_filename)
        n_cols=ref_data.shape[2]

        aspect = numpy.true_divide(n_rows,n_cols)
        fig_size = (n_cols*1.3, n_cols*aspect)    
        
        title = k
        fig = pylab.figure(figsize=fig_size)
        fig.suptitle(k)
        G = pylab.matplotlib.gridspec.GridSpec(n_rows, n_cols+1, wspace=0.0, hspace=0.0)   
        print(('\n'+'-'*80+'\n'+title))
    
        row_idx = 0
        for row in rows:
            row_conf = dict(config.items(row))
            lbl =row_conf.pop('label')

            try:
                lbl_scan_name = exp.patients[k][lbl]
            except KeyError:
                lbl_scan_name = ''

            subtitle = row_conf.pop('subtitle','')

            # added to account for different Bboxes for diff rois (E.g., IM and SC in OEP8)
            alternate_bbox = row_conf.pop('alternate_bbox','bbox')
            #try:
            #    alternate_bbox = row_conf.pop(args.alternate_bbox,'bbox')
            #except AttributeError:
            #    alternate_bbox = 'bbox'

            ax = fig.add_subplot(G[row_idx, 0])
            pylab.axis('off')
            ax.text(-.5,0.5,'\n'.join([lbl]), 
                     horizontalalignment='center', 
                     verticalalignment='center',
                     rotation='vertical',
                     transform=ax.transAxes)

            ax.text(-0.325,0.5,'\n'.join([subtitle]), 
                     horizontalalignment='center', 
                     verticalalignment='center',
                     rotation='vertical',
                     transform=ax.transAxes,
                     size='x-small')            

            ax.text(-0.15,0.5,'{0}'.format(lbl_scan_name), 
                    horizontalalignment='center', 
                    verticalalignment='center',
                    rotation='vertical',
                    transform=ax.transAxes,
                    size='xx-small')

            if lbl_scan_name == '':
                pylab.text(0.85,0.5,'Data not available',
                           horizontalalignment='center')
    
                row_idx += 1
                continue
    
            # allowed types can be:
            # img, plot, vtc, histo
    
            if row_conf.get('type', None) == 'img':
    
                # Needed to flip colormaps for rois...default is False
                colormapFlip = False
                
                #TODO: this statement is here because of **row_conf in imshow
                row_conf.pop('type', None)
                vmin = row_conf.pop('vmin',None)
                vmax = row_conf.pop('vmax',None)
                scn = sarpy.Scan(lbl_scan_name)
                adata_key = row_conf.pop('adata', None)

                    # try:
                    #     vmin = -0.2
                    #     vmax = 0.2
                    # except KeyError:
                    #     vmin = -10
                    #     vmax = 10

                if vmin is not None:
                    vmin = numpy.float(vmin)
                elif adata_key is not None and 'OEdraft' in adata_key: # Add limits specific to per scan for OE maps
                    vmin = -0.2

                if vmax is not None:
                    vmax = numpy.float(vmax)
                elif adata_key is not None and 'OEdraft' in adata_key: # Add limits specific to per scan for OE maps
                    vmax = 0.2
                ## Add an option to mask out the background 
                roi_mask = row_conf.pop('roi_adata', None)
                bat_adata = row_conf.pop('bat_adata', None)
                bat_threshold = row_conf.pop('bat_threshold', None)

                print(('\t {0}, {1}, {2}'.format(lbl, lbl_scan_name,adata_key)))

                if adata_key is not None:
    
                    if re.search('roi',adata_key):
                        colormapFlip = True                
                    try:
                        if roi_mask is None:
                            data = scn.adata[adata_key].data
                        else:
                            if bat_threshold is None:
                                data = scn.adata[adata_key].data * scn.adata[roi_mask].data
                            else:
                                BAT = scn.adata[bat_adata].data
                                bat_threshold = float(bat_threshold)
                                masked_threshold = sarpy.analysis.analysis.h_make_binary_mask(BAT,0,bat_threshold)
                                data = scn.adata[adata_key].data * scn.adata[roi_mask].data * masked_threshold
                    except KeyError:
                        pylab.text(0.85,0.5,'Data not available',
                           horizontalalignment='center')
    
                        row_idx += 1
                        print(('Adata Scan failed for {0},{1} \n \n'.format(k,adata_key)))
                        continue
                    except OSError:
                        pylab.text(0.85,0.5,'Data not available',
                           horizontalalignment='center')                        
                        # to deal with: sarpy.Scan('HerS10Do09.l11/4').adata['bbox']
                        row_idx += 1
                        print(('Adata Scan failed for {0},{1} \n \n'.format(k,adata_key)))
                        continue                        
                    # Set the image limits for adata
                    if (vmin is None) and (vmax is None):
                        (vmin, vmax) = sarpy.analysis.getters.get_image_clims(data)
                    else: #TODO This else statement (maybe) does NOTHING, get rid of it?
                        (vmin, vmax) = (vmin,vmax)

                    # When vmin is negative, it makes sense to invert the colorbar so that 
                    # the more negative values appear red and the values that are 0 appear blue
                    #
                    # See Github issue: https://github.com/DrSAR/SARlabpy/issues/275

                    #if vmin <= 0:
                    #    cm = row_conf.pop('type', 'jet')
                    #    cm = cm + '_r'

                else:
                    data = scn.pdata[0].data
                                    
                
                if len(data.shape)>3:
                    print('THIS SCAN DATA WAS chopped off; wrong dimensions, check!')
                    xdata = data[:,:,:,0]
                else:
                    xdata = data
                xdata_slices = sarpy.analysis.getters.get_num_slices(scn.shortdirname)
    
                # Used masked arrays to show nan values as black
    
                xdata_mask = numpy.ma.array(xdata,mask=numpy.isnan(xdata))
    
                resample_flag = row_conf.pop('resample', False) 
                if resample_flag and config.getboolean(row,'resample'):
                    raise NotImplementedError('please fix resampling')
                    print(('resampling {0}\n{1}\n onto {2}'.format(scn,data,ref_data)))
                #src_filename = tempfile.mktemp(suffix='.nii')
                #data.export2nii(src_filename)            
                #xdata, xdata_sitk_image = sarpy.ImageProcessing.resample_onto.resample_onto(src_filename, ref_filename)
                #xdata = sarpy.ImageProcessing.resample_onto.resample_onto_pdata(data, ref_data)
                #os.remove(src_filename)
                #find out where the 0 1 etc end up.
                #print(numpy.mean(numpy.mean(xdata, axis=0),axis=0))
    
                for col_idx in range(min(n_cols, xdata_slices)):
                    fig.add_subplot(G[row_idx,col_idx+1])
                    # Get the bbox as an adata 
                    try: 
                        bbox = scn.adata[alternate_bbox].data
                    except KeyError:
                        print('bbox import failed')
                        bbox = [0,scn.pdata[0].data.shape[0],0,scn.pdata[0].data.shape[1]]

                    if xdata_slices >1:     
                        if adata_key is not None and re.search('roi',adata_key):
                            # This image is to impose an ROI on top of an image
                            roi2 = scn.adata[adata_key].data
                            xdata_mask = scn.pdata[0].data

                            if len(roi2.shape)>3 or len(xdata_mask.shape)>3 :
                                print('THIS SCAN DATA WAS chopped off; wrong dimensions, check!')
                                xdata_mask = xdata_mask[:,:,:,0]

                            # Stefan's NEW method
                            sx = ndimage.sobel(numpy.nan_to_num(roi2[:,:,col_idx]), axis=0, mode='constant')
                            sy = ndimage.sobel(numpy.nan_to_num(roi2[:,:,col_idx]), axis=1, mode='constant')

                            sob = (numpy.hypot(sx, sy))[bbox[0]:bbox[1],bbox[2]:bbox[3]]

                            img_bgd = xdata_mask[bbox[0]:bbox[1],bbox[2]:bbox[3],col_idx]
                            img_norm = (img_bgd - numpy.min(img_bgd)) / (numpy.max(img_bgd)-numpy.min(img_bgd))
                            img_bgd_RGB = numpy.repeat(img_norm.reshape(img_norm.shape+(1,)), 3, axis=2)

                            red_channel = img_bgd_RGB[:,:,0]
                            red_channel[numpy.where(sob>2)] =1
                            img_bgd_RGB[:,:,0] = red_channel

                            t = pylab.imshow(img_bgd_RGB)    
                        else:
                            mapsToShow = xdata_mask[bbox[0]:bbox[1],bbox[2]:bbox[3],col_idx]

                            t=pylab.imshow(mapsToShow,**row_conf)
                            scalebar = ScaleBar(dx=scn.method.PVM_SpatResol[0]*1000,
                                                units='um',
                                                fixed_value=1,
                                                fixed_units='mm',
                                                location='lower left',label_formatter = lambda x, y:'',
                                                frameon=False,
                                                color='w',sep=-10,height_fraction=0.02)
                            pylab.gca().add_artist(scalebar)
                    else:
                        t=pylab.imshow(xdata_mask[bbox[0]:bbox[1],
                                         bbox[2]:bbox[3]],
                                         **row_conf)
                        scalebar = ScaleBar(dx=scn.method.PVM_SpatResol[0]*1000,
                                            units='um',
                                            fixed_value=1,
                                            fixed_units='mm',
                                            location='lower left',label_formatter = lambda x, y:'',
                                            frameon=False,
                                            color='w',sep=-10,height_fraction=0.02)
                        pylab.gca().add_artist(scalebar)
                    # Use black to show nan values
                    # Source: http://stackoverflow.com/questions/2578752/how-can-i-plot-nan-values-as-a-special-color-with-imshow-in-matplotlib
                    newcmap = copy.copy(t.get_cmap())
    
                    if colormapFlip:
                        newcmap.set_bad('w',1.)
                    else:
                        newcmap.set_bad('k',1.)
                    t.set_cmap(newcmap)         
                    t.set_clim(vmin, vmax)
    
                    pylab.axis('off')
                    
                    if row_idx == 0:
                        pylab.title('{0}'.format(col_idx+1))
                        
                row_idx += 1
            
            elif row_conf.get('type', None) == 'vtc':

                sepFiles = False
                
                scn = sarpy.Scan(lbl_scan_name)
                adata_key = row_conf.pop('adata', None)                
                print(('\t {0}, {1}, {2}'.format(lbl, lbl_scan_name,adata_key)))
                vtc_min = row_conf.pop('vtc_min',None)
                vtc_max = row_conf.pop('vtc_max',None)
                
                try:
                    data = scn.adata[adata_key]
                except KeyError:
                    pylab.text(0.85,0.5,'Data not available',
                       horizontalalignment='center')     
                    row_idx += 1
                    print(('Something failed in the vtc, cant get adata for scan {0}'.format(scn)))
                    continue

                reps = scn.pdata[0].data.shape[-1]
                # Set the image limits for the vtcs
                if (vtc_min is None) and (vtc_max is None):
                    vtc_min = 0
                    vtc_max = 1.5
                else:
                    vtc_min = numpy.float(vtc_min)
                    vtc_max = numpy.float(vtc_max)
    
                for col_idx in range(min(n_cols, xdata_slices)):

                    if xdata_slices ==1:
                        dat=scn.pdata[0].data[:,:,col_idx,:]
                        vtcdata = data.data[:,:,col_idx]

                    else: 
                        dat=scn.pdata[0].data[:,:,:]
                        vtcdata = data.data[:,:]

                    # Get the bbox as an adata 
                    bbox = scn.adata['bbox'].data

                    imgdata = numpy.mean(dat,axis=2)
                    axs=fig.add_subplot(G[row_idx,col_idx+1])

                                    
                    axs.imshow(imgdata[bbox[0]:bbox[1],\
                                       bbox[2]:bbox[3]],\
                                       cmap='gray', interpolation='none')
                    axs.set_axis_off()
                    fig.canvas.draw()
                    
                    box = axs.get_position().bounds
                    
                    height = box[3] / (bbox[1]-bbox[0])
                   
                    for i in range(bbox[0], bbox[1]):
                        tmpax = fig.add_axes([box[0], 
                                             box[1]+box[3]-(i-bbox[0]+1)*height, 
                                             box[2], height])
                        tmpax.set_axis_off()
                        tmpax.plot(vtcdata[i,(bbox[2]*reps):(bbox[3]*reps)],
                                           color='g', linewidth=.01)
                        pylab.ylim([vtc_min,vtc_max])
                row_idx += 1
                                           
            elif row_conf.get('type', None) == 'plot':
    
                scn = sarpy.Scan(lbl_scan_name)
                adata_key = row_conf.pop('adata', None)
                print(('\t {0}, {1}, {2}'.format(lbl, lbl_scan_name,adata_key)))
                if adata_key is not None:
                    
                    try:
                        data = scn.adata[adata_key]
                    except KeyError:
    
                        pylab.text(0.85,0.5,'Data not available',
                           horizontalalignment='center')
    
                        row_idx += 1                    
                        print(('Adata Scan failed for {0},{1} \n \n'.format(k,adata_key)))
                        continue                
                else:
                    data = scn.pdata[0]
                    
                xdata = data.data

                if ('OEdraft' in adata_key) and ('EC' in adata_key) :
                    fig.add_subplot(G[row_idx,1:])
                    pylab.plot(xdata)
                    pylab.xlabel('Reps')
                    pylab.ylabel('a.u')

                else:
    
                    for col_idx in range(min(n_cols, xdata.shape[0])):
                        fig.add_subplot(G[row_idx,col_idx+1])
                        pylab.plot(xdata[col_idx,0,:],xdata[col_idx,1,:])  
        
                        pylab.xlim(0,numpy.nanmax(xdata[:,0,:]))
                        pylab.ylim(0,numpy.nanmax(xdata[:,1,:]*1.3))
                        pylab.locator_params(axis='both',which='both',
                                             nbins=3) 
                        pylab.xlabel('Time (ms)',size='xx-small')
                        pylab.ylabel('Enhancement (au)',size='xx-small')                
                        
                        ax = pylab.gca()
                        
                        if row_idx == 0:
                            pylab.title('{0}'.format(col_idx+1))
                        if col_idx != 0:
                            ax.axes.get_yaxis().set_visible([])
                            ax.axes.get_xaxis().set_visible([])
                        
                    #TODO: finda way to put the axis labels on the INSIDE of the plot
                        
                        for label in ax.get_xticklabels():
                            label.set_fontsize(15)
        
                        for label in ax.get_yticklabels():
                            label.set_fontsize(15)
    
                row_idx += 1
    
            if row_conf.get('type', None) == 'text':
    
                scn = sarpy.Scan(lbl_scan_name)
                adata_key = row_conf.pop('adata', None)

                print(('\t {0}, {1}, {2}'.format(lbl, lbl_scan_name,adata_key)))
                
                if adata_key is not None:                
                    try:
                        data = scn.adata[adata_key]
                    except KeyError:
                        print(('Adata Scan failed for {0},{1} \n \n'.format(k,adata_key)))
                        continue
                else:
                    raise IOError('No text adata find')
                                    
                xdata = data.data
                
                subtitle = numpy.round(scipy.stats.nanmean(xdata))
    
                pylab.text(-0.15,0.5,'Avg: {0}'.format(subtitle), 
                                 horizontalalignment='center', 
                                 verticalalignment='center',
                                 rotation='vertical')            
    
                for col_idx in range(min(n_cols, xdata.shape[0])):
                    fig.add_subplot(G[row_idx,col_idx+1])
                    pylab.axis('off')
                    pylab.text(0.3,0.5,'{0}'.format(numpy.round(xdata[col_idx])))
                        
                row_idx += 1

        testPDF.savefig(fig)

        if sepFiles:
            # Saving Figure    
            filename = os.path.expanduser(os.path.join('~/sdata',exp_name,args.output,k + '.png'))
            pylab.savefig(filename, bbox_inches='tight', dpi=300)
            pylab.close('all')
    
    testPDF.close()
    pylab.clf()
    import time
    print('File is opened and ready to use')
    print((time.strftime("%c")))

    
if __name__ == "__main__":
    # we are being run from the commandline
    conf_parser = argparse.ArgumentParser(
        # Turn off help, so we print all options in response to -h
            description=__doc__)
    conf_parser.add_argument("-c", "--conf_file",
                             help="Specify full path to config file", metavar="FILE")
    conf_parser.add_argument("-t", "--test",default=False, action="store_true",
                        help='dry run of mriBoard creation')

    #parse_known_args does not produce an error for unknown args
    args = conf_parser.parse_args()

    # at this point args is a 'Namespace' which is defined in the argparse module
    # to get access to its attributes you can simply say, e.g., args.conf_file
    # to treat it like a dictionary you have to use d=vars(args). Then you can
    # read and write to it: d['newattribute']=42 at which point you can verify
    # that args.newatttribute exists and equals 42. The magic of shallow copies...
    generate(**vars(args))


