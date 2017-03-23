####### Imports #######

import ipywidgets
from ipywidgets import interact
import sarpy
import pylab
import numpy
import sys
sys.path.append('/home/fmoosvi/Desktop/Dropbox/code/python-cest')
import cest
import cest.analysis
####### Visualization Functions #######

def browse_MRimages(data):
	'''
	This function takes in the data and simply plots it with some slider bars for
	x,y,slice and puts a cursor location at the pixel.
	'''

	import pylab

	datashape = data.shape
	def view_image(x, y, i):
		pylab.figure(figsize=(16,5))

		pylab.subplot(121)
		pylab.imshow(data[:,:,i,0])
		pylab.axvline(x,color='w')
		pylab.axhline(y,color='w')

		pylab.subplot(122)
		pylab.plot(data[y,x,i,:],'--')
		pylab.title('Curve', fontsize=18)

	interact(view_image, x=ipywidgets.IntSlider(description='X-axis:',min=0,max=datashape[1]-1,step=1),
                         y=ipywidgets.IntSlider(description='Y-axis:',min=0,max=datashape[0]-1,step=1),
                         i=ipywidgets.IntSlider(description='Slice:',min=0,max=datashape[2]-1,step=1))

def browse_CEST(scn_to_view,adata_label):
    
    scn = sarpy.Scan(scn_to_view)
    cestParams = scn.adata[adata_label].data
    dataShape = cestParams.shape
        
    try:
        bbox = scn.adata['bbox'].data
    except:
        bbox = [0,dataShape[0]-1,0,dataShape[1]-1]    
    
    def view_image(x, y):
        pylab.figure(figsize=(16,5))
        
        #### Create the image
        pylab.subplot(121)
        pylab.imshow(scn.pdata[0].data[:,:,0],interpolation='None')
        pylab.axvline(y,color='w')
        pylab.axhline(x,color='w')
        
        pylab.ylim(bbox[1],bbox[0])
        pylab.xlim(bbox[2],bbox[3])
        
        #### Plot the CEST Spectrum
        pylab.subplot(122)
        # Raw Data
        freqdata,zdata = cest.analysis.cest_spectrum(scn.shortdirname,
                                                      x,
                                                      y,
                                                      shift_water_peak = True,
                                                      normalize=True,
                                                      normalize_to_ppm = 66.6,
                                                      ppm_limit_min = -50,
                                                      ppm_limit_max = 50,
                                                      exclude_ppm = 66.6,
                                                      pdata_num = 0)
        
        
        pylab.plot(freqdata,zdata,'--',label='raw')
        pylab.xlim(10,-10)
        pylab.title('Curve', fontsize=18)
        pylab.axvline(2.2,ymin=0.6,label='2.2 amine',color='y', alpha=0.4)
        pylab.axvline(3.5,ymin=0.6,label='3.5 amide',color='r', alpha=0.4)
        pylab.axvline(1.5,ymin=0.6,label='1.5 OH',color='b', alpha=0.4)
        pylab.axvline(-3.25,ymin=0.6,label='-3.25 aliphatic',color='g', alpha=0.4)
        #pylab.axvline(-3.0,label='-3.0 ppm-amine')
        pylab.legend()
        # Fit Data
        w = numpy.arange(-20.,20.,0.01) # Sample frequencies
        
        fitteddata = cest.analysis.h_zspectrum_New(cestParams[x,y],w)
        pylab.plot(w,fitteddata,label='Fit')
        pylab.xlim(15,-15)
        
        
        
    interact(view_image, x=ipywidgets.IntSlider(description='X-axis:',min=bbox[0],max=bbox[1],step=1),
                         y=ipywidgets.IntSlider(description='Y-axis:',min=bbox[2],max=bbox[3],step=1))	