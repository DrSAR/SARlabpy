####### Imports #######

import ipywidgets
from ipywidgets import interact
import sarpy
import sarpy.analysis.cest
import pylab
import numpy
import sys

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

def browse_CEST(scn_to_view,adata_label,doFit = False, displayFit = False):
    
    scn = sarpy.Scan(scn_to_view)
    dataShape = scn.pdata[0].data.shape
        
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
        freqdata,zdata = sarpy.analysis.cest.process_cest(scn.shortdirname,x,y)
        
        
        pylab.plot(freqdata,zdata,'--',label='raw')
        pylab.xlim(5,-5)
        pylab.title('Curve', fontsize=18)

        pylab.legend()
        # Fit Data
        
        if displayFit is True:

            if doFit is True: 
                cestParams, acqfreqs, data = sarpy.analysis.cest.fit_px_cest(scn_to_view,x,y)
            else:
                cestParams = scn.adata[adata_label].data
            sarpy.analysis.cest.plotCestPeaks(cestParams,x,y)
          
    interact(view_image, x=ipywidgets.IntSlider(description='X-axis:',min=bbox[0],max=bbox[1],step=1),
                         y=ipywidgets.IntSlider(description='Y-axis:',min=bbox[2],max=bbox[3],step=1))	