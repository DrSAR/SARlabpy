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

def browse_CEST(scn_to_view,adata_label=None,doFit = False, displayFit = False):
    
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
        pylab.axvline(0,linestyle='--')

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

def browse_LLfit(llscn_name,filtereddata):
    
    scn = sarpy.Scan(llscn_name)
    dataShape = scn.pdata[0].data.shape

    try:
        bbox = scn.adata['bbox'].data
    except:
        bbox = [0,dataShape[0]-1,0,dataShape[1]-1]    
    
    def view_image(x, y,slc):
        pylab.figure(figsize=(16,5))
        
        #### Create the image
        pylab.subplot(121)
        pylab.imshow(scn.pdata[0].data[:,:,slc,0],interpolation='None')
        pylab.axvline(y,color='w')
        pylab.axhline(x,color='w')
        
        pylab.ylim(bbox[1],bbox[0])
        pylab.xlim(bbox[2],bbox[3])
        
        #### Fit fit the data then plot
        pylab.subplot(122)
        
        # do the fit
        infodict,mesg,ier,fit_params, T1,t_data = sarpy.analysis.analysis.h_fitpx_T1_LL_FAind(
                                                                  scn_to_analyse=llscn.shortdirname,
                                                                  y_data=filtereddata[x,y,slc,:],
                                                                  slc=slc)        
        # Stored Fit Data
        storedDict = scn.adata['T1_LL_fit'].data
        stored_fitparams = [storedDict['params']['a'][x,y,slc],
                            storedDict['params']['b'][x,y,slc],
                            storedDict['params']['T1_eff'][x,y,slc],
                            storedDict['params']['phi'][x,y,slc]]
        stored_T1 = scn.adata['T1_LL'].data[x,y,slc]

        #### Plot the LL data
        
        # Raw Data      
        ydata = numpy.real(scn.fftfid[x,y,slc,:])
        filtered_y = numpy.real(filtereddata[x,y,slc,:])
        pylab.plot(t_data,ydata,'o',label='raw data')
        pylab.plot(t_data,numpy.real(filtered[x,y,slc,:]),'x',label='filtered data')
        pylab.title('Traces', fontsize=18)

        # New Fit Data
        filtered_fit = sarpy.analysis.analysis.h_func_T1_FAind(fit_params,t_data)
        pylab.plot(t_data,filtered_fit,label='filtered fit')
        pylab.axvline(stored_T1)
        
        # Old Fit Data
        stored_fit = sarpy.analysis.analysis.h_func_T1_FAind(stored_fitparams,t_data)
        pylab.plot(t_data,stored_fit,label='Original fit')
        pylab.axvline(T1)
        pylab.legend()        
                
          
    interact(view_image, x=ipywidgets.IntSlider(description='X-axis:',min=bbox[0],max=bbox[1],step=1),
                         y=ipywidgets.IntSlider(description='Y-axis:',min=bbox[2],max=bbox[3],step=1),
                         slc=ipywidgets.IntSlider(description='Slice:',min=0,max=dataShape[2],step=1))