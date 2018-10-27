# -*- coding: utf-8 -*-
# Copyright (C) 2012-2013 Stefan A Reinsberg and SARlab members
# full license details see LICENSE.txt
"""Collection of analysis routines for Bruker data
test

"""
import numpy
import sarpy
import sarpy.analysis.getters as getters
import scipy.integrate
import scipy.optimize
import scipy.fftpack
import scipy.stats
import copy
import os
import json
import datetime
import collections
import imp
import lmfit

from lmfit.models import ConstantModel, LorentzianModel



################################################

####### Fitting & Processing Functions #######    

################################################

# def h_old_lorentzian(A,w,p,freqs):
#     return numpy.divide(A,(1+4*((freqs-p)/w)**2))

def h_lorentzian(amplitude,sigma,center,freqs=None):
    if freqs is None:
        freqs=numpy.arange(-50,50,0.1)
    return freqs,(amplitude/numpy.pi)*numpy.divide(sigma,(freqs-center)**2 + sigma**2)

# def h_zspectrum_N(params,freqs):
#     ''' Updated Zspectrum N function to now require a Params object'''

#     arr = numpy.zeros_like(freqs)

#     # First get the global DC offset
#     DCoffset =  params['DCoffset']

#     # Now get the other peaks as regular lorentzians
#     for i in numpy.arange(0,int(len(params.values())/3)):
#         A = params['A_%i' % i]
#         w = params['w_%i' % i]
#         p = params['p_%i' % i]

#         newpeak = h_lorentzian(A,w,p,freqs)

#         if numpy.isnan(numpy.sum(newpeak)):
#             raise ValueError('Values for Lorentzian gives problems {} {} {}'.format(A,w,p))
#         arr += newpeak

#     return 100-arr + DCoffset

def fit_water_peak(data,offset_freqs,shiftWaterPeak = False):

    """
    Fits a Lorentzian to the data, and 
    returns the water offset frequency
    """
    ## First set the weights so that the ones around the peak are the most important
    # Weights 
    # Get the index at the maximum value which represents 0ppm roughly
    pk_1_idx = [i for i, j in enumerate(data) if j == max(data)][0]

    # create an array of 0s the right shape
    waterweights = numpy.array(offset_freqs.copy())*0

    # Set the 6 points around 0 to be utmost importance 100, everything else to be 1
    waterweights[pk_1_idx-3:pk_1_idx+3] = 100
    waterweights[pk_1_idx-10:pk_1_idx+10] = 1

    # In case you only want to give it a subset of data:
    xtofit = offset_freqs[pk_1_idx-10:pk_1_idx+10]
    ytofit = data[pk_1_idx-10:pk_1_idx+10]
    waterweightstofit = waterweights[pk_1_idx-10:pk_1_idx+10]

    # Set up the model
    lmodel = ConstantModel()+LorentzianModel(prefix='p1_')

    # Set up the initial parameters
    pars = lmodel.make_params(c=5)

    isigma = [1, 1]
    sigma_min = [0,0,]
    sigma_max = [2000,2000]

    icentre = [-.5, 0]
    centre_min = [-2,-4]
    centre_max = [2,4]

    iamplitude = [300, 130]
    amplitude_min = [0,0]
    amplitude_max = [5000,5000]

    for i in range(1):
        pref = 'p%i_' % (i+1)
        amp, cen, wid = iamplitude[i], icentre[i], isigma[i]
        amp_min, cen_min, wid_min = amplitude_min[i], centre_min[i], sigma_min[i]
        amp_max, cen_max, wid_max = amplitude_max[i], centre_max[i], sigma_max[i]
        
        pars['%scenter' % pref].set(value=cen, min=cen_min, max=cen_max)    
        pars['%samplitude' % pref].set(value=amp, min=amp_min, max=amp_max)
        pars['%ssigma'% pref].set(value=wid, min=wid_min, max=wid_max)
    
    # Perform the fit    
    out  = lmodel.fit(ytofit, pars, x=xtofit,weights = waterweightstofit)
    pkloc = out.best_values['p1_center']

    if shiftWaterPeak is False:
        return pkloc
    else:
        shifted_data = numpy.interp(numpy.array(offset_freqs), numpy.array(offset_freqs)-pkloc, numpy.array(data))
        return pkloc, shifted_data


def fit_px_cest(scn_to_analyse, xval, yval):

    # Get the processed data

    freqdata,zdata,wateroffset = sarpy.analysis.cest.process_cest(scn_to_analyse,xval,yval,shiftWaterPeak=True)

    # Set up the model

    lmodel =  (LorentzianModel(prefix='p1_')
              + LorentzianModel(prefix='p2_')  + LorentzianModel(prefix='p3_')
              + LorentzianModel(prefix='p4_')  + LorentzianModel(prefix='p5_'))

    pars = lmodel.make_params()

    icentre =    [-3.2,-1.5,0,2.0,3.6]
    centre_min = [ic - 0.4 for ic in icentre]
    centre_max = [ic + 0.4 for ic in icentre]

    isigma =    [2,25,0.5,1, 1,1.5]
    sigma_min = [0,0,0,0,0]
    sigma_max = [8,100,5,2,2]

    iamplitude =    [50,500,150,30,20]
    amplitude_min = [0,10,0,0,0]
    amplitude_max = [200,1000,200,100,100]

    for i in range(5):
        pref = 'p%i_' % (i+1)
        amp, cen, wid = iamplitude[i], icentre[i], isigma[i]
        amp_min, cen_min, wid_min = amplitude_min[i], centre_min[i], sigma_min[i]
        amp_max, cen_max, wid_max = amplitude_max[i], centre_max[i], sigma_max[i]
        
        pars['%scenter' % pref].set(value=cen, min=cen_min, max=cen_max)    
        pars['%samplitude' % pref].set(value=amp, min=amp_min, max=amp_max)
        pars['%ssigma'% pref].set(value=wid, min=wid_min, max=wid_max)

    out  = lmodel.fit(zdata, pars, x=freqdata,fit_kws={'maxfev': 50000})

    return out

def fit_cest(scn_to_analyse=None,roi_adata_label = None):

    scan_object = sarpy.Scan(scn_to_analyse)

    # Data
    cestdata = scan_object.pdata[0].data

    x_size = cestdata.shape[0]
    y_size = cestdata.shape[1]
    
    # Deal with bounding boxes
    ## Check for bbox traits and create bbox_mask to output only partial data

    try:
        bbox = scan_object.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])

    if roi_adata_label:
        roi = scan_object.adata[roi_adata_label].data
    else:
        roi = numpy.empty_like(cestdata[:,:,0])*0 + 1

    # Fit multiple peaks, need some empty arrays
    p1_amplitude = numpy.empty_like(roi) + numpy.nan
    p1_center = numpy.empty_like(roi) + numpy.nan
    p1_sigma = numpy.empty_like(roi) + numpy.nan

    p2_amplitude = numpy.empty_like(roi) + numpy.nan
    p2_center = numpy.empty_like(roi) + numpy.nan
    p2_sigma = numpy.empty_like(roi) + numpy.nan

    p3_amplitude = numpy.empty_like(roi) + numpy.nan
    p3_center = numpy.empty_like(roi) + numpy.nan
    p3_sigma = numpy.empty_like(roi) + numpy.nan

    p4_amplitude = numpy.empty_like(roi) + numpy.nan
    p4_center = numpy.empty_like(roi) + numpy.nan
    p4_sigma = numpy.empty_like(roi) + numpy.nan

    p5_amplitude = numpy.empty_like(roi) + numpy.nan
    p5_center = numpy.empty_like(roi) + numpy.nan
    p5_sigma = numpy.empty_like(roi) + numpy.nan

    fit_redchi = numpy.empty_like(roi) + numpy.nan
    fit_report = numpy.empty_like(roi,dtype='object')

    newstruct = numpy.zeros(roi.shape, dtype=[
       ('p1_amplitude', 'float64'),('p1_sigma', 'float64'),('p1_center', 'float64'),
       ('p2_amplitude', 'float64'),('p2_sigma', 'float64'),('p2_center', 'float64'),
       ('p3_amplitude', 'float64'),('p3_sigma', 'float64'),('p3_center', 'float64'),
       ('p4_amplitude', 'float64'),('p4_sigma', 'float64'),('p4_center', 'float64'),
       ('p5_amplitude', 'float64'),('p5_sigma', 'float64'),('p5_center', 'float64')])

    # Nan the array so there are no zeroes anywhere
    newstruct[:] = numpy.nan

    for xval in range(0,x_size):    
        for yval in range(0,y_size):
            if roi[xval,yval] == 1:
                
                try:
                    output = fit_px_cest(scn_to_analyse,xval,yval)
                except (ValueError,TypeError) as e:
                    continue
                p1_amplitude[xval,yval] = output.best_values['p1_amplitude']
                p1_sigma[xval,yval] = output.best_values['p1_sigma']
                p1_center[xval,yval] = output.best_values['p1_center']

                p2_amplitude[xval,yval] = output.best_values['p2_amplitude']
                p2_sigma[xval,yval] = output.best_values['p2_sigma']
                p2_center[xval,yval] = output.best_values['p2_center']

                p3_amplitude[xval,yval] = output.best_values['p3_amplitude']
                p3_sigma[xval,yval] = output.best_values['p3_sigma']
                p3_center[xval,yval] = output.best_values['p3_center']

                p4_amplitude[xval,yval] =  output.best_values['p4_amplitude']
                p4_sigma[xval,yval] = output.best_values['p4_sigma']
                p4_center[xval,yval] =  output.best_values['p4_center']

                p5_amplitude[xval,yval] =  output.best_values['p5_amplitude']
                p5_sigma[xval,yval] =  output.best_values['p5_sigma']
                p5_center[xval,yval] =  output.best_values['p5_center']
                
                fit_redchi[xval,yval] = output.redchi

                fit_report[xval,yval] = output.fit_report()

    # Save the data as a structured array

    newstruct['p1_amplitude'] = p1_amplitude
    newstruct['p1_sigma'] = p1_sigma
    newstruct['p1_center'] = p1_center
    newstruct['p2_amplitude'] = p2_amplitude
    newstruct['p2_sigma'] = p2_sigma
    newstruct['p2_center'] = p2_center
    newstruct['p3_amplitude'] = p3_amplitude
    newstruct['p3_sigma'] = p3_sigma
    newstruct['p3_center'] = p3_center
    newstruct['p4_amplitude'] = p4_amplitude
    newstruct['p4_sigma'] = p4_sigma
    newstruct['p4_center'] = p4_center
    newstruct['p5_amplitude'] = p5_amplitude
    newstruct['p5_sigma'] = p5_sigma
    newstruct['p5_center'] = p5_center

    return {'':newstruct,'fit_redchi':fit_redchi,'fit_report':fit_report}

def h_fitoutput_to_struct(fitoutput):

    newstruct = numpy.zeros((1), dtype=[
   ('p1_amplitude', 'float64'),('p1_sigma', 'float64'),('p1_center', 'float64'),
   ('p2_amplitude', 'float64'),('p2_sigma', 'float64'),('p2_center', 'float64'),
   ('p3_amplitude', 'float64'),('p3_sigma', 'float64'),('p3_center', 'float64'),
   ('p4_amplitude', 'float64'),('p4_sigma', 'float64'),('p4_center', 'float64'),
   ('p5_amplitude', 'float64'),('p5_sigma', 'float64'),('p5_center', 'float64')])

    newstruct['p1_amplitude'] = fitoutput.best_values['p1_amplitude']
    newstruct['p1_sigma'] = fitoutput.best_values['p1_sigma']
    newstruct['p1_center'] = fitoutput.best_values['p1_center']

    newstruct['p2_amplitude'] = fitoutput.best_values['p2_amplitude']
    newstruct['p2_sigma'] = fitoutput.best_values['p2_sigma']
    newstruct['p2_center'] = fitoutput.best_values['p2_center']

    newstruct['p3_amplitude'] = fitoutput.best_values['p3_amplitude']
    newstruct['p3_sigma'] = fitoutput.best_values['p3_sigma']
    newstruct['p3_center'] = fitoutput.best_values['p3_center']

    newstruct['p4_amplitude'] =  fitoutput.best_values['p4_amplitude']
    newstruct['p4_sigma'] = fitoutput.best_values['p4_sigma']
    newstruct['p4_center'] =  fitoutput.best_values['p4_center']

    newstruct['p5_amplitude'] =  fitoutput.best_values['p5_amplitude']
    newstruct['p5_sigma'] =  fitoutput.best_values['p5_sigma']
    newstruct['p5_center'] =  fitoutput.best_values['p5_center']

    return newstruct

def process_cest(scn_to_analyse, xval, yval, shiftWaterPeak = False, pdata_num = 0):
    
    scn = sarpy.Scan(scn_to_analyse)
    
    # Defining parameters
    freq_list = scn.method.CEST_FreqListPPM

    ppm_limit_min = -55
    ppm_limit_max = 55
    exclude_ppm = 200
    normalize_to_ppm = 66.6

    possibleNormalizations = [i for i, x in enumerate(freq_list) if numpy.abs(x - normalize_to_ppm) <1E-4]
    normalizeTo = numpy.nonzero([scn.pdata[0].data[xval,yval,pn] for pn in possibleNormalizations])[0][-1] # Return the last non-zero value to make sure  data is sensible

    # Get only the frequencies within the ppm_limit
    ppm_filtered = [f for f in freq_list if ppm_limit_max > f > ppm_limit_min]

    # Exclude the dummy frequencies at the beginning (66.6 ppm). Sort for prettier plotting
    ppm_filtered = sorted([n for n in ppm_filtered if n!= exclude_ppm])

    # Get the index of the good frequencies relative to the original list 
    ppm_filtered_ind = [freq_list.index(c) for c in ppm_filtered]  

    # Normalize the data
    tmp = 100*(1 - scn.pdata[0].data[xval,yval,:][ppm_filtered_ind] / scn.pdata[0].data[xval,yval,normalizeTo])

    if shiftWaterPeak is False:
        return ppm_filtered, tmp
    else:
        waterOffset, shiftedData = fit_water_peak(tmp,ppm_filtered,shiftWaterPeak = True)
        return ppm_filtered,shiftedData,waterOffset

def find_nearest(array, value):
    array = numpy.asarray(array)
    idx = (numpy.abs(array - value)).argmin()
    return idx,array[idx]
################################################

####### Displaying and Plotting Functions #######    

################################################

def fitsanity(out):
    
    properPeaks = ['NOE (-3.2)','MT (-1.5)','Water (0)','Amine (2.0)','Amide (3.6)']
   
    # first get all the centres so we can figure out which peak is which
    names = ['p{0}'.format(i)+'_' for i in range(1,6)]
    actualcentres = [-3.2,-1.5,0,2,3.6]
    centres = [out.values[n+'center'] for n in names]

    # this gets the actual idx of the peaks in case they move around
    actualidx = []
    for ctr in actualcentres:
        a,b = find_nearest(centres,ctr)
        actualidx.append(a+1)

    # Redo getting the names and centres so that the peak order is maintained
    names = ['p{0}'.format(i)+'_' for i in actualidx]
    centres = [out.values[n+'center'] for n in names]
    amps = [out.values[n+'amplitude'] for n in names]
    sigmas = [out.values[n+'sigma'] for n in names]                  

    for i,pk in enumerate(properPeaks):
                   
        print('=======',pk,'=======')
        print('Centre: {0:.2f} \nsigma {1:.2f} \t \nAmplitude {2:.2f}'.format(centres[i],sigmas[i],amps[i]))
    print('\n')

def plotPeaks(output,peaksToPlot = [1,2,3,4,5]):

    import pylab
    
    print(peaksToPlot)
    sumpk = numpy.array([0]*1000,dtype=object)
    for count in peaksToPlot:       
        peaknum = 'p'+str(count)+'_'

        freqs,tmppk = sarpy.analysis.cest.h_lorentzian(output[peaknum+'amplitude'],
                                                 output[peaknum+'sigma'],
                                                 output[peaknum+'center'])
        pylab.plot(freqs, tmppk,label=peaknum[0:-1])

        sumpk +=tmppk

    pylab.axvline(2.0,ymin=0,label='2.0 amine',color='y', alpha=0.4)
    pylab.axvline(3.6,ymin=0,label='3.6 amide',color='r', alpha=0.4)
    pylab.axvline(-3.2,ymin=0,label='-3.2 NOE',color='g', alpha=0.4)
    pylab.axvline(-1.5,ymin=0,label='MT')    
    pylab.legend()
    pylab.plot(freqs,sumpk,label='Sum')
    #pylab.axvline(1.5,ymin=0.6,label='1.5 OH',color='b', alpha=0.4)
    return freqs,sumpk

def generate_offset_list(additionalDict = None,
                         manuallyInsertedOffsets = None,
                         manuallyInsertedPositions = None,
                         alternateFreqs = True):

    if additionalDict is None:
        additionalDict = collections.OrderedDict([
                               ('dummy',[-60,60,10]),
                               ('start',[-5.1,5.1,0.1]),
                               ('baseline',[-60,60,10]),       
                               ('2.5',[2.2,2.5,0.01]),
                               ('3.4',[3.3,3.5,0.01])
                              ])

    offsetList = []

    for k,v in list(additionalDict.items()):

        offsetList.extend(numpy.round(numpy.arange(v[0],
                                                   v[1],
                                                   v[2]),3))

    # Reorder the list so it's alternating
    # (largest -> smallest -> next largest -> next smallest ->)
    # Gotten from: http://stackoverflow.com/questions/17436870/python-alternating-elements-of-a-sorted-array
    # Of course :-)

    if alternateFreqs is True:
        offsetList = list(sum(list(zip(reversed(offsetList), offsetList)), ())[:len(offsetList)])

    # Now manually insert offsets and frequency
    if manuallyInsertedOffsets is None:
        print([numpy.float("{:.3f}".format(off)) for off in offsetList])
        return numpy.round(offsetList,3)
    else:

        assert len(manuallyInsertedPositions) == len(manuallyInsertedOffsets), "List lengths not the same, please check input lists"
        
        for off,pos in zip(manuallyInsertedOffsets,manuallyInsertedPositions):

            offsetList.insert(pos,off)

        print([numpy.float("{:.3f}".format(off)) for off in offsetList])
        return numpy.round(offsetList,3)          