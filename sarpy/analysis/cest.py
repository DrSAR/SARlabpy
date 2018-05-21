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

################################################

####### Fitting & Processing Functions #######    

################################################

def h_lorentzian(A,w,p,freqs):
    return numpy.divide(-A,(1+4*((freqs-p)/w)**2))

def h_zspectrum_N(params,freqs):
    ''' Updated Zspectrum N function to now require a Params object'''

    arr = numpy.zeros_like(freqs)

    # First get the global DC offset
    DCoffset =  params['DCoffset']

    # Now get the other peaks as regular lorentzians
    for i in numpy.arange(0,int(len(params.values())/3)):
        A = params['A_%i' % i]
        w = params['w_%i' % i]
        p = params['p_%i' % i]

        newpeak = h_lorentzian(A,w,p,freqs)

        if numpy.isnan(numpy.sum(newpeak)):
            raise ValueError('Values for Lorentzian gives problems {} {} {}'.format(A,w,p))
        arr += newpeak

    return 100-arr + DCoffset

def h_zspectrum_N_residuals(params,freqs, data):
    return h_zspectrum_N(params, freqs) - data


def h_peak_N(params,freqs,peak):
    ''' Zspectrum N function to only return one peak'''

    # First get the water peak as a super lorentzian
    #DCoffset =  params['DCoffset']
    DCoffset =  params['DCoffset']
    A = params['A_%i' % peak]
    w = params['w_%i' % peak]
    p = params['p_%i' % peak]

    ########################################
    #  Need to fix this.
    return h_lorentzian(A,w,p,freqs) + DCoffset   

def fit_water_peak(data,offset_freqs,allParams = False):

    """
    Fits a Super Lorentzian to the data, and 
    returns the water offset frequency

    """
    params = lmfit.Parameters()

    # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)    
    params.add_many(('DCoffset', 1, True, None, None, None),
                    ('A_0', 0.8, True, None, None, None),
                    ('w_0', 1.3, True, None, None, None),
                    ('p_0', 0.01, True, None, None, None),
    )        

    try:
        out = lmfit.minimize(h_zspectrum_N_residuals, params, args=(numpy.array(offset_freqs), data))
        return numpy.array(out.params['p_0'].value)
    except ValueError: # when the data given is nonsense
        return numpy.array(0)

def fit_px_cest(scn_to_analyse, xval, yval, pdata_num = 0):
    scn = sarpy.Scan(scn_to_analyse)
    params = lmfit.Parameters()
    params.add('DCoffset',  value = 1)

    ampls = [.84, .11, .09, .17, .05]
    amplsmax = [0.9, 0.3,0.3,0.3,0.3]
    width = [1.1, 1.98, 1.34, 4.2, 1.2]

    peaks =    [.05, 2.2, 3.6, -3.67, -3.86]
    peaksmin = [-1.5, 1.8, 2.5, -3.8, -4]
    peaksmax = [1, 2.5, 5, -2, -3]

    for i in numpy.arange(0,5):
            # add with tuples: (NAME VALUE VARY MIN  MAX  EXPR  BRUTE_STEP)    
            params.add_many(('A_%i' % i, ampls[i], True, 0, amplsmax[i], None),
                            ('w_%i' % i, width[i], True, None, None, None),
                            ('p_%i' % i, peaks[i], True, None, None, None),
            )    

    acqfreqs, data = process_cest(scn.shortdirname,xval,yval)

    fitoutput = lmfit.minimize(h_zspectrum_N_residuals, params, args=(numpy.array(acqfreqs), data))
    
    return fitoutput,acqfreqs, data 

def fit_5_peaks_cest(scn_to_analyse):

    scn = sarpy.Scan(scn_to_analyse)
    pdata_num = 0

    try:
        roi = scn.adata['roi'].data
    except KeyError:
        roi = scn.pdata[0].data[:,:,0]*0+1

    # Get the bbox so that the whole image isn't fit
    try:
        bbox = scn.adata['bbox'].data
    except KeyError:       
        bbox = numpy.array([0,x_size-1,0,y_size-1])   

    datashape = roi.shape
    roi_reshaped = numpy.reshape(roi,[roi.shape[0],roi.shape[1],1])
    cestscan_roi = scn.pdata[0].data * roi_reshaped       

    # Fit multiple peaks, need some empty arrays
    offst = numpy.empty_like(roi) + numpy.nan
    pk1_amp = numpy.empty_like(roi) + numpy.nan
    pk1_pos = numpy.empty_like(roi) + numpy.nan
    pk1_width = numpy.empty_like(roi) + numpy.nan

    pk2_amp = numpy.empty_like(roi) + numpy.nan
    pk2_pos = numpy.empty_like(roi) + numpy.nan
    pk2_width = numpy.empty_like(roi) + numpy.nan

    pk3_amp = numpy.empty_like(roi) + numpy.nan
    pk3_pos = numpy.empty_like(roi) + numpy.nan
    pk3_width = numpy.empty_like(roi) + numpy.nan

    pk4_amp = numpy.empty_like(roi) + numpy.nan
    pk4_pos = numpy.empty_like(roi) + numpy.nan
    pk4_width = numpy.empty_like(roi) + numpy.nan

    water_amp = numpy.empty_like(roi) + numpy.nan
    water_pos = numpy.empty_like(roi) + numpy.nan
    water_width = numpy.empty_like(roi) + numpy.nan

    fit_quality = numpy.empty_like(roi) + numpy.nan

    newstruct = numpy.zeros(roi.shape, dtype=[('DCoffset', 'float64'),
       ('A_0', 'float64'),('w_0', 'float64'),('p_0', 'float64'),
       ('A_1', 'float64'),('w_1', 'float64'),('p_1', 'float64'),
       ('A_2', 'float64'),('w_2', 'float64'),('p_2', 'float64'),
       ('A_3', 'float64'),('w_3', 'float64'),('p_3', 'float64'),
       ('A_4', 'float64'),('w_4', 'float64'),('p_4', 'float64')])

    # Nan the array so there are no zeroes anywhere
    newstruct[:] = numpy.nan

    for xval in range(bbox[0],bbox[1]):    
        for yval in range(bbox[2],bbox[3]):
            
            output, acqfreqs, data = fit_px_cest(scn.shortdirname,xval,yval)

            offst[xval,yval] = output.params['DCoffset'].value

            water_amp[xval,yval] =  output.params['A_0'].value
            water_width[xval,yval] =  output.params['w_0'].value
            water_pos[xval,yval] =  output.params['p_0'].value

            pk1_amp[xval,yval] = output.params['A_1'].value
            pk1_width[xval,yval] = output.params['w_1'].value
            pk1_pos[xval,yval] = output.params['p_1'].value

            pk2_amp[xval,yval] = output.params['A_2'].value
            pk2_width[xval,yval] = output.params['w_2'].value
            pk2_pos[xval,yval] = output.params['p_2'].value

            pk3_amp[xval,yval] = output.params['A_3'].value
            pk3_width[xval,yval] = output.params['w_3'].value
            pk3_pos[xval,yval] = output.params['p_3'].value

            pk4_amp[xval,yval] =  output.params['A_4'].value
            pk4_width[xval,yval] = output.params['w_4'].value
            pk4_pos[xval,yval] =  output.params['p_4'].value
            
            fit_quality[xval,yval] = output.chisqr

    # Save the data as a structured array

    newstruct['DCoffset'] = offst
    newstruct['A_1'] = pk1_amp
    newstruct['w_1'] = pk1_width
    newstruct['p_1'] = pk1_pos
    newstruct['A_2'] = pk2_amp
    newstruct['w_2'] = pk2_width
    newstruct['p_2'] = pk2_pos
    newstruct['A_3'] = pk3_amp
    newstruct['w_3'] = pk3_width
    newstruct['p_3'] = pk3_pos
    newstruct['A_4'] = pk4_amp
    newstruct['w_4'] = pk4_width
    newstruct['p_4'] = pk4_pos
    newstruct['A_0'] = water_amp
    newstruct['w_0'] = water_width
    newstruct['p_0'] = water_pos

    return {'':newstruct,'fit_quality':fit_quality} 

def process_cest(scn_to_analyse, xval, yval, pdata_num = 0):
    
    scn = sarpy.Scan(scn_to_analyse)
    
    # Defining parameters
    freq_list = scn.method.CEST_FreqListPPM

    ppm_limit_min = -50
    ppm_limit_max = 50
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
    
    # get the freqs that'll be used for water fit
    water_fit_freqs = [f for f in ppm_filtered if (numpy.abs(f)< 3.)]
    water_fit_freqs_ind = sorted([ppm_filtered.index(c) for c in water_fit_freqs])
    
    # First do the water fit and shift the data so water is at 0  
    waterShift = fit_water_peak(tmp[water_fit_freqs_ind],water_fit_freqs)

    # Interpolating the Y-data so that it gets shifted to the acquired offsets!
    s_shifted_back = scipy.interp(ppm_filtered, ppm_filtered-waterShift, tmp)
    #print('watershift is: {0}'.format(waterShift))

    if ~numpy.isfinite(numpy.sum(s_shifted_back)):
        print('\n \n Look here: {} \n \n',scn_to_analyse, xval,yval)
    
    return ppm_filtered, s_shifted_back


################################################

####### Displaying and Plotting Functions #######    

################################################

def plotCestPeaks(cestParams,x,y):
    ''' This function takes in a full fit and returns a plot of the individual peaks plotted on the flipped axis'''
    import pylab

    freqs = numpy.arange(-20,20,0.1)

    try:

        for i in numpy.arange(0,int(len(cestParams[x,y])/3)-1):
            pylab.plot(freqs,1-sarpy.analysis.cest.h_peak_N(cestParams[x,y],freqs,i))

    except:
        for i in numpy.arange(0,int(len(cestParams.params)/3)):
            pylab.plot(freqs,1-sarpy.analysis.cest.h_peak_N(cestParams.params,freqs,i))

    pylab.axvline(2.2,ymin=0,label='2.2 amine',color='y', alpha=0.4)
    pylab.axvline(3.5,ymin=0,label='3.5 amide',color='r', alpha=0.4)
    pylab.axvline(-3.25,ymin=0,label='-3.25 aliphatic',color='g', alpha=0.4)
    pylab.axvline(-3.0,ymin=0,label='-3.0')    
    #pylab.axvline(1.5,ymin=0.6,label='1.5 OH',color='b', alpha=0.4)


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