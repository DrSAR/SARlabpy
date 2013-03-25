#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fitting Pharmaco-Kinetic models to dynamic contrast enhanced data

Included models: 
    - Tofts
    - Extended Tofts (with vascular term)
    - 2CXM
    
    
Created on Wed Jun  1 10:45:00 2011

@author: Stefan Reinsberg
"""

import numpy 
from numpy import exp, sqrt
from scipy import optimize

import matplotlib.pyplot as plt

def AIF_factory(model='Lyng'):
    '''Arterial Input Function factory
    
    Sets up a function that will correspond to one of several possible 
    (literature) AIFs. Model (default='Lyng') determines the choice of AIF:
    
    Parameters
    ----------
    Model : String
        'Lyng'          ... from Lyng 1998 - mice, same Gd
        
        'Checkley'      ... from Checkley ???
        
        'Pickup'        ... from Pickup 2004 - mice
        
        'Fritz-Hansen'  ... not sure and hence not implemented FIXME

    Returns
    -------
    (f, vfunc) : 2-tupel of functions
        A function that accepts 1D float parameters (time).
        and a function that accepts a vector of floats
        
        Will raise a TypeError exception on unkown AIFs
        
    Examples
    --------
    >>> a=PKfit.AIF_factory(model='Lyng')
    choosing Lyng AIF    
    >>> a[0](0)
    6.5
    >>> a[1](array([0,1,2])
    array([ 6.5       ,  6.0893048 ,  5.70760869])
    '''
    
    if model == 'Lyng':
        print('choosing Lyng AIF')
        func = lambda t: 5.8*exp((-4.4/60.)*t) + 0.7*exp((-0.05/60.)*t) 
    elif model == 'Checkley':
        print('choosing Checkley AIF')
        func = lambda t: .3*11.95*exp((-.0195)*t) + .3*4.67*exp((-0.0009)*t)
    elif model == 'Pickup':
        print('choosing Pickup AIF')
        func = lambda t: 0.19*exp(-0.069*t) + 0.1*exp(-0.00105*t) 
    elif model == 'Fritz-Hansen': # FIXME
        print('choosing Fritz-Hansen AIF')
        raise TypeError('Parameters not yet coded - pls look up and fix') 
        func = lambda t: exp(-t/10) + exp(-t/100)
    else:
        raise TypeError('Unknown AIF choice')
        
    vfunc = numpy.vectorize(func)
    
    return (func, vfunc)
    
    
def params_2CXM_to(params):
    '''mapping (Em, Kp, Km, Fp) -> (Vp, Ve, Fe, Fp)'''
    (Em, Kp, Km, Fp) = params
    Tb = 1./(Kp-Em*(Kp-Km))
    Te = 1./(Tb*Kp*Km)
    Tp = 1./(Kp+Km-1./Te)
    
    Vp = Fp*Tb
    Fe = Fp*(Tb/Tp-1)
    Ve = Fe * Te
    
    return (Vp, Ve, Fe, Fp)
    
def params_2CXM_fro(params):
    '''mapping (Vp, Ve, Fe, Fp) -> (Em, Kp, Km, Fp)'''
    (Vp, Ve, Fe, Fp) = params

    Tp = Vp / (Fe+Fp)
    Te = Ve / Fe
    Tb = Vp / Fp
    
    Kp = .5*(1./Tp + 1./Te + 
        sqrt((1./Tp+1./Te)**2 - 4./(Te*Tb)))
    Km = .5*(1./Tp + 1./Te - 
        sqrt((1./Tp+1./Te)**2 - 4./(Te*Tb)))
    Em = (Kp-1./Tb) / (Kp-Km)    
            
    return (Em, Kp, Km, Fp)
    
def TwoCXM_factory(time_axis, t0, AIFmodel=None):
    '''2CXM 
    
    Sourbron provides an analytical solution in 
    `Quantification of cerebral blood flow, cerebral blood volume, and 
    blood-brain-barrier leakage with DCE-MRI
    Sourbron et al, Magnetic Resonance in Medicine
    Volume 62, Issue 1, pp 205
    <http://onlinelibrary.wiley.com/doi/10.1002/mrm.22005/abstract>`_

    Parameters
    ----------
    time_axis: float array
        time points for which to calculate concentrations 
        
    t0: int
        time of contrast arrival - AIF needs to be shifted by that much
        
    AIFmodel: string(None)
        used by AIF_factory
        
    params: 4-element float array 
        [Em, Kp, Km, Fp] 
        
        !! not supplied to this factory function but needed in calls to the 
        returned function handle - just thought we should mention this here
        
    Returns
    -------
    conc: function handle 
        2XCM function evaluated at the locations of time_axis, t, with params as
        parameters (see comment above):
            conc(params, t) = Fp * numpy.convolve(res(t), AIF(t), mode='valid')
    '''

    # set up AIF vector
    AIFfunc, AIFvfunc = AIF_factory(model=AIFmodel)
    # prefix AIF with time_axis.shape[0]+t0 zeros 
    AIF = numpy.hstack((numpy.zeros(time_axis.shape[0]+t0-1), 
                       (AIFvfunc(time_axis))[0:-t0]))

    #set the four parameters governing the tissue residual 
    #(Em, Kp, Km, Fp) = params
    #res = exp(-time_axis*Kp) + Em*(exp(-time_axis*Km) - exp(-time_axis*Kp))
        
    # perform convolution to get concentration time curve
    # mode='full' will return a convolution whose length is the sum of length
    # of res and AIF, we would need to discard points which sounds costly
    # mode='valid' appears to give the correct results in terms of length if 
    # the length of the input vectors is tuned appropriately: We have to
    # prepend N zeros to the AIF, where N is the number of points in residue
    # function res
    # conc = Fp * numpy.convolve(res, AIF, mode='valid')

    TwoCXM = lambda (Em, Kp, Km, Fp), time: Fp * numpy.convolve(
        exp(-time_axis*Kp)+Em*(exp(-time_axis*Km)-exp(-time_axis*Kp)),
        AIF,
        mode='valid')
                                    
    
    return TwoCXM
    
def PKfitfunction(params):
    '''cost function (sum of least squares?) in the process of fitting one 
    of any PK model functions
    
    Parameters
    ----------
    params: tupel of floats
        will be passed to the function calculating the model concentration curve
        
    Returns
    -------
    res: scalar
        sum of least squares or some such useful measure of goodness of fit
    '''
    
    return 0    
    
def ToftsModel(t):
    '''Toft's model
    '''
    pass


def PKfit():
    '''wrapper function to read data, choose model and AIF, fit model to data,
    and show results
    
    PS: consider the `Cookbook/FittingData <http://www.scipy.org/Cookbook/FittingData>`_
    for help
    '''
    
    #read data
    #C=conc(time_axis)
    Npts = 400
    time=numpy.linspace(0, 30., num=Npts)
    DataGenerator = TwoCXM_factory(time, 10, AIFmodel='Lyng')
    pinput = numpy.array([0.5, .15, .2, .01])
    conc = DataGenerator(pinput, time) + numpy.random.randn(Npts)/50
    #plot data
    plt.plot(time, conc, "ro")   
    
    #setup fit by choosing AIF and model
    PKfitfunc = TwoCXM_factory(time, 10, AIFmodel='Lyng')
    errfunc = lambda p, t, c: PKfitfunc(p, t) - c # distance to data
    
    #perform fit
    p0=pinput+pinput*numpy.random.randn(4)/10
    p1, success = optimize.leastsq(errfunc, p0[:], args=(time, conc))
    #show results
    print(p1, pinput)
    print((p1-pinput)/pinput)
    print(params_2CXM_to(p1))
    plt.plot(time, PKfitfunc(p1, time), "b-")
    plt.show()
    