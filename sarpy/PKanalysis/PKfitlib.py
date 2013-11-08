#!/usr/bin/env python
''' Functions around modelling of contrast agent uptiake in tissue
mostly inspired by the 2-Comparment eXchange Model (2CXM) paper by 
Sourbron and Buckley paper.
'''

from numpy import exp, sqrt, log, pi
import numpy
import scipy.optimize as optimize

def heaviside(time, tau):
    """Heaviside function for use in AIF model."""
    out = numpy.zeros_like(time)
    out[time>tau]=1
    return out

def AIF_model(parms, time):
    ''' Create AIF with linear upslope and single exponential decay.'''
#    tau = vol/rate
    alpha, beta, vol, c_f = parms[0], parms[1], parms[2], parms[3]
    tau = numpy.sqrt(2*vol/alpha)
    up = ( 1 - heaviside(time,tau) ) * alpha * time
    down = heaviside(time,tau) *\
        ( (alpha*(tau)-c_f) *\
        numpy.exp(-beta*(time-tau) ) + c_f  )
    return up+down
    
def create_model_AIF(parms, time_axis, 
                     modify_rate = False,
                     shift_h = False, 
                     shift_v = False,
                     shift_arriv = False,
                     scale = False,
                     epsilon = None):
    '''Create AIF for different rates (alphas), and fixed parameters
    from fit of data with AIF_model function. 
    
    coopyright Tammo Rukat, 2013'''
        
    if modify_rate == True:
        parms = numpy.hstack([parms[0]*(1+epsilon),parms[1:5]])
        
    alpha, V1, V2, tf, cf = parms[0], parms[1], parms[2], parms[3], parms[4]    
    Vtotal = V1+V2
    
    tau = sqrt( (2*V1)/alpha )

    if shift_h == True:
        alpha_old = alpha
        tau_old = sqrt( (2*V1)/alpha )
        
        const = alpha_old * tau_old
        tau_new = tau_old * (1 + epsilon)
        alpha_new = const/tau_new        
        V1_new = (const*tau_new)/2
        V2_new = Vtotal - V1_new
        
        V1 = V1_new
        V2 = V2_new
        alpha = alpha_new
        tau = tau_new
    
    if shift_v == True:
        alpha_old = alpha
        const = (2*V1)/alpha
        tau_old = sqrt(const) 
        
#        epsilon *= alpha_old * tau_old - cf
        alpha_new = alpha_old*(1 + epsilon) - (epsilon*cf/tau_old)
        V1_new = (const*alpha_new)/2
        V2_new = Vtotal - V1_new
        
        V1 = V1_new
        V2 = V2_new
        alpha = alpha_new
        tau = tau_old
        
    beta = (alpha*tau - cf) / (V2 - cf*(tf-tau))
    AIF = AIF_model([alpha, beta, V1, cf], time_axis)
    if shift_arriv == True:
        AIF = AIF_model([alpha, beta, V1, cf], time_axis-epsilon)
        AIF[AIF<0] = 0
    if scale == True:
        AIF *= 1 + epsilon
    return AIF

def aif(t, t0=0, aifchoice=None, zeropad=False):
    '''arterial input function by Lyng
    
    Parameters
    ----------
    t: scalar or vector
        time in min
        
    t0: scalar (default = 0)
        shift (+ to the left) of start, units of t

    aifchoice : String
        'Lyng'          ... from Lyng 1998 - mice, same Gd
        
        'Checkley'      ... from Checkley ???
        
        'Pickup'        ... from Pickup 2004 - mice
        
        'Fritz-Hansen'  ... not sure and hence not implemented FIXME

    zeropad: boolean (default = False)
        return zero-padded AIF to the left with len(t)-1 zeros
        
    Results
    -------
    returns a tupel of the AIF and (if requested) a zeropadded AIF
    '''
    
    AIF = numpy.zeros(len(t))
    if aifchoice == 'Lyng':
        AIF[t>=t0] = 5.8*exp((-4.4)*(t[t>=t0]-t0)) + 0.7*exp((-0.05)*(t[t>=t0]-t0))
    elif aifchoice == 'Moroz':
        # load parameters to simulate aif
        aif_parms = [5.18134715e-02, 9.72041215e+00, 
                     4.59268952e+01, 3.59900000e+02, 
                     1.20795007e-01]
        # create aif
        AIF[t>=t0] = create_model_AIF(aif_parms, t[t>=t0]-t0)
    elif aifchoice == 'Checkley':
        AIF[t>=t0] = (.3*11.95*exp((-.0195)*(t[t>=t0]-t0)) +
                      .3*4.67*exp((-0.0009)*(t[t>=t0]-t0)))
    elif aifchoice == 'Pickup':
        AIF[t>=t0] = 0.19*exp(-0.069*(t[t>=t0]-t0)) + 0.1*exp(-0.00105*(t[t>=t0]-t0)) 
    elif aifchoice == 'HPG':
        AIF[t>=t0] = (5.8*exp((-4.4)*(t[t>=t0]-t0)) + 0.7*exp((-0.05)*(t[t>=t0]-t0)))/5
    elif aifchoice == 'Parker': # model from DOI: 10.1002/mrm.21066
        # this is a population average (23 Px in 113 studies), 
        # gadodiamide .1 mmol/kg body weight of Omniscan 0.5mmol/ml
        # automatically in descending aorta or iliac arteries
        A1 = .809 # mmol.min
        A2 = .330 # mmol.min
        T1 = .17046 # min
        T2 = .365   # min
        s1 = .0563  # min
        s2 = .132   # min
        alpha = 1.050 # mmol 
        beta  = .1685 # min^-1
        s     = 38.078# min^-1
        tau   = .483  # min

        AIF[t>=t0] =(
            A1/s1/sqrt(2*pi)*exp(-(t-T1)**2/(2*s1**2))  +
            A2/s2/sqrt(2*pi)*exp(-(t-T2)**2/(2*s2**2)) + 
            alpha * exp(-beta*t) / (1+exp(-s*(t-tau-t0))))        
    else:
        raise TypeError('unidentified AIF choice')
        
    if zeropad:
        return numpy.concatenate((numpy.zeros(len(t)-1), AIF))
    else:
        return AIF

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
    >>> a=AIF_factory(model='Lyng')
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

def residual_Tofts(t, params):
    '''simple Tofts tissue response function - see equation [37] and [38]i
 
    params = (Ktrans, ve)         
    '''
    
    (Ktrans, ve) = params
    
    res = Ktrans*exp(-t*Ktrans/ve)
    
    return res

def conc_Tofts(t, modelparams, ca, dt):
    '''concentration time curve (see eqn [1]) from simple Tofts model
    
    Parameters
    ----------
    t: time (vector)
    
    modelparams: (Ktrans, ve) (2-tupel)

    ca: arterial input function (vector), should have length of 2*len(t)-1 
    to allow 'valid' convolution with impulse response function

    Results
    -------
    concentration time curve for time points t
    '''
    res = residual_Tofts(t,modelparams) # see eqn [15]
    
    c = numpy.convolve(ca,res, mode='valid')*dt
    
    return c


def conc_XTofts(t, modelparams, ca, dt):
    '''concentration time curve (see eqn [2]) from extended Tofts model. The 
    ad-hoc vascular term is included in addition to the convolution with the 
    Tofts residue, tissue response function.
    
    Parameters
    ----------
    t: time (vector)
    
    modelparams: (Ktrans, ve, vp) (3-tupel)

    ca: arterial input function (vector), should have length of 2*len(t)-1 
    to allow 'valid' convolution with impulse response function

    Results
    -------
    concentration time curve for time points t
    '''    
    (Ktrans, ve, vp) = modelparams
        
    c = vp * ca[-len(t):] + conc_Tofts(t, (Ktrans, ve), ca, dt)
    return c


def residual_2CXM(t, params):
    '''tissue response function
    params = (Fp, Fm, Kp, Km)'''
    
    (Fp, Fm, Kp, Km) = params
    
    res = Fp*exp(-t*Kp) + Fm*exp(-t*Km)
    return res

def conc_2CXM(t, modelparams, ca, dt):
    '''concentration time curve (see eqn [8]) from 2 parameter exchange model
    
    Parameters
    ----------
    t: time (vector)
    
    modelparams: (PS, Fpl, ve, vp) (4-tupel)

    ca: arterial input function (vector), should have length of 2*len(t)-1 
    to allow 'valid' convolution with impulse response function

    Results
    -------
    concentration time curve for time points t
    '''
    
        #set up tissue response function
    params = paramconv_modeltosoln(modelparams)
    res = residual_2CXM(t,params) # see eqn [15]
    
    
    c = numpy.convolve(ca,res, mode='valid')*dt
    
    return c
    
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
    
def paramconv_modeltosoln(params):
    ''' convert models used to set up the model (PS, Fpl, ve, vp) into 
    parameters more useful in the solution (Fp, Fm, Kp, Km)
    
    function paramconv_solntomodel() is the inversion'''
   
    (PS, Fpl, ve, vp) = params
    
    e = ve / (vp + ve)  # extravascular fraction of the extracellular volume
    
    E = PS / (PS + Fpl)  # extraction fraction
    
    root = sqrt(1-4*E*e*(1-E)*(1-e)/((E-E*e+e)**2))
    taup = (E - E*e +e)/(2*E) * (1 + root)
    taum = (E - E*e +e)/(2*E) * (1 - root)
    
    Kp = Fpl / ((vp + ve)*taum)
    Km = Fpl / ((vp + ve)*taup)
    
    
    Fp =  Fpl * (taup - 1)/(taup - taum)
    Fm = -Fpl * (taum - 1)/(taup + taum)
    
    return (Fp, Fm, Kp, Km)
    
    
def paramconv_solntomodel(params):
    ''' convert models used to set up the model (PS, Fpl, ve, vp) into 
    parameters more useful in the solution (Fp, Fm, Kp, Km)
    
    this is the inverse of paramconf_modeltosoln()'''
    
    (Fp, Fm, Kp, Km) = params
    
    taum = (Fp/(Fm+1)) / (Fp/Fm + Kp/Km)
    taup = taum * Kp/Km
    
    Fpl = Fp * (taup- taum) / (taup - 1)
    
    E = (taup + taum -1 - taup*taum) / ((taup+taum)*(taup+taum-1)-taup*taum)
    
    e = (taup + taum -1 - taup*taum)  / (taup+taum-1)
    
    ve = e * Fpl / (Kp * taum)
    
    vp =  Fpl / (Kp * taum) - ve
    
    PS = Fpl * E / (1-E)

    return (PS, Fpl, ve, vp)
    
def bounds_penalty(parms, bounds):
    '''implementation of a check to see whether a list of parameters is within the 
    given bounds. A bound is a 2D array
    
    Interesting observation that the awkward loop is faster (by a lot!) than
    if numpy.any(parms<bounds[0]) or numpy.any(parms>bounds[1]):
        return 1e9
    else:
        return 0.0
    The performance is even worse if you start messing with numpy.sign
    '''
    for idx, p in enumerate(parms):
        if p > bounds[idx][1] or p < bounds[idx][0]:
            return 1e30
    return 0.0

def fit_generic(t, data, model, p0, *pargs, **kwargs):
    '''fit a generic function
    
    Parameters
    ----------
    t: time (vector)
    data: concentration time curve (vector)
    model: function handle
    p0: starting guess
    
    *pargs
    ca: AIF (vector) note that len(ca)=2*len(t)-1 for 'valid' convolution
    e.g.    ca = aif(t,aifchoice='Lyng')
            ca = numpy.concatenate((numpy.zeros(len(t)-1), ca))
    **kwargs
    fit_subset: vector (or slice) to define the subset of model results to
                be evaluated (default: slice(None))
    Results
    -------
    modelparameters: tupel depending on model
    success:         as handed back from optimize.leastsq()
    fit:             best fit
    Tip: 
    ----
    Always consult docs: http://wiki.scipy.org/Cookbook/FittingData
    '''

    bounds = kwargs.get('bounds', None)
    fit_subset = kwargs.get('fit_subset', slice(None)) 
    if bounds is None:
        errfunc = lambda p, x, pargs, y: model(x, p, *pargs)[fit_subset] - y
    else:
        errfunc = lambda p, x, pargs, y: numpy.r_[
                                            model(x, p, *pargs)[fit_subset] - y,
                                            bounds_penalty(p, bounds)]
                                                 
    p1, success = optimize.leastsq(errfunc, p0[:], args=(t, pargs, data))

    fit = model(t, p1, *pargs)     
    ss_err = sum(errfunc(p1, t, pargs, data)**2)
    ss_tot=((data-data.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)

    return (p1, success, fit, ss_err, rsquared)

def fit_generic_array(t, data, model, p0, *pargs, **kwargs):
    '''wrapper for fit_generic but data is a 4D array instead of a 
    one-dimensional vector
    :mask: describes a 3D mask
    :data: is a 4D array with the first three dimensions matched by
           the data
    '''
    mask = kwargs.pop('mask', None)
    result = numpy.ndarray(mask.shape+p0.shape, dtype=float) + numpy.NaN
    fit = numpy.ndarray(mask.shape+t.shape, dtype=float) + numpy.NaN
    ss = numpy.ndarray(mask.shape, dtype=float) + numpy.NaN
    rsq= numpy.ndarray(mask.shape, dtype=float) + numpy.NaN
    success = numpy.ndarray(mask.shape, dtype=int) + numpy.NaN

    for idx, val in numpy.ndenumerate(numpy.isnan(mask)):
        if not val:
            (res, scs, ft, sos, rs) = fit_generic(t, 
                                             data[idx],
                                             model, p0, *pargs, **kwargs)
            result[idx] = res
            fit[idx] = ft
            success[idx] = scs
            ss[idx] = sos
            rsq[idx] = rs
                                             
    return {'result':result, 'fit':fit, 'ss':ss, 'rsq': rsq,
            'success':success}     
            
            
def fit_2CXM(t, ca, data, p0):
    '''fit the 2CXM with its four parameters to some data
    
    Parameters
    ----------
    t: time (vector)
    
    ca: AIF (vector) note that len(ca)=2*len(t)-1 for 'valid' convolution
        e.g.     
            ca = aif(t,aifchoice='Lyng')
            ca = numpy.concatenate((numpy.zeros(len(t)-1), ca))
 
    data: concentration time curve (vector)
    
    p0: starting guess (4-tupel)
    
    Results
    -------
    modelparameters: (ve, vp, PS, Fpl) (4-tupel'''


    # distance to the target function    
    errfunc = lambda p, x, aif, y: conc_2CXM(x, p, aif) - y
    
    p1, success = optimize.leastsq(errfunc, p0[:], args=(t, ca, data))
            
    return (p1, success)
    
    
def AIC_from_SSE(SSE=None, k=None, N=None):
    '''
    `AIC is formally defined <http://en.wikipedia.org/wiki/Akaike_information_criterion>`_
    as 2k - 2ln(L)
    
    Assumption of gaussian noise around measured values leads to minimizing SSE 
    as the best strategy to maximize the likelihood of obtaining the optimal 
    parameters in nonlinear regression.
    
    Hence, AIC becomes:
        AIC = N ln (SSE/N) + 2k
        N ... number of data points
        SSE ... sum square of errors
        k ... number of fit parameters + 1 (since you're also estimating SSE)
    
    Note: this seems to imply that L = (SSE/N)^(-N/2) for Gaussian error        
        
    The 2nd-order corrected AIC (for finite sample sizes)
        AICc= AIC + 2K(K+1) / (N-K-1)
    Note: You need at least 3 more data points than parameters.'''

    AIC = N * log(SSE/N) + 2*k
    AICc = AIC + 2*k*(k+1)/(N-k-1)  # AICc = AIC for large N
    return AICc

def rel_prob_from_AIC(AICs, axis=-1):
    '''
    Relative probabilities for models whose AIC values are given

    :param ndarray AICs: 
        array that contains the calculated AIC for every model under
        consideration
    :param int axis: 
        axis along which the model AIC are organized in the AICs array 
        default: -1 (last value)

    The model with a minimum AIC will have a rel likelihood of 1. The other
    models will have a likelihood indicating the probability (relative to the
    minimum one) to minimize the information loss due to using a candidate model
    to represent the "true" model.
    '''
    AICmin = numpy.min(AICs, axis=axis, keepdims=True)
    rel_prob = numpy.exp((AICmin-AICs)/2)
    return rel_prob

def hierarchical_fit_2CXM(t, ca, data, p0):
    ''' Follow the suggestion by Sourbron & Buckley wrt fitting 
    models of increasing complexity'''
    
    # our choice of nested models up to the full 2CXM
    models = {
    'Tofts'  : (conc_Tofts, [0.01, 0.01]),
    'X-Tofts': (conc_XTofts, [0.01, 0.01]),
    '2CXM'   : (conc_2CXM, [0.01, 0.01])
    }
    
    aic = {}
    fitres = {}
    SSE_arr = numpy.empty(0, dtype=float)
    #k_arr = numpy.empty(0, dtype=int)
    
    for modellabel, model in models.iteritems():
        fitres = fit_generic(t, data, model[0], model[1], ca)
        # find Akaike's Information Criterion
        SSE_arr = [SSE_arr, 1]
        #k_arr = len(model[1])+1
        
        aic[modellabel] = AIC_from_SSE(fitres) 

    return (fitres, aic)
    

def fit_double_2CXM(t, ca1, ca2, data1, data2, p0):
    '''fit the 2CXM with its four parameters to some data that has been 
    acquired with two(!) contrast agents
    
    Parameters
    ----------
    t: time (vector)
    
    ca1, ca2: AIF for Gd-DTPA and HPG-Gd
    data1: concentration measured (Gd-DTPA)
    
    data2: concentration measured (HPG-Gd)
    
    p0: starting guess (4-tupel)
    
    Results
    -------
    modelparameters: (ve, vp, PS, Fpl) (4-tupel)'''


    # distance to the target function    
    errfunc = lambda p, x, aif1, aif2, y1, y2: numpy.r_[
                    conc_2CXM(x, p[0:4], aif1) - y1,
                    conc_2CXM(x, numpy.r_[p[0:2],p[4],p[3]], aif2) - y2
                    ]
        
    p1, success = optimize.leastsq(errfunc, p0[:], args=(t, ca1, ca2, data1, data2))
            
    return (p1, success)

    
