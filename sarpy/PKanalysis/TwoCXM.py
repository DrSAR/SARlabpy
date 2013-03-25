#!/usr/bin/env python
''' Function around the 2-Comparment eXchange Model (2CXM). Mostly inspired by 
Sourbron and Buckley paper.
'''

from numpy import exp, sqrt, log, pi
import numpy
import scipy.optimize as optimize
import matplotlib.pyplot as plt

def aif(t, t0=0, aifchoice=None, zeropad=False):
    '''arterial input function by Lyng
    
    Parameters
    ----------
    t: scalar or vector
        time in min
        
    t0: scalar (default = 0)
        shift (+ to the left) of start, units of t
    
    zeropad: boolean (default = False)
        return zero-padded AIF to the left with len(t)-1 zeros
        
    Results
    -------
    returns a tupel of the AIF and (if requested) a zeropadded AIF
    '''
    
    AIF = numpy.zeros(len(t))
    if aifchoice == 'Lyng':
        AIF[t>=t0] = 5.8*exp((-4.4)*(t[t>=t0]-t0)) + 0.7*exp((-0.05)*(t[t>=t0]-t0))
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


'''
vp+ve = 30/100
Fp = 50ml / min / 100ml
PS = 25ml / min / 100ml
'''

def residual_Tofts(t, params):
    '''tissue response function - see equation [37] and [38]'''
    
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
    
    modelparams: (vp, Ktrans, ve) (3-tupel)

    ca: arterial input function (vector), should have length of 2*len(t)-1 
    to allow 'valid' convolution with impulse response function

    Results
    -------
    concentration time curve for time points t
    '''    
    (vp, Ktrans, ve) = modelparams
        
    c = vp * ca[-len(t):] + conc_Tofts(t, (Ktrans, ve), ca, dt)
    return c


def residual_2CXM(t, params):
    '''tissue response function'''
    
    (Fp, Fm, Kp, Km) = params
    
    res = Fp*exp(-t*Kp) + Fm*exp(-t*Km)
    return res

def conc_2CXM(t, modelparams, ca, dt):
    '''concentration time curve (see eqn [8]) from 2 parameter exchange model
    
    Parameters
    ----------
    t: time (vector)
    
    modelparams: (ve, vp, PS, Fpl) (4-tupel)

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
    
def paramconv_modeltosoln(params):
    ''' convert models used to set up the model (ve, vp, PS, Fpl) into 
    parameters more useful in the solution (Fp, Fm, Kp, Km)
    
    function paramconv_solntomodel() is the inversion'''
    
    (ve, vp, PS, Fpl) = params
    
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
    ''' convert models used to set up the model (ve, vp, PS, Fpl) into 
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

    return (ve, vp, PS, Fpl)

def fit_generic(t, data, model, p0,  *pargs):
    '''fit a generic function
    
    Parameters
    ----------
    t: time (vector)

    data: concentration time curve (vector)
    
    model: function handle
    
    p0: starting guess
    
    *pargs
    ca: AIF (vector) note that len(ca)=2*len(t)-1 for 'valid' convolution
    e.g.     
            ca = aif(t,aifchoice='Lyng')
            
            ca = numpy.concatenate((numpy.zeros(len(t)-1), ca))
 
    
    Results
    -------
    modelparameters: tupel depending on model
    
    success:         as handed back from optimize.leastsq()
    
    fit:             best fit'''

    # distance to the target function    
    errfunc = lambda p, x, pargs, y: model(x, p, *pargs) - y
    
    p1, success = optimize.leastsq(errfunc, p0[:], args=(t, pargs, data))

    fit = model(t, p1, *pargs)     
        
    return (p1, success, fit)

    
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
    Note: You need at least 3 more data points that parameters.'''

    AIC = N * log(SSE/N) + 2*k
    AICc = AIC + 2*k*(k+1)/(N-k-1)  # AICc = AIC for large N
    return AICc

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
    k_arr = numpy.empty(0, dtype=int)
    
    for modellabel, model in models.iteritems():
        fitres = fit_generic(t, data, model[0], model[1], ca)
        # find Akaike's Information Criterion
        SSE_arr = [SSE_arr, 1]
        k_arr = len(model[1])+1
        
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

    
def parameter_stability(N=100):

    #create some data
    t=numpy.arange(400.)/60
    
    # set up AIFs
    aif_Gd = aif(t,aifchoice='Lyng')
    ca_Gd = numpy.concatenate((numpy.zeros(len(t)-1), aif_Gd))

    aif_HPG = aif(t,aifchoice='HPG')
    ca_HPG = numpy.concatenate((numpy.zeros(len(t)-1), aif_HPG)) 

    # create and fit intermediate case
    print('Gd fit with intermediate vascularization')
    p0 = [5./100, 5./100, 25./100, 50./100]
    print(p0)
    param1_i=[]
    param2_i=[]
    for i in numpy.arange(N):
        c_ivasc = conc_2CXM(t, p0, ca_Gd) + numpy.random.randn(len(t))/4
        p1, success = fit_2CXM(t, ca_Gd, c_ivasc, p0)
        param1_i.append(p1[1])
        param2_i.append(p1[2])

    
    # create and fit high vasc case
    print('Gd fit with high vascularization')
    p0 = [1./100, 9./100, 25./100, 50./100]
    print(p0)
    param1_h=[]
    param2_h=[]
    for i in numpy.arange(N):
        c_hvasc = conc_2CXM(t, p0, ca_Gd) + numpy.random.randn(len(t))/4
        p1, success = fit_2CXM(t, ca_Gd, c_hvasc, p0)
        #print(p1[paramnr])
        param1_h.append(p1[1])
        param2_h.append(p1[2])
    
    #double fit
    print('double fit')
    p0.append(0.01/100)
    print(p0)
    param1_d=[]
    param2_d=[]
    for i in numpy.arange(N):
        c_hvasc = conc_2CXM(t, p0[0:4], ca_Gd) + numpy.random.randn(len(t))/4
        cHPG_hvasc = conc_2CXM(t, numpy.r_[p0[0:2],p0[4],p0[3]], ca_HPG
                          ) + numpy.random.randn(len(t))/4
        p1, success = fit_double_2CXM(t, ca_Gd, ca_HPG, 
                            c_hvasc, cHPG_hvasc, p0)
        #print(p1[paramnr])
        param1_d.append(p1[1])            
        param2_d.append(p1[2])            

    plt.subplot(222)
    plt.hist(param1_h,color='b')
    plt.hist(param1_i,color='r')
    plt.hist(param1_d,color='g')
    plt.plot(sorted(param1_h),numpy.arange(len(param1_h))/2,'b-')
    plt.plot(sorted(param1_i),numpy.arange(len(param1_i))/2,'r-')
    plt.plot(sorted(param1_d),numpy.arange(len(param1_d))/2,'g-')

    plt.title(r'Histogram of $v_p$')
    plt.xlabel(r'plasma volume, $v_p$')
    
    plt.subplot(224)
    plt.hist(param2_h,color='b')
    plt.hist(param2_d,color='g')
    plt.hist(param2_i,color='r')
    plt.plot(sorted(param2_h),numpy.arange(len(param2_h)),'b-')
    plt.plot(sorted(param2_d),numpy.arange(len(param2_d)),'g-')
    plt.plot(sorted(param2_i),numpy.arange(len(param2_i)),'r-')
    plt.title(r'Histogram of $PS$')
    plt.xlabel(r'permeability surface product, $PS$')
    plt.subplots_adjust(hspace=0.3)
        
    print(numpy.mean(param1_i), numpy.std(param1_i))
    print(numpy.mean(param1_h), numpy.std(param1_h))
    print(numpy.mean(param1_d), numpy.std(param1_d))
    print(numpy.mean(param2_i), numpy.std(param2_i))
    print(numpy.mean(param2_h), numpy.std(param2_h))
    print(numpy.mean(param2_d), numpy.std(param2_d))

    #    return (param_i, param_h, param_d)
    
    
def test_fit():

    #create some data
    t=numpy.arange(400.)/60
    p0 = [5./100, 5./100, 25./100, 50./100]
    
    # set up AIFs
    aif_Gd = aif(t,aifchoice='Lyng')
    ca_Gd = numpy.concatenate((numpy.zeros(len(t)-1), aif_Gd))

    aif_HPG = aif(t,aifchoice='HPG')
    ca_HPG = numpy.concatenate((numpy.zeros(len(t)-1), aif_HPG)) 

    # create and fit intermediate case
    c_ivasc = conc_2CXM(t, p0, ca_Gd) + numpy.random.randn(len(t))/4
    
    print('Gd fit with intermediate vascularization')
    print(p0)
    p1, success = fit_2CXM(t, ca_Gd, c_ivasc, p0)
    print(p1)
    print((p1-p0)/p0*100)
    fit_ivasc = conc_2CXM(t, p1, ca_Gd)    

    # create and fit highly vascular case    
    p0 = [1./100, 9./100, 25./100, 50./100, 0.01/100]
    c_hvasc = conc_2CXM(t, p0[0:4], ca_Gd) + numpy.random.randn(len(t))/4
    print('Gd fit with high vascularization')
    print(p0)
    p1, success = fit_2CXM(t, ca_Gd, c_hvasc, p0[0:4])
    print(p1)
    print((p1-p0[0:4])/p0[0:4]*100)
    fit_hvasc = conc_2CXM(t, p1[0:4], ca_Gd)
    
    # create HPG case (PS=0)
    cHPG_hvasc = conc_2CXM(t, numpy.r_[p0[0:2],p0[4],p0[3]], ca_HPG
                     ) + numpy.random.randn(len(t))/4
    print('Gd fit with high vascularization and PS=0')
    print(p0)
    p1, success = fit_2CXM(t, ca_HPG, cHPG_hvasc, numpy.r_[p0[0:2],p0[4],p0[3]])
    print(p1)
    print((p1-p0[0:4])/p0[0:4]*100)
    fitHPG_hvasc = conc_2CXM(t, p1, ca_HPG)
    
    #double fit
    print('double fit')
    print(p0)
    p1, success = fit_double_2CXM(t, ca_Gd, ca_HPG, c_hvasc, cHPG_hvasc, p0)
    print(p1)
    print((p1-p0)/p0*100)
    fitd_hvasc = conc_2CXM(t, p1[0:4], ca_Gd)
    fitHPGd_hvasc = conc_2CXM(t, numpy.r_[p0[0:2],p0[4],p0[3]], ca_HPG)

    #plt.subplot(121)
    plt.title('simulated data: 2CXM')
    plt.plot(t, c_ivasc,'ro')
    plt.plot(t, fit_ivasc, 'r-',label=r'intermediate vasc, $v_p$=5 ml/100ml')
    plt.plot(t, c_hvasc,'bo')
    plt.plot(t, fit_hvasc, 'b-',label='highly vasc, $v_p$=9 ml/100ml')
    '''plt.plot(t, cHPG_hvasc,'gx', label='high-molecular weight agent')
    plt.plot(t, fitHPG_hvasc, 'g-',label=r'highly vasc and slow exchange, $PS=0$')'''
    plt.xlabel('time / min')
    plt.ylabel('concentration / a.u.')
    plt.legend()
    
    '''plt.subplot(122)
    plt.title('simulated data - simultaneous fit')
    plt.plot(t, c_hvasc,'bo',label='highly vasc - Gd-DTPA')
    plt.plot(t, fitd_hvasc, 'b-')
    plt.plot(t, cHPG_hvasc,'go',label=r'highly vasc and slow exchange')
    plt.plot(t, fitHPGd_hvasc, 'g-')
    plt.legend()    
    '''
    plt.show()

    
    
def dryrun():
    t=numpy.arange(400.)/60
    
    #prefix zeros so that we can run the convolution in the valid mode 
    #requiring total overlap of the two input vectors
    aif_Gd = aif(t,aifchoice='Lyng')
    ca_Gd = aif(t,aifchoice='Lyng', zeropad=True)
    
    aif_HPG = aif(t,aifchoice='HPG')
    ca_HPG = aif(t,aifchoice='HPG', zeropad=True)
    
    # intermediately vascualarized
    ve = 20./100
    vp = 10./100
    PS = 25./100
    Fpl = 50./100
    print(paramconv_modeltosoln((ve, vp, PS, Fpl)))
    c_ivasc = conc_2CXM(t, (ve, vp, PS, Fpl), ca_Gd)
        
    # intermediately vascularized and slow tracer exchange regime
    ve =  20./100
    vp =  10./100
    PS =   0.01/100
    Fpl = 50./100
    cHPG_ivasc = conc_2CXM(t, (ve, vp, PS, Fpl), ca_HPG)
        
    # weakly vascularized
    ve = 29./100
    vp = 1./100
    PS = 25./100
    Fpl = 50./100
    c_wvasc = conc_2CXM(t, (ve, vp, PS, Fpl), ca_Gd)
        
    # weakly vascularized and slow tracer exchange regime
    ve =  29./100
    vp =   1./100
    PS =   0.01/100
    Fpl = 50./100
    cHPG_wvasc = conc_2CXM(t, (ve, vp, PS, Fpl), ca_HPG)

    # highly vascularized
    ve = 1./100
    vp = 29./100
    PS = 25./100
    Fpl = 50./100
    c_hvasc = conc_2CXM(t, (ve, vp, PS, Fpl), ca_Gd)

    # highly vascularized and low tracer exchange
    ve = 1./100
    vp = 29./100
    PS = 0.01/100
    Fpl = 50./100
    cHPG_hvasc = conc_2CXM(t, (ve, vp, PS, Fpl), ca_HPG)
        
    fig=plt.figure()
    plt.subplot(121)
    plt.plot(t, c_wvasc, label=r"weakly vasc $v_p\to 0$")
    plt.plot(t, c_ivasc, label='intermediately vasc')
    plt.plot(t, c_hvasc, label=r'highly vasc $v_e\to 0$')
    plt.plot(t, aif_Gd, 'r--')
    plt.xlabel('time / min')
    plt.title(r"Typical Tracer Exchange: $PS>0$")
    plt.ylabel('tracer concentration / a.u.')
    plt.legend()
    plt.subplot(122)
    plt.plot(t, cHPG_wvasc, label=r"weakly vasc $v_p\to 0$")
    plt.plot(t, cHPG_ivasc, label='intermediately vasc')
    plt.plot(t, cHPG_hvasc, label=r"highly vasc $v_e \to 0$")
    plt.xlabel('time / min')
    plt.ylabel('tracer concentration / a.u.')
    plt.plot(t, aif_HPG, 'r--')
    plt.title(r"Slow Tracer Exchange: $PS=0$")
    plt.legend()
    plt.show()


def Equations_2CXM():
    plt.figure()
    plt.figtext(0,.92,r"$v_p \frac{dC_p}{dt}=F_p c_a + PS\, c_e - (F_p +\, PS)\,c_p$", fontsize=20)
    plt.figtext(0,.82,r"$v_e \frac{dC_e}{dt}=PS\, c_p - PS\, c_e$", fontsize=20)
    plt.figtext(0,.67,r"$C(t)=I(t)\, \circ\ C_a(t)$", fontsize=20)
    plt.figtext(0,.57,r"$I(t)=F_+e^{-tK_+} + F_-e^{-tK_-}$", fontsize=20)
    
    plt.figure()
    plt.figtext(0,.92,r"$C(t)=I(t)\, \circ\ C_a(t)$   where   $I(t)=K^{trans}e^{-\frac{K^{trans}}{v_e}t}$", fontsize=20)

def Reproduce_OnTheScopeOfTofts():
    
    dt=1./60                     # sample at 1s time resolution
    t=numpy.arange(0,6, dt)      # time running from 0 to 6min
    ca = aif(t, aifchoice='Parker',zeropad=True)
     # intermediately vascularized
    vp = 10./100 # ml/ml
    ve = 20./100 # ml/ml
    PS = 25./100 # ml/min/ml
    Fpl = 50./100 # ml/min/ml
    ''' consider equation [32] and [33] to obtain a Ktrans,
    according to [32] it is Ktrans = 17  # ml/min/100ml '''
    Ktrans = PS*Fpl / (PS+Fpl) # ml/min/ml - eqn [32]
  

    c_2CXM = conc_2CXM( t, (ve, vp, PS, Fpl), ca, dt) 
    c_XTofts =  conc_XTofts(t,(vp, Ktrans, ve), ca, dt)
    c_Tofts =  conc_Tofts(t,(Ktrans, ve), ca, dt)
    plt.subplot(331)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,c_XTofts,'b:')
    plt.plot(t,c_Tofts,'g--')
    plt.axis([0,6,0,.4])

    '''now fit the data fit_generic(t, ca, data, p0, model):'''
    resTofts  = fit_generic(t, c_2CXM, conc_Tofts, (.17, 0.2), ca, dt)
    resXTofts = fit_generic(t, c_2CXM, conc_XTofts, (.1, .17, 0.2), ca, dt)
    res2CXM   = fit_generic(t, c_2CXM, conc_2CXM, (.2, .1, .25, .5), ca, dt)

    print resTofts[0]
    print resXTofts[0]
    print res2CXM[0]

    plt.subplot(332)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,resTofts[2],'g--')
    plt.plot(t,resXTofts[2],'b:')
    #plt.plot(t,res2CXM[2],'g-')


    # weakly vascularized
    vp = 1./100 # ml/ml
    ve = 29./100 # ml/ml
    PS = 25./100 # ml/min/ml
    Fpl = 50./100 # ml/min/ml
    Ktrans = PS*Fpl / (PS+Fpl) # ml/min/ml - eqn [32]
   
    c_2CXM = conc_2CXM(t,(ve, vp, PS, Fpl), ca, dt)
    c_XTofts =  conc_XTofts(t,(vp, Ktrans, ve), ca, dt)
    c_Tofts =  conc_Tofts(t,(Ktrans, ve), ca, dt)
    plt.subplot(334)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,c_XTofts,'b:')
    plt.plot(t,c_Tofts,'g--')
    plt.axis([0,6,0,.25])

    '''now fit the data fit_generic(t, ca, data, p0, model):'''
    resTofts  = fit_generic(t, c_2CXM, conc_Tofts, (.17, 0.2), ca, dt)
    resXTofts = fit_generic(t, c_2CXM, conc_XTofts, (.1, .17, 0.2), ca, dt)
    res2CXM   = fit_generic(t, c_2CXM, conc_2CXM, (.2, .1, .25, .5), ca, dt)

    print resTofts[0]
    print resXTofts[0]
    print res2CXM[0]
    
    plt.subplot(335)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,resTofts[2],'g--')
    plt.plot(t,resXTofts[2],'b:')

    # highly vascularized
    vp = 29./100 # ml/ml
    ve = 1./100 # ml/ml
    PS = 25./100 # ml/min/ml
    Fpl = 50./100 # ml/min/ml
    Ktrans = PS*Fpl / (PS+Fpl) # ml/min/ml - eqn [32]

    c_2CXM = conc_2CXM(t,(ve, vp, PS, Fpl), ca, dt)
    c_XTofts =  conc_XTofts(t,(vp, Ktrans, ve), ca, dt)
    c_Tofts =  conc_Tofts(t,(Ktrans, ve), ca, dt)
    plt.subplot(337)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,c_XTofts,'b:')
    plt.plot(t,c_Tofts,'g--')
    plt.axis([0,6,0,.4])

    '''now fit the data fit_generic(t, ca, data, p0, model):'''
    resTofts  = fit_generic(t, c_2CXM, conc_Tofts, (.17, 0.2), ca, dt)
    resXTofts = fit_generic(t, c_2CXM, conc_XTofts, (.1, .17, 0.2), ca, dt)
    res2CXM   = fit_generic(t, c_2CXM, conc_2CXM, (.2, .1, .25, .5), ca, dt)

    print resTofts[0]
    print resXTofts[0]
    print res2CXM[0]
    
    plt.subplot(338)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,resTofts[2],'g--')
    plt.plot(t,resXTofts[2],'b:')


    return (resTofts, resXTofts, res2CXM)
    
