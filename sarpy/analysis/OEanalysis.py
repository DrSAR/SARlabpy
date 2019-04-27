import numpy
import sarpy
import scipy.integrate
import scipy.optimize
import scipy.fftpack
import scipy.stats
import copy
import os
import json
import datetime
import collections
import sarpy.analysis.getters as getters
import imp
import pylab
import lmfit

def definite(x,A,mean,sigma):
    # This is the lognormal definite and indefinite integral
    # https://www.wolframalpha.com/input/?i=integrate+(A_0%2F2)*erfc(+(log(z%2Fx)+-+u)%2F(s_0*sqrt(2)))+dz
    # this is the result of eq'n 8 integrated dz. This is from the john hudson paper : 
    # QUANTIFICATION OF FLOW USING ULTRASOUND AND MICROBUBBLES: A DISRUPTION REPLENISHMENT MODEL BASED ON PHYSICAL PRINCIPLES
    # bear://x-callback-url/open-note?id=F0BE0087-CB51-4236-BCB0-F8DC97A53F06-331-000013EDCE1953C1

    def indefiniteIntegral(x,A,mean,sigma,z):
        
        part1 = z*scipy.special.erfc((numpy.log(z/x)-mean)/(numpy.sqrt(2)*sigma))
        part2 = x*numpy.exp(mean+sigma**2/2)*scipy.special.erf((numpy.log(z/x)-mean-sigma**2)/(numpy.sqrt(2)*sigma))
        # P.S. according to wolfram alpha, the negative sign on part 2 is in a different place inside the erf. but erf is an odd function so erf(-x) = -erf(x)
        # https://en.wikipedia.org/wiki/Error_function
    
        return A/2*(part1+part2)
    
    res = indefiniteIntegral(x,A,mean,sigma,1) - indefiniteIntegral(x,A,mean,sigma,0)
    
    # to try and avoid fit failures
    return numpy.nan_to_num(res)

def fitParams_lognormal(parameterObject,timeVector):

    from lmfit.models import ExponentialModel,ConstantModel
    from lmfit import Model

    fit_model = Model(definite)+ConstantModel()

    return fit_model.eval(parameterObject,x=timeVector)


def fit_dOEMRI_EC(scn_to_analyze,
                  adata_label,
                  secondsToAnalyse = 90,
                  viz=True,
                  fitStart=None,
                  printFitResult=True,
                  weights=None,
                  startingParamDict = None):

    
    from lmfit.models import ExponentialModel,ConstantModel
    from lmfit import Model
    
    scn = sarpy.Scan(scn_to_analyze)

    # Specify the data being used
    time_per_rep = sarpy.analysis.getters.get_time_per_rep(scn_to_analyze)
    numpoints = int(secondsToAnalyse/time_per_rep)
    
    data = scn.adata[adata_label].data
    if fitStart is None:
        fitStart = sarpy.helpers.find_o2_inj_point(data)
    yd = data[fitStart:fitStart+numpoints]
    td = numpy.arange(0,numpoints)*time_per_rep/60+1E-4 # 1E-4 is added so it's not evaluated at 0

    ################## Fitting
    
    modelt = numpy.linspace(0,numpoints*time_per_rep,100)/60

    # Set up the model first
    fit_model = Model(definite)+ConstantModel()
    pars = fit_model.make_params()
    
    if startingParamDict is None:
        startingParamDict={'A':0.4,
                           'mean':0.2,
                           'sigma':0.3,
                           'c':0.2}
    
    pars['A'].set(startingParamDict['A'], min=-1E5, max=1E5)
    pars['mean'].set(startingParamDict['mean'], min=0,max=50)
    pars['sigma'].set(startingParamDict['sigma'], min=0,max=50)
    pars['c'].set(yd[0], min=yd[0]*0.9,max=yd[0]*1.1,vary=False)
    
    # Weight the first 5 points 5x more
    if weights is not None:
        raise NotImplementedError('this feature has not yet been implemented')
    else:    
        weights = numpy.ones(numpoints)
        #weights[0:5] = 2*weights[0:5]

    #Grab the initial fit
    init = fit_model.eval(pars,x=modelt)

    # Do the fit
    result = fit_model.fit(yd,pars,x=td,weights=weights)
    
    if viz:
        display_dOEMRI_fit(result,time_per_rep)
      
    if printFitResult:
        print(result.fit_report())
    
    return result

def display_dOEMRI_fit(fitresult,time_per_rep):
    
    modelt = numpy.linspace(0,len(fitresult.data),100)*time_per_rep/60+1E-4 # 1E-4 is added so it's not evaluated at 0
    time_data = numpy.arange(0,len(fitresult.data))*time_per_rep/60+1E-4 # 1E-4 is added so it's not evaluated at 0
    #ax = pylab.subplots(111)
    
    pylab.figure(figsize=(8,6))
    # Plot the results
    pylab.plot(time_data,fitresult.data,marker='x',linewidth=0,label='raw') #raw data
    pylab.plot(modelt,fitresult.model.eval(fitresult.params,x=modelt), 'r-',label='fit')    
    pylab.plot(modelt,fitresult.model.eval(fitresult.init_params,x=modelt), 'g--',label='init')    

    pylab.xlabel('Time (min)',fontsize=18)
    pylab.ylabel('Normalized Component \n Strength',fontsize=18)
    pylab.title('Fitting Oxygen replenishment',fontsize=24)
    
    try:
        pylab.text(0.5,-0.05,'$\mu$ = {0:.2f}$\pm${1:.2f} mm/min'.format(fitresult.params['mean'].value,fitresult.params['mean'].stderr),fontsize=14)
        pylab.text(0.5,-0.07,'$\sigma_f$ = {0:.2f}$\pm${1:.2f} mm/min'.format(fitresult.params['sigma'].value,fitresult.params['sigma'].stderr),fontsize=14)
        pylab.text(0.5,-0.09,'A = {0:.2f}$\pm${1:.2f}'.format(fitresult.params['A'].value,fitresult.params['A'].stderr),fontsize=14)

    except TypeError:
        print('Uncertainty could not be calculated for one or more fit params')
    pylab.legend()

    return None