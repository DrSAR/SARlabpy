# -*- coding: utf-8 -*-
"""
function library for fit_single_voxel

Created on Thu Sep 12 16:56:41 2013

@author: tammo
"""

from __future__ import print_function
from __future__ import division
import numpy as np
from scipy import optimize
import math
import pdb
import matplotlib
from matplotlib import rcParams
import matplotlib.pyplot as plt
import scipy.special as spsp



def create_model_AIF(params, time_axis, 
                     modify_rate = False,
                     shift_h = False, 
                     shift_v = False,
                     shift_arriv = False,
                     scale = False,
                     epsilon = None):
    ''' Create AIF for different rates (alphas), and fixed parameters
        from fit of data with AIF_model function. '''
                         
    def AIF_model(params, time):
        ''' Create AIF with linear upslope and single exponential decay.'''
        
        def heaviside(time, tau):
            """Heaviside function for use in AIF model."""
            out = np.zeros(len(time))
            i = 0
            for t in time:
                if t >= tau:
                    out[i] = 1
                i += 1
        
            return out
    #    tau = vol/rate
        alpha, beta, vol, c_f = params[0], params[1], params[2], params[3]
    
        tau = np.sqrt(2*vol/alpha)
        up = ( 1 - heaviside(time,tau) ) * alpha * time
    
        down = heaviside(time,tau) *\
            ( (alpha*(tau)-c_f) *\
            np.exp(-beta*(time-tau) ) + c_f  )
    
        return up+down     

        
    alpha, V1, V2, tf, cf = params[0], params[1], params[2], params[3], params[4]    
    
    tau = np.sqrt( (2*V1)/alpha )

        
    beta = (alpha*tau - cf) / (V2 - cf*(tf-tau))
    
    AIF = AIF_model([alpha, beta, V1, cf], time_axis)
    
    return AIF


def PKfit(aif, targetfunc, time, model='tofts', initial_params=[1,1],\
            method = 'leastsq', bounds = None, brute_ranges = None):
    """ Fits concentration curves derived by one of three different
        models (i.e. TM, ETM 2CXM) and an AIF to concentration data."""

    if method == 'leastsq':
        def errorfunc(params, time, aif):
            #print(params)
            if within_bounds(params, bounds) or bounds == None:
#                print(targetfunc - PKmodel(params, time, aif, model))
                return targetfunc - PKmodel(params, time, aif, model)
            else:
                return 1e9
                
        fit_params, cov, infodict, mesg, sucess = optimize.leastsq(errorfunc, \
                                                initial_params, \
                                                args=(time,aif),
                                                full_output=True,
                                                maxfev=200000)
                                                
#        print('\n\n\t',sucess)
#        print('\n\n\t',mesg)
#        print('\n\n\t',infodict['nfev'])
#        average_residual = sum((targetfunc - PKmodel(fit_params, time, aif, model))**2)
#        average_residual /= len(time)
                                                
    if method == 'anneal':
        def errorfunc(params, time, aif):
            #print(params)
            if within_bounds(params, bounds):
#                print(targetfunc)
#                print(sum(abs(targetfunc - PKmodel(params, time, aif, model))**2))
#                print('in')
                return sum(abs(targetfunc - PKmodel(params, time, aif, model))**2)
            else:
#                print('out')
                return 1e9
                
            #fast or cauchy are available
        fit_params, Jmin, T, feval, iters, accept, retval = optimize.anneal(\
                                            func = errorfunc,\
                                            x0 = initial_params, \
                                            args =(time,aif),\
                                            schedule='fast',\
                                            full_output=True,\
                                            T0=None,\
                                            Tf=1e-10,\
                                            maxeval=None,\
                                            maxaccept=None,\
                                            maxiter=500,\
                                            boltzmann=.5,\
                                            learn_rate=0.1,\
                                            feps=1e-6,\
                                            quench=1.0,\
                                            m=1.0, n=1.0,\
                                            lower=[i[0] for i in bounds],\
                                            upper=[i[1] for i in bounds],\
                                            dwell=1000)
#        print(feval, iters, retval, T)
        
        
    elif method == 'brute':
        
        def errorfunc(params, time, aif):
            out = sum(abs(targetfunc - PKmodel(params, time, aif, model))**2)
            return out            
            
        fit_params = optimize.brute(errorfunc,\
                            ranges = brute_ranges,\
                            args=(time,aif))
    
    elif method == 'fmin_tnc':
        if bounds == None:
            raise TypeError('Method fmin_tnc needs parameter bounds.')

        errorfunc = lambda params, time, aif:\
        sum(abs(targetfunc - PKmodel(params, time, aif, model)))
        
        fit_params, nfeval, infodict = optimize.fmin_tnc(func = errorfunc, \
                                                x0 = initial_params, \
                                                fprime = None,\
                                                args=(time,aif),\
                                                approx_grad=True,\
                                                bounds = bounds,\
                                                epsilon = 1e-10,\
                                                accuracy = 1e-12,\
                                                maxfun = 20000)
                                                
        average_residual = sum((targetfunc - PKmodel(fit_params, time, aif, model))**2)
        average_residual /= len(time)

        
    best_fit = PKmodel(fit_params, time, aif, model)        
    
    # Calculate Coefficient of Determination
    ss_err=((best_fit - targetfunc)**2).sum()
    ss_tot=((targetfunc - np.array(targetfunc).mean())**2).sum()
    rsquared= 1 - (ss_err/ss_tot)
    
    N = len(best_fit)
    K = len(fit_params)
    AIC = N*np.log(ss_err/N) + 2 * K + ( (2*K*(K+1)) / (N-K-1) )
    BIC = N*np.log(ss_err/N) + K * np.log(N)
    
#    fit_params = initial_params
        
    return fit_params, rsquared, AIC, BIC, best_fit, ss_err
    

def within_bounds(params, bounds):
    if params is None:
        return False
    if bounds is None:
        return True
    else:
        for i in xrange(len(params)):
            if params[i] < bounds[i][0] or params[i] > bounds[i][1]:
                #print('PARAMETER OUT OF BOUDNS')
                return False
        return True
        


def PKmodel(params, time, aif, model='tofts', vasc=2):    
    out = 1
    if model == 'tofts':
        """ Calculate irf and convolve with aif according to TM. """
        
        # convert K_trans to 1/s
        K_trans, v_e = params[0]/60, params[1]
            
        irf = K_trans * np.exp( ( - time * K_trans) / v_e )
        out = convolve_aif ( aif, irf, time )

    elif model == 'etofts':
        K_trans, v_e, v_p = params[0]/60, params[1], params[2]
        irf = K_trans * np.exp( ( - time * K_trans) / v_e )    
        out = ((v_p * aif) + convolve_aif (aif, irf, time))
        
    elif model == 'aath':
#        E, F_p, v_e, t_c = params[0], params[1]/60, params[2], params[3]
        # parameterize with v_p for easier comparsion with 2cxm
        E, F_p, v_e, v_p = params[0], params[1]/60, params[2], params[3]
        
        t_c = v_p/F_p        
        
        irf = np.zeros(len(time))
        
        for i in xrange(len(time)):
            if time[i] > t_c:
                irf[i] = E*F_p * np.exp( ( - (time[i] - t_c) * E*F_p) / v_e )
            elif time[i] <= t_c and time[i] > 0:
                irf[i] = F_p
            elif time[i] < 0:
                irf[i] = 0
            i += 1
        
        out = convolve_aif ( aif, irf, time )


    elif model == '2cxm':
        ''' Make sure E is not to close to one. Will cause Nans! '''
        
#        raise TypeError('2cxm not working yet. FIX LATER!')
        #===== Unpack Paramtes ====
        E, F_P, V_E, V_P = params[0], params[1]/60, params[2], params[3]
        #print(params)
        #==== Define variables ====
        e = V_E / ( V_P + V_E )
        
        tau_factor = ( E + e - ( E * e ) ) / ( 2 * E )
        
        tau_frac_numerator   = E * e * ( 1 - E ) * ( 1 - e )
        tau_frac_denominator = np.power( E + e - ( E * e ), 2)
        tau_frac             = tau_frac_numerator / tau_frac_denominator
        
        tau_plus  = tau_factor * ( 1 + np.sqrt ( 1 - 4 * tau_frac ) )
        tau_minus = tau_factor * ( 1 - np.sqrt ( 1 - 4 * tau_frac ) )
        
        v_tot   = V_P + V_E
        K_plus  = F_P / ( v_tot * tau_minus )
        K_minus = F_P / ( v_tot * tau_plus  )
        
        tau_diff = tau_plus - tau_minus
        
        F_Plus   =   ( F_P * ( tau_plus  - 1 ) )  / ( tau_diff )
        F_minus  = - ( F_P * ( tau_minus - 1 ) )  / ( tau_diff )
        
        #==== Compute irf ====
        if vasc == 2:
            irf = F_Plus  * np.exp ( - time * K_plus  ) + \
                  F_minus * np.exp ( - time * K_minus )
        
        elif vasc == 3:
            irf = F_P * np.exp ( -time * ( F_P / V_P ) )
        
        elif vasc == 1:
            PS          = ( F_P * E ) / ( 1 - E )
            weak_factor = ( F_P * PS ) / ( F_P + PS )
            irf         = weak_factor * np.exp ( ( - time * weak_factor ) / V_E )
        
        else:
            print ("No valid VASCULARIZATION type given.\n \
                    Try VASCULARIZATION = 'intermediate', 'high' or 'weak'.\n \
                    Intermediate was used here to avoid runtime error.")
            irf = F_Plus  * np.exp ( - time * K_plus  ) + \
                  F_minus * np.exp ( - time * K_minus )
        
        #high_vasc = F_P * np.exp ( - time * ( F_P / V_P ) )
        #weak_vasc = ( ( F_P * PS ) / ( F_P + PS ) ) *\
        #np.exp( - time * ( ( F_P * PS ) / ( (F_P + PS) * V_E ) ) )
        out = convolve_aif ( aif, irf, time )
        
    elif model == '2cum':
        
        E, F_P, V_P = params[0], params[1]/60, params[2]
        
        K = E * F_P        
        
        PS = ( E * F_P ) / (1 - E)
        
        Tp = (V_P * E) / PS        
        
        irf = F_P * np.exp(-time/Tp) + K * ( 1 - np.exp(-time/Tp) )
        
        out = convolve_aif ( aif, irf, time )
        
    elif model == 'pat':
        K, v_p = params[0]/60, params[1]
        
        irf = [K]*len(time)
        
        out = ((v_p * aif) + convolve_aif (aif, irf, time))
        
    elif model == 'distr':
        
        E, F_p, v_e, v_p = params[0], params[1]/60, params[2], params[3]
        
        t_c = v_p/F_p
        
        PS = - F_p * np.log(1 - E)
#        PS = F_p*E
        
        irf = np.zeros(len(time))
        
                
        for i in xrange(len(time)):
            
            if time[i] > t_c:
                
                
                delta_time = time[1] - time[0]

                integral = 0
                
                kappa = PS**2 / (v_e*F_p)
                
                for tau in time[1:len(time)]:
                    
                    if tau > time[i]-t_c:
                        break;
                                        
                    integral += delta_time*\
                        (np.exp(-(PS*tau)/v_e)*np.sqrt( kappa / tau ) *\
                         spsp.i1(2*np.sqrt(kappa*tau)))
                
                
                irf[i] = F_p*(1-(np.exp(-PS/F_p)*(1+integral)))
#                irf[i] = np.exp(-E) ( 1 + Integral )
#                    
                    
            elif time[i] <= t_c and time[i] >= 0:
                irf[i] = F_p
                
            elif time[i] < 0:
                irf[i] = 0
                            
            i += 1

        out = convolve_aif ( aif, irf, time )

    elif model == 'PFM':
        V, F_p = params[0], params[1]/60, 
        
        irf = np.zeros(len(time))
        
        for i in xrange(len(time)):
          
          if V/F_p > time[i]:
                irf[i] = F_p
          else:
                irf[i] = 0
                
        out = convolve_aif ( aif, irf, time )
    
    return out
    
    
def SPGR_sgnl_from_ctr(sgnl, exp_params = [4.99e-3, 100, 2000, 25*2*np.pi/360]):
    
    ''' Convert CA concetrnation to SPGR/FLASH signal (alpha in rad).'''
    r = exp_params[0]
    TR = exp_params[1]
    T10 = exp_params[2]
    alpha = exp_params[3]
    
    E = ( np.sin(alpha) - sgnl ) / ( np.sin(alpha) - sgnl*np.cos(alpha) )
    T1 = - TR / np.log(E)
    ctr = (1/r) * ( (1/T1) - (1/T10) )
    
    return ctr



def convolve_aif ( aif, irf, time ):
    """ Calculate a discrete convolution
        (adds zeros to arrays to user np.convolve. """

    T_STEP = ( max(time) - min(time) ) / ( len(time) - 1 )

    aif = np.insert(aif,np.zeros(len(irf)-1),0)
# Stefans solution
#    AIF = numpy.hstack((numpy.zeros(time_axis.shape[0]+t0-1), 
#                       (AIFvfunc(time_axis))[0:-t0]))    
    return np.convolve(irf,aif,mode='valid')*T_STEP
