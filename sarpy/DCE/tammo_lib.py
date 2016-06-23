# -*- coding: utf-8 -*-
"""
This library provides functions to analyze MRI data.
Created on Wed Sep 19 19:19:42 2012 @author: tammo
"""
from __future__ import print_function
from __future__ import division
import numpy as np
import sys
import csv
from scipy import optimize
import matplotlib.pyplot as plt
#import time
import math

def movingaverage(x, N):
    '''
    this replaces movingaverage from movingaverage to cut down on external
    dependencies
    '''
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / N

def AIF_from_TDM( injection, time, parms ):
    '''Calculates AIF from a injection curve and a Tracer distribution
        model H_TDM(parms) by convolution.'''

    ### fraction of max time added (and later cutoff from AIF)
    cutoff_frac = 2

    ### Calculate time step
    d = time[1]-time[0]

    ### double time array, cut off later
    time_long = np.arange(time[0], cutoff_frac*max(time), d)

    freq = np.fft.fftfreq( len(time_long), d)
	#AIF = convolve_aif( abs(np.fft.ifft(H_TDM(freq, parms))), injection ,time)
    AIF = convolve_aif( abs(np.fft.ifft(H_TDM(freq, parms))), injection, time_long )
    AIF = AIF[:len(time)]
    return AIF #, freq


def convolve_aif ( aif, irf, time ):
    """ Calculate a discrete convolution
        (adds zeros to arrays to user np.convolve. """

    T_STEP = ( max(time) - min(time) ) / ( len(time) - 1 )

    aif = np.insert(aif,np.zeros(len(irf)-1),0)
    return np.convolve(irf,aif,mode='valid')*T_STEP


def H_TDM (freq, parms):
    ''' Calculate the TDM (Tracer distr. model) transfer function. '''
    omega = 2*np.pi*freq
    A0, A1, A2, B0, B1, B2, B3 = parms[:]
    num = B3*omega**3 + B2*omega**2 + B1*omega + B0
    den = omega**3 + A2*omega**2 + A1*omega + A0
    return num/den


def Trapezoid (xdata, parms):
    ''' Trapezoid shaped function to simluate CA injection profiles.
    parms = volume, xmax, slope - of the injection.
    Interpolates values for timespacing smaller than 1. This
    might be fixed later.'''

    delta_x_final = xdata[1]-xdata[0]
    xdata_final = xdata
    xdata = range(int(math.ceil(np.max(xdata))))
    v, x, s = parms[:]
    s /= delta_x_final

    # Determine maximum injection rate
    y = ( ( x*s ) / 2 ) - np.sqrt( ( ( x*s ) / 2 )**2 - ( s*v ) )
    out = np.zeros(len(xdata))

    # Determine beginning and end of plateau
    x_l = y/s
    x_r = x - x_l

    i = 0
    while xdata[i] <= x: # Iterate all data points during injection
        if xdata[i] < x_l:
            out[i] = s*xdata[i]
        elif xdata[i] > x_r:
            out[i] = (x-xdata[i])*s
        elif xdata[i] >= x_l and xdata[i] <= x_r:
            out[i] = y
        i += 1

    out = np.interp(xdata_final, xdata, out)

    return out


def scale_aif(aif_factor, p_aif, t_aif):
    """ Scales AIF according to aif_fit_function, for given parameters
        and scaling parameter. """

    a_1 =       p_aif[0]
    sigma_1 =   p_aif[1]#*aif_factor
    mu_1 =      p_aif[2]
    a_2 =       p_aif[3]
    sigma_2 =   p_aif[4]
    mu_2 =      p_aif[5]
    alpha =     p_aif[6]
    beta =      p_aif[7]*aif_factor # adjust final concentration
    s =         p_aif[8]
    tau =       p_aif[9]

    # scale aif width by scaling factor
    t_aif =     t_aif*(1/(aif_factor*5))

    p_aif = [ a_1, sigma_1, mu_1, a_2, sigma_2, mu_2, alpha, beta, s, tau ]
    aif = aif_fit_function(t_aif, p_aif)

    aif_shift = aif[len(aif)-1]

    ## shift aif to have the final concentration equal zero
    aif -= aif_shift

    ## sacle aif height by scaling factor
    aif /= aif_factor

    ## roll aif values smaller than zero to the end
    i = 0
    while aif[i] < 0:
        i += 1
    aif = np.roll(aif, -i)

    ## set negative aif values at the end to final concentration
    i = len(aif) - 1
    while aif[i] < 0:
        i -= 1
    aif[i:] = aif[i-1]

    ## shift aif back up
    aif += aif_shift

    return aif


def import_AIF(path, aif_time_min = 0, aif_time_spacing = 1, aif_dump = 0):
    ''' Load AIF data an create time array.
        Delete data points in the beginning. '''
    aif = import_csv(path)

    aif_time_max = len(aif)*aif_time_spacing

    aif_time = np.arange( aif_time_min, aif_time_max, aif_time_spacing )

    aif_time, aif = aif_time[:-aif_dump], aif[aif_dump:]

    return aif, aif_time



def T1_from_FLASH (signal, T1_0, TR, alpha, r = 4.3e-3):
    ''' Derive (DELTA) T1 values from FLASH signal.
        Input:  
            signal - array with (relative) signal values
            T1_0 - baseline T1 values, same dimension as signal
            TR - Repetition time
            alpha - Flip angle in grade
            r - Factor (here for Gd)
        Output: 
            T1 array, with same dimensions as signal array.
    '''

    deg_to_rad = 2*np.pi/360
    alpha *= deg_to_rad
    angle_factor = np.cos(alpha)
    E_0 = np.exp(-TR/T1_0)

    chi = signal*(E_0 - 1)
    xi = (E_0*angle_factor) - 1

    E1 = ( xi - chi ) / ( xi - chi*angle_factor )
    T1 = -TR / np.log(E1)

    c_t = (1/r) * ((1/T1) - (1/T1_0))

    return T1, c_t


def plot_image(data, time_points):
    """ Plots images for one ore several time points. """

    images = len(time_points)

    if images == 1:
        plt.imshow(data[:, :, time_points[0] ], interpolation='none')
        plt.colorbar()
        plt.savefig(get_current_date() + '_image' )
        plt.show()

    else:

        data = data[:, :, time_points]

        image_length = np.sqrt(images)
        fig = plt.figure()
        for i in np.arange(data.shape[2]):
            fig.add_subplot(image_length, image_length, i).imshow(\
                data[:,:,i], interpolation='none',\
                vmin = np.nanmin(data[:, :, :]),\
                vmax = np.nanmax(data[:, :, :]) )

        plt.savefig(get_current_date() + '_image' )
        plt.show()



def plot_signal_time_map(signal, slice, marker='-'):
    """Plots a pseudo-3D pixel-wise time-evolution map
        from a 4D array ( [y, x, slice/angle, time] ) for a chosen slice """



    fig = plt.figure(facecolor='white')
    fig.canvas.set_window_title('T1(t)-map')
    ax_main = plt.axes(frameon=True, title ='signal(t)-map')
    ax_main.get_xaxis().set_visible(False)      # Hidex x-axis
    ax_main.get_yaxis().set_visible(False)      # Hide y-axis

    # Iterate of pixels (subplots)
    n=0
    for j in range(signal.shape[0]):
        for i in range(signal.shape[1]):
            ax1 = fig.add_subplot(signal.shape[0], signal.shape[1], n)
            ax1.plot(signal[j, i, slice, :], marker)
            ax1.set_frame_on(False)                     # Hide frame
            ax1.axes.get_xaxis().set_visible(False)     # Hide y-axis
            ax1.axes.get_yaxis().set_visible(False)     # Hide x-axis
            ax1.axes.get_xaxis().tick_bottom()          # Turn of ticks at top
            n += 1

    plt.savefig( get_current_date() + '_signal_time_map' )
    plt.show()



def divide_by_baseline(data, baseline_range = 10):
    baseline_range = range(baseline_range) # number of values considered as baseline (preinjection)
    for i in range(data.shape[0]):    # Iterate over pixels and slices
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                # Divide signal by baseline
                data[i,j,k,:] = data[i,j,k,:]/\
                                  np.mean( data[i, j, k, baseline_range] )

    return data



def timeline_moving_average(data, ave_window = 10):

    ave_signal = np.zeros([data.shape[0], data.shape[1],data.shape[2],\

                            data.shape[3]-ave_window+1])
    for i in range(data.shape[0]):    # Iterate over pixels and slices
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                            # calculate moving average
                ave_signal[i,j,k,:] = np.array(list(\
                                movingaverage(
                                data[i,j,k,:], ave_window)))
    return ave_signal






def filter_noise( data, bins=15, plots = 'Noisy', threshold = 1.10,\
                 figsize=(22,8), data_start = 0, data_end = 'full' ):
    """Sophisticated filter to get rid of noisy pixels.
        Histograms from data[i,j,k,:] (concentration time curves) with
        bins are created. Gaussians are fitted to the histogram.
        The gaussian mean is divided by the histograms median bar.
        If the ratio is smaller than 1.10 (less increase) the pixel
        is dismissed.
        Maybe it is interesting to work with the histograms to examine
        the type of noise and to take the chi2 of the gaussian fit as
        filter criterion.
        data_start, data_end allow you to take only parts of the C(t)
        curve into account. """

    def gaussian(parms, x):
        """"Gaussian function for fitting """
        mu = parms[0]
        sigma = parms[1]
        height = parms[2]
        return height * np.exp(-( (x-mu)/sigma )**2 )

    def err(parms, x, hist_data):
        """ Cost function for fitting """
        return(gaussian(parms, x) - hist_data)


    def plot_c_t_curve_and_hist():
        """ Plots concentration-time-curve and histogram with fit. """
        fig = plt.figure(1, figsize=figsize)
        ax1 = fig.add_subplot(121)
        ax1.plot(data_temp)

        ax2 = fig.add_subplot(122)
        ax2.hist(data_temp, len(hist_data))
        ax2.plot(x_hist, gaussian(p_fit, x_hist),'r-o')

        ax2.text(x_hist.min(), (9/10)*hist_data.max(), \
                'chi2 per bar = '+str(chi2/len(hist_data))+\
                '\n gauss_mean / hist_median = '+str(mean_med_ratio)+\
                '\n gauss_pars = '+str(p_fit))
        plt.show()


    if data_end == 'full':
        """ Set data end to maximal value. """
        data_end = data.shape[3]


    out = np.zeros(data.shape)  # Create output array (not necessary)

    # Iterate over all pixels and slices
    for i in np.arange(data.shape[0]):
            for j in np.arange(data.shape[1]):
                    for k in np.arange(data.shape[2]):

                        # Cut out time-interval of interest
                        data_temp = data[i, j, k, data_start:data_end]

                        # Estimate Gaussian fit paramters
                        # (Order of magn. guesses)
                        mu = np.median(data_temp)
                        sigma = np.median(data_temp)
                        height = np.max(data_temp)
                        p0 = [mu, sigma, height]

                        # Create histgramm of C-time-curve
                        hist_data, hist_edges = np.histogram(data_temp, bins=bins)

                        # x-values for histogram (only for plotting)
                        x_hist = np.array(list(\
                                    movingaverage(hist_edges, 2)))

                        # Fit gaussian to histogram
                        p_fit, sucess = optimize.leastsq(\
                            err, p0, args=(x_hist, hist_data), maxfev=800)

                        # Fitted gaussian
                        gaussian_result = gaussian(p_fit, x_hist)
                        chi2 = sum( (gaussian_result - hist_data)**2 )

                        # Calculate ratio mean over x-median
                        mean_med_ratio = p_fit[0]/np.median(x_hist)

                        # Plot commands for noisy c-t-curves
                        if mean_med_ratio <= threshold:

                            out[i,j,k,:] = float('NaN')

                            if plots == 'Noisy' or plots == 'All':
                                plot_c_t_curve_and_hist()


                        # Plot commands for non-noisy c-t-cuves
                        else:

                            out[i,j,k,:] = data[i,j,k,:]

                            if plots == 'All':
                                plot_c_t_curve_and_hist()

    return out



def fit(aif, targetfunc, time, model, initial_parms):
    """ Fits concentration curves derived by one of three different
        models (i.e. TM, ETM 2CXM) and an AIF to concentration data."""



    if model == 'two_compartment_xm':

        errorfunc = lambda parms, time, aif: \
            targetfunc - two_compartment_xm(parms, time, aif, \
                                            VASCULARIZATION = 'intermediate')
        p_fit, sucess = optimize.leastsq(errorfunc, \
                                            initial_parms, \
                                            args=(time,aif))

        chi2 = sum((targetfunc - two_compartment_xm(p_fit, time, aif))**2)
        chi2 = chi2/len(time)

    if model == 'tofts_model':

        errorfunc = lambda parms, time, aif: \
            targetfunc - tofts_model(parms, time, aif)

        p_fit, sucess = optimize.leastsq(errorfunc,\
                initial_parms, args=(time,aif))

#        p_fit = optimize.fmin_slsqp(errorfunc,\
#                initial_parms, args=(time,aif), iter = 10000)

#        sucess = 1

        chi2 = sum((targetfunc - tofts_model(p_fit, time, aif))**2)
        chi2 = chi2/len(time)

    elif model == 'extendet_tofts_model':

        errorfunc = lambda parms, time, aif: \
            targetfunc - extendet_tofts_model(parms, time, aif)

        p_fit, sucess = optimize.leastsq(errorfunc,\
                initial_parms, args=(time,aif))

        chi2 = sum((targetfunc - extendet_tofts_model(p_fit, time, aif))**2)
        chi2 = chi2/len(time)

    return p_fit, sucess, chi2



def aif_fit_function(time, parms = [.809, .0563, .17046, .330, .132, .365, 1.050, .1685, 38.078, .483] ):
    """ Calculates aif according to the fitfunction by Parker et. al. """
    a_1 = parms[0]
    sigma_1 = parms[1]
    mu_1 = parms[2]
    a_2 = parms[3]
    sigma_2 = parms[4]
    mu_2 = parms[5]
    alpha = parms[6]
    beta = parms[7]
    s = parms[8]
    tau = parms[9]

    out = a_1 * gaussian( time, sigma = sigma_1, mu = mu_1 ) +\
          a_2 * gaussian( time, sigma = sigma_2, mu = mu_2 ) +\
          alpha * np.exp(- beta * time ) / (1 + np.exp(-s * ( time - tau )))

    return out



def aif_fit_tammo( time, parms = [1, 1, 1, 1] ):
    """ Try own AIF model. """
    out = parms[0]*time*heaviside(-time, -parms[1]) +\
        parms[0]*parms[1]*np.exp(-(time-parms[1])*parms[3]) * heaviside(time, parms[1])
    return out


def heaviside(time, tau):
    """Heaviside function for use in AIF model."""
    out = np.zeros(len(time))
    i = 0
    for t in time:
        if t >= tau:
            out[i] = 1
        i += 1

    return out





def tofts_model ( parms, time, aif ):
    """ Calculate irf and convolve with aif according to TM. """

    # Unpack Paramtes ########### old ###############
    # E, F_P, V_E, V_P = parms[0], parms[1], parms[2], parms[3]
    # irf = E * F_P * np.exp( ( - time * E * F_P ) / V_E )
    #################################################

    K_trans, V_E = parms[0], parms[1]

    irf = K_trans * np.exp( ( - time * K_trans) / V_E )

    return convolve_aif ( aif, irf, time )


def extendet_tofts_model ( parms, time, aif ):
    """ Calculate irf and convolve with aif according to ETM. """

    #===== Unpack Paramtes ====
    E, F_P, V_E, V_P = parms[0], parms[1], parms[2], parms[3]


    irf = E * F_P * np.exp( ( - time * E * F_P ) / V_E )

    return ((V_P * aif) + convolve_aif (aif, irf, time))



def two_compartment_xm ( parms, time, aif, \
                         VASCULARIZATION = None ):

    """ Returns the irf for the 2CXM. """

    #===== Unpack Paramtes ====
    E, F_P, V_E, V_P = parms[0], parms[1], parms[2], parms[3]

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
    if VASCULARIZATION == 'intermediate':
        irf = F_Plus  * np.exp ( - time * K_plus  ) + \
              F_minus * np.exp ( - time * K_minus )

    elif VASCULARIZATION == 'high':
        irf = F_P * np.exp ( -time * ( F_P / V_P ) )

    elif VASCULARIZATION == 'weak':
        PS          = ( F_P * E ) / ( 1 - E )
        weak_factor = ( F_P * PS ) / ( F_P + PS )
        irf         = weak_factor * np.exp ( ( - time * weak_factor ) / V_E )

    else:
        print ("No valid VASCULARIZATION type given.\n \
                Try VASCULARIZATION = 'intermediate', 'high' or 'weak'.\n \
                Intermediate was used here to avoid runtime error.",\
                file=sys.stderr)
        irf = F_Plus  * np.exp ( - time * K_plus  ) + \
              F_minus * np.exp ( - time * K_minus )

    #high_vasc = F_P * np.exp ( - time * ( F_P / V_P ) )
    #weak_vasc = ( ( F_P * PS ) / ( F_P + PS ) ) *\
    #np.exp( - time * ( ( F_P * PS ) / ( (F_P + PS) * V_E ) ) )

    return convolve_aif ( aif, irf, time )


def import_csv (csv_file) :
    """ Imports single row csv to list. """

    reader = csv.reader( open( csv_file, 'rU' ) )
    out_y = []

    for row in reader:
        out_y.append(float(*row))

    return np.array(out_y)


def T1_from_lin_fit(p1, TR):
    """ Returns TR and M0 from a linear fit a pre inj data at different angles."""
    T1 = - TR / np.log(p1[0])
    M0 = p1[1] / (1 - p1[0])
    return (T1, M0)



def average_mri_4d(orig_data, ave_len = 4):
    """ Calculate array of averaged mri image with dimensions [x,y,slice,angle]."""

    #   ave_len: side length of area to average

    #   prepare new array
    ave_data = np.zeros([ orig_data.shape[0]/ave_len,\
                            orig_data.shape[1]/ave_len,\
                            orig_data.shape[2],\
                            orig_data.shape[3] ])

    #   Iterate over angles and slices
    for i in range(orig_data.shape[2]):
        for j in range(orig_data.shape[3]):
            #   Iterate over areas to average
            for x_area in range(orig_data.shape[0]//ave_len):
                for y_area in range(orig_data.shape[1]//ave_len):

                    # Creat indices
                    x_start, x_stop = x_area*ave_len, (x_area+1)*ave_len
                    y_start, y_stop = y_area*ave_len, (y_area+1)*ave_len

                    # Select data and write average to new array
                    ave = orig_data[ x_start:x_stop, y_start:y_stop, i, j]
                    ave_data[x_area, y_area, i, j] = np.average(ave)
    return ave_data


def get_current_date():
    """ Determine date-string for file-name. """

    date = time.localtime()

    current_date = str(date[0])+'_'+str(date[1])+'_'+str(date[2])+\
        '_'+str(date[3])+'_'+str(date[4])+'_'+str(date[5])

    return current_date



def gaussian(x, sigma=1, mu=0):
    """ Returns a Gaussian. """

    return (sigma * np.sqrt( 2 * np.pi ) )**(-1) * \
         np.exp( - ( x - mu )**2 /  ( 2* sigma**2 ) )



#    if x_max != len(irf):
#         print( "Arrays have different length", file=sys.stderr )



#==============================================================================
# def convolve_aif_old ( aif, irf, time ):
#     """ Calculate a discrete convolution:
#         C(i) = \sum\limits_{j=0}^i aif(j) irf(i-j) """
#     i = 0
#     C = np.zeros( len ( aif ) )
#     T_STEP = ( max(time) - min(time) ) / ( len(time) - 1 )
# #    print(T_STEP)
#
#     while i < len( aif ):
#
#         C_temp, j = 0, 0
#
#         while j <= i:
#
#             C_temp   +=  aif[ j ] * irf [ i - j ]
#             j        += 1
#
#         C [ i ] = C_temp * T_STEP
#        i = i + 1
#    return C
#==============================================================================



#==============================================================================
# ##class par:
# ##    def __init__(self, value):
# ##            self.value = value
# ##
# ##    def set(self, value):
# ##            self.value = value
# ##
# ##    def __call__(self):
# ##            return self.value
# ##
# ##def fit(function, parameters, y, x = None):
# ##
# ##    def errorfunc(params):
# ##
# ##        i = 0
# ##        for p in parameters:
# ##            p.set(params[i])
# ##            i += 1
# ##        return y - function(x)
# ##
# ##    if x is None: x = arange(y.shape[0])
# ##
# ##    p = [param() for param in parameters]
# ##    print(p)
# ##    return optimize.leastsq(errorfunc, p)
#==============================================================================
#
