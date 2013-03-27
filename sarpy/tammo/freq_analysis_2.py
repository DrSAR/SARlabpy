# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 16:31:55 2012

TODO
- AIFs are constant in the end, should slowly decrease!
- Reproduce Aerts et al.
- FFTSHIFT ?!? (like in Matlab)

@author: tammo
"""

from __future__ import print_function, division
import tammo_lib as tl
from scipy import optimize
import scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import norm
import matplotlib.mlab as mlab


def generate_AIFs():
    ''' Generates AIFs from a Tracer distr. model, assuming
        different injection curves. Output is written to file
        AIFs_from_TDM.npy '''

    ## Choose/Load time array. Injection curve is sensitive to resolution. TO BE FIXED!
    time = np.arange(0,580.9,.1) # This needs to be exactly the same as used in freq_fit. why?!?
    # time = np.load('aif_time.npy')

    ## Load TDM parameter (from fit in 'freq_fit.py')
    p_TDM = np.load('p1_TDM.npy')


    ## Set Injection curve parameters
    x_max_list = np.arange(10,300,30) # original x_max = 21.4
    slope = .5
    volume = 20
    idx = 0

    ## Prepare array to write AIFs
    AIFs_from_TDM = np.zeros( [len(time), len(x_max_list)] )

    ## Itearte injection curves and create AIFs
    for x_max in x_max_list:
        injection = tl.Trapezoid(time, [ volume, x_max, slope] )

        AIF = tl.AIF_from_TDM( injection , time,  p_TDM)

        plt.figure(1)
        plt.plot(time,injection,'b')
        # plt.figure(2)
        plt.plot(time,AIF,'k')

        AIFs_from_TDM[:,idx] = AIF

        idx += 1
    plt.show()


    ## Save array with generated AIFs (1st index: time, 2nd index: AIF count)
    np.save('AIFs_from_TDM', AIFs_from_TDM)
    np.save('time_AIFs_from_TDM', time)
    np.save('x_max_list', x_max_list)



def determine_IRF():
    ''' Derives IRF parameters for a given concentration curve and AIF.
        The same IRF is used throughout the whole bootstrap procedure.
        Concentration curves are then simulated by using:
        concentration = tl.tofts_model(IRF_parameters, time_array, aif)'''

    ## Load concentration arrays (created in 'derive_T1.py')
    c_t = np.load('c_t.npy')
    c_t_time = np.load('time.npy')

    ## Load Parker AIF
    p_aif = np.load('p_aif.npy')
    parker_aif = tl.aif_fit_function(c_t_time, p_aif)

    ## Load AIF from linear fit in freq domain
    aif = np.load('AIFs_from_TDM.npy')[:,1]
    aif_time = np.load('time_AIFs_from_TDM.npy')

    ## interpolate aif to concentr. time
    aif = np.interp(c_t_time, aif_time, aif)
    plt.plot(parker_aif,'b')
    plt.plot(aif,'k')
    plt.show()

    ## Choose concentration curve
    c_t = c_t[15:18,12:15,0:1,:] # select ROI
    c_t = tl.average_mri_4d(c_t, ave_len = 3) # average
    c_t = c_t[0,0,0,:]

    p0 = np.array((0.00329732, 0.42715405))

    p_c_t, sucess, chi2_temp = tl.fit(\
        aif = aif,\
        targetfunc = c_t,\
        time = c_t_time,\
        model = 'tofts_model',\
        initial_parms = p0)


    print(sucess, p_c_t)
    plt.plot(c_t_time, c_t, 'x')
    plt.plot(c_t_time, tl.tofts_model(p_c_t, c_t_time, aif),'.')
    plt.show()

    np.save('tofts_concentration_parms', p_c_t)

#    plt.plot(c_t_time, aif)
#    plt.show()



def bootstrap():
    ''' Bootstrap procedure for different AIFs and Tofts Model.
        Loads aifs_from_TDM.npy, tofts_concentration_parms.npy.'''

    rand_iterations = 5000
    noise_amplitude = .2

    ## Load AIF from linear fit in freq domain and time
    aifs = np.load('AIFs_from_TDM.npy')

    ## Create time array (crucial for program performance)
    aif_time_original = np.load('time_AIFs_from_TDM.npy')
    aif_time = np.arange(500)
    # aif_time = aif_time_original

    ## Interpolate AIF accordingly
    aifs_interpolated = np.zeros((len(aif_time), len(aifs[1,:])))
    for k in range(len(aifs[0,:])):
        aifs_interpolated[:,k] = np.interp(aif_time, aif_time_original, aifs[:,k])
    aifs = aifs_interpolated


    ## Load concentration parms
    p_c_t = np.load('tofts_concentration_parms.npy')


    ## Prepare array to save fit parameters
    parms_distribution = np.zeros(( len(aifs[0,:]), rand_iterations ))

    ## Iterate AIFs
    for i in range(len(aifs[0,:])):
        print(i)
        aif = aifs[:,i]

        c_t = tl.tofts_model(p_c_t, aif_time, aif) #(1)

        ## Iterate noise
        for rand_index in np.arange(rand_iterations):
            #print(rand_index)
            c_t_temp = c_t[:] + (scipy.rand(len(c_t))-0.5)*\
                                   (noise_amplitude)

            #plt.plot(aif_time, c_t_temp)
            #plt.show()
            #plt.plot(aif_time, aifs[:,i])
            #plt.show()

            p0 = np.array((0.003, 0.5))
            #p0 = np.array((0, 0))

            p_fit, sucess, chi2_temp = tl.fit(\
                        aif = aif,\
                        targetfunc = c_t_temp,\
                        time = aif_time,\
                        model = 'tofts_model',\
                        initial_parms = p0)

            ## plot fits
            #plt.plot(aif_time, tl.tofts_model(p_fit, aif_time, aif))
            #plt.plot(aif_time, c_t_temp, 'x')
            #plt.show()

            ## Write fit parameters in array
            parms_distribution[i, rand_index] = p_fit[0]

    ## Examine distributions
    #plt.hist(parms_distribution[0,:])
    #plt.show()

    np.save('parms_distribution', parms_distribution)



def fit_parameter_distributions():
    ''' Fit Gaussians to parameter distributions.
        Loads parms_distribution.npy, x_max_list.npy '''

    ## Load parameters array [aif_count, rand_index]
    parms_distribution = np.load('parms_distribution.npy')
    #rand_iterations = parms_distribution.shape[1]
    aif_realizations = parms_distribution.shape[0]


    ## Prepare arrays to write gaussian parameters
    mu = []
    sigma = []
    bins = []

    plt.figure(2)
    for i in range(aif_realizations):

        (mu_temp, sigma_temp) = norm.fit(parms_distribution[i, :])
        hist, bins_temp, patch = plt.hist(parms_distribution[i, :],
                                normed=1, histtype = 'step', color = 'black', bins = 30 )

        mu.append(mu_temp)
        sigma.append(sigma_temp)
        bins.append(bins_temp)

    for i in range(aif_realizations):
        y = mlab.normpdf( bins[i], mu[i], sigma[i] )
        plt.plot(bins[i], y, linewidth=2 )

    sigma_rel = [(1000*a/b) for a, b in zip(sigma, mu)]
    np.save('sigma_rel', sigma_rel)

    plt.figure(2)
    plt.tick_params(axis='y', which='both', labelleft = 'off')
    plt.title(r'$K_{trans}$ distributions')
    plt.xlabel(r'$K^{trans}$ [1/s]')
    #plt.savefig('/home/tammo/Master/python/figures_nov_11/gaussians.eps', format = 'eps')
    plt.show()



def plot_bootstrap_results():
    ''' Loads sigma_rel, the relative standard deviations from bootstrapping. '''

    sigma_rel = np.load('sigma_rel.npy')

    plt.figure(3)
    #plt.plot(aif_factors, mu)
    plt.plot((np.load('x_max_list.npy')), sigma_rel, 'kx')
    #plt.legend('mu', 'sigma')
    #plt.title(r'$K_{trans}$ uncertainty as a function of AIF width')
    plt.xlabel(r'injection duration')
    plt.ylabel(r'relative standard deviation [permil]')
    plt.show()




def main():

    generate_AIFs()

#    determine_IRF()

#    bootstrap()

#    fit_parameter_distributions()

    plot_bootstrap_results()




if __name__ == '__main__':
    main()