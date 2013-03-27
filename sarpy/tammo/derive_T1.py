# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 20:31:10 2012
@author: tammo

This script derives a T1 matrix from time dependent preinjection scans
at a single flip angle. Output is written to 'T1_file.npy'.
"""

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import tammo_lib as tl
from scipy import optimize


def derive_T1():
    # Set path and filename for Matlab data file
    path = '/home/tammo/Master/python/Firas/'
    file = 'HerS7Ss01.fy1-DCE.mat'

    # Import post injection data with dimension (Y, X, SLICE, TIME (Z) )
    signal = scipy.io.loadmat(path+file)
    signal = signal['HerS7Ss01DCE']

    # ==========================================================================
    # ==== Derive signal over baseline signal and moving average in time =======
    # ==== Derive T1 and concentration arrays ==================================
    # ==========================================================================
    t_dump = 10 # number of data points at the beginning to be deleted
    c_t_time_spacing = 3.42 # firas thinks its 2.1 (but told tammo 3.42 before?)
    #    c_t_time_spacing = .01  # FOR AIF SIMULATION IN FREQ DOMAIN

    #signal = signal[66:67,53:54,0:1,:]     #SET ROI (50:100; 40:64)
    signal = signal[50:100, 40:64, 0:1, :]

    #    signal = tl.average_mri_4d(signal, ave_len=3)   #Average pixels

    #signal = tl.filter_noise(signal, bins=10, plots='None',\
    #            threshold = 1.10, data_start=0, data_end=100) # filter noise

    signal = tl.divide_by_baseline( signal, baseline_range = 10) # divide by baseline

    #signal = tl.timeline_moving_average (signal, ave_window = 1) # moving average for timelines

    T1_0 = np.zeros(signal.shape)
    T1_0[:] = 2000  # Set T1_0 values to 2000ms

    # Derive T1 and concentration time curve from signal
    T1, c_t = tl.T1_from_FLASH( signal, T1_0, TR = 400, alpha = 40, r = 4.3e-3 )

    # Cut first data points (time)
    T1, c_t = T1[:, :, :, t_dump:], c_t[:, :, :, t_dump:]

    # Create time array
    time = np.arange(0, len(c_t[0,0,0,:])*c_t_time_spacing-1, c_t_time_spacing)
    # ==========================================================================
    return time, T1, c_t


def import_AIF(AIFpath):
# ==========================================================================
# ======================= Import AIF =======================================
# ==========================================================================

    aif_time_spacing = 0.1
    aif_dump = 190 # number of data points to skip at beginning of AIF

    aif, aif_time = tl.import_AIF(path = AIFpath,\
                                    aif_time_spacing = aif_time_spacing,\
                                    aif_dump = aif_dump)

    plt.plot(aif_time, aif)
    plt.xlabel('Time [s]')
    plt.ylabel('Concentration [some units]')
    plt.show()
    # =========================================================================


# ==========================================================================
# ================== Fit parker AIF-Function to AIF ========================
# ==========================================================================
    def aif_err_func(p, time, aif):
        err = aif - tl.aif_fit_function( time, parms = p)
        return err

    p0 = [ -1e+05, -1e+02, -3e+02, -7e+00, 2e+00, -4e-02,
            8e-01, 3e-04, 1.5e+00, 7e+00]


    p_aif, sucess = optimize.leastsq(aif_err_func, p0, args=(aif_time, aif))
    print(sucess)

    np.save('p_aif', p_aif)

    p_aif = np.load('p_aif.npy')
    # Create according AIF
    parker_aif = tl.aif_fit_function(aif_time, p_aif)

    # Plot AIF fit
    #plt.plot(aif_time, aif,'x')
    #plt.plot(aif_time, tl.aif_fit_function(aif_time, p_aif),'-')
    #plt.show()
    # ========================================================================
    return aif_time, aif, p_aif, parker_aif




def main():
    ''' Se functions for description. '''
#    time, T1, c_t = derive_T1()

#    print(len(time), len(c_t))
    short_aif_time, short_aif, short_parker_aif_parms, short_parker_aif = import_AIF(\
        AIFpath = '/home/tammo/Master/python/AIFs/aif_old/popaveAIF.csv')

    #Interpolate aif to c(t)-time
    #aif = np.interp(time, aif_time, aif)

#    np.save('c_t', c_t)
#    np.save('T1', T1)
#    np.save('time', time)
#
#    np.save('aif_time', aif_time)
#    np.save('aif', aif)
#
#    np.save('p_aif', p_aif)
#    np.save('parker_aif', parker_aif)

    return 0


if __name__ == '__main__':
	main()





#========================== CREATE VIDEO ===================================
# for j in range(len(data[1,1,:,1])):
#     for i in range(len(data[1,1,1,:])):
#         plt.imshow(data[50:80,40:64,j,i],interpolation='none')
#         plt.title('slice'+str(j)+'@'+str(i)+'s')
#         plt.savefig('animation'+'%0.3d' %(i)+'.png',format='png')
#
#     os.system('convert -delay 20 -loop 0 animation*.png animation'+str(j)+'.gif')
#===========================================================================

#pdb.set_trace()
#print ("Ich bin eine Fehlernachricht", file=sys.stderr)
