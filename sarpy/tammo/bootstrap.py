# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 14:41:27 2012

Bootstrip method to derive parameter distribution for various AIFs

@author: tammo
"""

from __future__ import print_function
from __future__ import division

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['CM']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

import tammo_lib as tl
from scipy import optimize
import scipy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import norm
import matplotlib.mlab as mlab

matplotlib.rcParams.update({'font.size': 22})

##### Determine IRF parameters ###########################################
aif = np.load('parker_aif.npy')

c_t = np.load('c_t.npy')
c_t_time = np.load('time.npy')


c_t = c_t[15:18,12:15,0:1,:] # select ROI
c_t = tl.average_mri_4d(c_t, ave_len = 3) # average

#c_t = c_t[16:17,13:14,0:1,:]
c_t = c_t[0,0,0,:]


p0 = np.array((0.01,0.5))

p_c_t, sucess, chi2_temp = tl.fit(\
    aif = aif,\
    targetfunc = c_t,\
    time = c_t_time,\
    model = 'tofts_model',\
    initial_parms = p0)

#plt.plot(c_t_time, c_t, 'x')
#plt.plot(c_t_time, tl.tofts_model(p_c_t, c_t_time, aif))
#plt.show()
#
#plt.plot(c_t_time, aif)
#plt.show()
##########################################################################

#p_aif =  [1,100,1,.001]
#plt.plot(c_t_time, tl.aif_fit_tammo(c_t_time, p_aif))
#plt.show()





###### Determine AIF parameters ###########################################

def aif_err_func(p, time, aif):
    err = aif - tl.aif_fit_function(time, parms = p)
    return err

p0 = [ -1.08e+05,  -9.80e+01,  -3.25e+02,  -6.95e+00,
      2.12e+00, -3.94e-02,   7.81e-01,   2.92e-04,
      1.51e+00,   7.18e+00]

#p0 = [1,100,1,.001]

p_aif, sucess = optimize.leastsq(aif_err_func, p0, args=(c_t_time, aif))
print(p_aif)
##########################################################################

delta_time = 1

time = np.arange(0, 1400, 1)

#plt.plot(c_t_time, aif, 'x')
#plt.plot(c_t_time, tl.aif_fit_function(c_t_time, p_aif))
#plt.show()

aif_factors = np.arange(0.1,5.000001,0.1)

#for aif_factor in aif_factors:
#    aif = tl.scale_aif(aif_factor, p_aif, time)
#
#    area = sum(aif)
#    print(area)
#    plt.plot(time, aif)
#
#plt.show()
#
#
#m SystemExit

############################# Iterate Noise ################################

print(5*np.mean(c_t)/np.nanmax(c_t))

rand_iterations = 3000
noise_amp = np.nanmax(c_t)*0.2
K_trans = np.zeros( [aif_factors.shape[0], rand_iterations] )


aif_index = 0
for aif_factor in aif_factors:
    print(aif_factor)

    # aif normalized to maximum of resulting concentration
    # aif = tl.scale_aif(aif_factor, p_aif, time) / np.max(tl.tofts_model(p_c_t, time, tl.scale_aif(aif_factor, p_aif, time)))

    aif = tl.scale_aif(aif_factor, p_aif, time)


    plt.figure(0)
    plt.plot(time, aif)
    plt.xlabel('time [s]')
    plt.ylabel('concentration [mMol]')
#    plt.title('AIFs')

    c_t = tl.tofts_model(p_c_t, c_t_time, aif)
#    c_t /= max(c_t)
    plt.figure(1)
    plt.plot(time, c_t)
    plt.xlabel('time [s]')
    plt.ylabel('concentration [mMol]')
#    plt.title('contrast agent concentration in tissue')

    for rand_index in np.arange(rand_iterations):

        c_t_temp = c_t[:] + (scipy.rand(len(c_t))-0.5)*\
                                   (noise_amp)

        p0 = np.array((0.01,0.5))

        p1, sucess, chi2_temp = tl.fit(\
            aif = aif,\
            targetfunc = c_t_temp,\
            time = c_t_time,\
            model = 'tofts_model',\
            initial_parms = p0)


#        plt.figure(1)
#        plt.plot(c_t_time, c_t_temp)
#        plt.plot(c_t_time, tl.tofts_model(p1, c_t_time, aif))
#        plt.show()
        K_trans[aif_index, rand_index] = (p1[0]*60)

    aif_index += 1


plt.figure(0)
plt.savefig('/home/tammo/Master/python/figures_nov_11/AIFs.eps', format = 'eps')


plt.figure(1)
plt.savefig('/home/tammo/Master/python/figures_nov_11/concentrations.eps', format = 'eps')


##############################


mu = []
sigma = []
bins = []

plt.figure(2)
for i in range(aif_factors.shape[0]):

    (mu_temp, sigma_temp) = norm.fit(K_trans[i, :])
    hist, bins_temp, patch = plt.hist(K_trans[i, :],
                            normed=1, histtype = 'step', color = 'black', bins = 30 )

    mu.append(mu_temp)
    sigma.append(sigma_temp)
    bins.append(bins_temp)

for i in range(aif_factors.shape[0]):
    y = mlab.normpdf( bins[i], mu[i], sigma[i] )
    plt.plot(bins[i], y, linewidth=2 )


plt.figure(2)
plt.tick_params(axis='y', which='both', labelleft = 'off')
#plt.title(r'$K_{trans}$ distributions')
plt.xlabel(r'$K^{trans}$ [1/min]')
plt.savefig('/home/tammo/Master/python/figures_nov_11/gaussians.eps', format = 'eps')

sigma_rel = [(1000*a/b) for a, b in zip(sigma, mu)]

np.save('mu_for_different_aifs', mu )
np.save('sigma_for_different_aifs', sigma )
np.save('sigma_rel_for_different_aifs', sigma_rel )
np.save('aif_factors', aif_factors)



plt.figure(3)
#plt.plot(aif_factors, mu)
plt.plot(aif_factors, sigma_rel, 'x')
#plt.legend('mu', 'sigma')
#plt.title(r'$K_{trans}$ uncertainty as a function of AIF width')
plt.xlabel(r'AIF peak width factor')
plt.ylabel(r'relative standard deviation [permil]')

plt.savefig('/home/tammo/Master/python/figures_nov_11/result.eps', format = 'eps')
#plt.show()


raise SystemExit
















#========
mu = []
sigma = []
sigma_rel = []
bins = []

plt.figure(figsize=(9,4))

for i in np.arange(K_trans_distr.shape[2]):

    (mu_temp, sigma_temp) = norm.fit(K_trans_distr[0,0,i,:])

    mu.append(mu_temp)
    sigma.append(sigma_temp)
    sigma_rel.append(sigma_temp/mu_temp)

    print('\n aif ', i, ': ', 'mu = ', mu_temp, 'sigma = ', sigma_temp,\
        'sigma_rel = ', sigma_temp/mu_temp)

    hist, bins_temp, bla = plt.hist(K_trans_distr[0,0,i,:],\
        normed=1, bins=20, histtype = 'step', color = 'black')
    bins.append(bins_temp)


for i in np.arange(len(bins)):
    y = mlab.normpdf( bins[i], mu[i], sigma[i] )
    l = plt.plot(bins[i], y, linewidth=2 )

plt.xlabel('K_trans [1/min]')
plt.ylabel('PD [a.u.]')
plt.show()
#========

