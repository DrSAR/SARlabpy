# -*- coding: utf-8 -*-
# Copyright (C) 2012-2013 Andrew Yung and UBC 7T MRI lab
"""
Analysis routines for Bloch-Siegert B1 and flip angle mapping on Bruker data. 

Ref:  B1 mapping by Bloch-Siegert shift. Sacolick LI, Wiesinger F, Hancu I, 
Vogel MW.  Magn Reson Med. 2010 May;63(5):1315-22. 

Brief description of method:  off-resonance pulse applied after spin excitation, far away from Larmor frequency that it does not perturb magnetization, but strong enough to perturb direction of magnetic field.  Resultant phase shift is related to the peak B1 of the off-resonance pulse.  This B1 measurement can be used to map B1 and flip angle of other pulses (with knowledge of relative transmit power levels).

The technique requires at least two scans, repeated at positive and negative offset centered around the Larmor frequency.  The phase difference is proportional to B1, which is independent of B0 inhomogeneity up to the first order.  The method requires an additonal two scans with zero-power off-resonant pulses to determine the phase shift originating from switching the frequency synthesizer for the Bloch-Siegert pulse.  As of April 2013, these 4 scans are encapsulated into one expno as separate scan repetitions if BSB1mapFLASH is used.

=====================

Assumptions:
    
1) The user is interested in using the peak B1 from the Bloch-Siegert experiment to determine the peak B1 or flip angle of a RF pulse of interest ("POI" from a different scan).  

2) The experimenter has set appropriate Bloch-Siegert offsets and powers which produce phase difference maps which have at most one phase wrap.  The Bloch-Siegert pulse parameters therefore depend on the choice of transmit coil and its B1 inhomogeneity. 

3) The Bloch-Siegert input data is a phase image reconstructed from Paravision.

4) The pulse shapes are available in ~/wave (use symbolic links if they are somewhere else).
   
5) Sample data used in test functions are available in ~/data/<anyusername>/nmr

=====================

Several function wrappers exist for calculation of flip angle maps of a pulse of interest (POI) which differ by their input type:
    
1) BS_flipangle_wrt_Bruker2dseq: file path of Bloch-Siegert expno and POI expno 

2) BS_flipangle_wrt_ScanObject: sarlabpy Scan objects

3) calc_BS_flipangle:  arrays for the individual Bloch-Siegert acquired scans

write_BS_flipangle_2dseq creates a 2dseq and visu_pars file for the calculated flip angle map, thus enabling direct visualization of the results in Paravision image display. 

"""

from __future__ import division

import numpy
import os
import pylab
import sarpy
import sarpy.utils.plt_utils
import scipy.integrate
import scipy.interpolate
import scipy.optimize
import scipy.fftpack
import math
import re

debug=False



def BS_test(studyname,BS_expno,BS_procno,POI_expno,offConsole=True):
    '''
    Test function for Bloch-Siegert data with visualization of flip angle map
    
    :param string studyname:  directory name of study
    :param int BS_expno:  expno of Bloch-Siegert input data
    :param int BS_procno:  procno of Bloch-Siegert input data phase
    :param int POI_expno:  expno which contains the RF pulse of interest (can be acquired or just ready to scan)
    :param bool offConsole:  is the program being run on or off the scanner computer?  (affects our choice of path for RF pulse shapes)
    '''
    
    BS_path = os.path.join(studyname,str(BS_expno))
    POI_path = os.path.join(studyname,str(POI_expno))

    if offConsole:
        shapefile_path = os.path.expanduser(os.path.join('~','wave'))
    else:
        shapefile_path = 'opt/PV5.1/exp/stan/nmr/lists/wave'

    B1peakmap, flipanglemap = BS_B1_flipangle_wrt_BrukerPath(BS_path,BS_procno,POI_path,shapefile_path)
    
    #visualize the flip angle of POI in the middle slice of the BS slice pack
    i_slice = int(math.floor(flipanglemap.shape[2]/2))
    pylab.imshow(flipanglemap[:,:,i_slice],interpolation='none')
    pylab.title('flip angle map of pulse inside ' + POI_path)
    pylab.colorbar()
    pylab.show()


    
def BS_B1_flipangle_wrt_BrukerPath(BS_path,BS_procno,POI_path,shapefile_path,abspath=False):
    '''
    Return the peak B1 and flip angle map for a RF pulse of interest, given the expno paths of 
    
    :param str BS_path: expno path of Bloch-Siegert scan (all data in one expno)
    :param int BS_procno: procno number with Bloch-Siegert phase data
    :param str POI_path: expno path containing RF pulse of interest (POI)
    :param str shapefile_path: path of RF shape files
    :param bool abspath: 
        if False: paths given as 'studydirname/expno' and root path is ~/data/<anyusername)/nmr
        if True: paths given as absolute paths
    :return:
        image sets of peak B1 and flip angle for pulse of interest
    '''
    
    BS_scan = sarpy.io.BRUKER_classes.Scan(BS_path, abspath)
    POI_scan = sarpy.io.BRUKER_classes.Scan(POI_path, abspath)
    
    B1peak, flipangle = BS_B1_flipangle_wrt_ScanObject(BS_scan, BS_procno, 
                                                       POI_scan, shapefile_path)

    return B1peak, flipangle


    
def BS_B1_flipangle_wrt_ScanObject(BS_scan, BS_procno, POI_scan, shapefile_path):
    '''
    Return peak B1 and flipangle for a RF pulse of interest, given some sarpy.io.BRUKER_classes Scan objects

    :param Scan BS_scan: Scan object containing all Bloch-Siegert data
    :param int BS_procno: procno number with Bloch-Siegert phase data
    :param str POI_scan: Scan object containing RF pulse of interest (POI)
    :param str shapefile_path: path of RF shape files
    :return:
        image sets of peak B1 and flip angle for pulse of interest
    '''
    
    BS_method = BS_scan.method.Method
    assert BS_method == 'bSB1mapFLASH' or BS_method == 'bSB1mapMSME', \
                        'wrong scan type'    
     
    try:
        freqoffset = BS_scan.method.BSFreqOffset
        dBpwr_BS = BS_scan.method.BSPulse[3]
        width_BS = BS_scan.method.BSPulse[0]*1e-3
        BS_shape = BS_scan.method.BSPulse[8]
        BS_shape = \
            re.findall(r'<\s*(\S+)\s*>',BS_shape)[0]  #regex to strip off <>      

        POI_parx =  {}
        POI_parx['dBpwr'] = POI_scan.method.ExcPulse[3]
        POI_parx['integralratio'] = POI_scan.method.ExcPulse[10]
        POI_parx['width'] = POI_scan.method.ExcPulse[0]*1e-3
        POI_parx['name'] = \
                re.findall(r'<\s*(\S+)\s*>',POI_scan.method.ExcPulse[8])[0]

        acq_ndim = POI_scan.acqp.ACQ_dim
        if acq_ndim == 3:
            POI_parx['nPE2'] = POI_scan.acqp.ACQ_size[2]
        
    except:
        print('something wrong happened with pulse parameter query')
        
    try:
        assert BS_scan.pdata[BS_procno-1].reco.RECO_image_type == 'PHASE_IMAGE', \
                        'BS data is not a phase image'
        on_BSplus_phase = BS_scan.pdata[BS_procno-1].data[:,:,:,0]
        on_BSminus_phase = BS_scan.pdata[BS_procno-1].data[:,:,:,1]
        off_BSplus_phase = BS_scan.pdata[BS_procno-1].data[:,:,:,2]
        off_BSminus_phase = BS_scan.pdata[BS_procno-1].data[:,:,:,3]
    except:
        print('error in retrieving phase images')

    if BS_method=='bSB1mapFLASH':
        seqtype = 'gradientecho'
    elif BS_method=='bSB1mapMSME':
        seqtype = 'spinecho'
        
    KBS = calc_KBS(freqoffset, width_BS, BS_shape, shapefile_path, seqtype)

    B1peak_BS = calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, 
                               off_BSminus_phase, off_BSplus_phase)

    B1peak_POI = calc_POI_B1peak(B1peak_BS, dBpwr_BS, POI_parx['dBpwr'])
    
    flipangle_POI = calc_POI_flipangle(B1peak_POI, POI_parx, shapefile_path,
                                       on_resonance = False)
        
    return B1peak_POI, flipangle_POI



def calc_POI_B1peak(B1peak_BS, dBpwr_BS, dBpwr_POI):
    '''
    Return peak B1 of pulse of interest in Tesla, based on Bloch-Siegert peak B1 and 
    pulse powers 
    '''    
    B1peak_POI = B1peak_BS * (math.pow(10,(dBpwr_BS-dBpwr_POI)/20))
    return B1peak_POI
 
 
  
def calc_POI_flipangle(B1peak_POI, POI_parx, shapefile_path,
                          on_resonance=False, acq_ndim=2):
    '''
    return flip angle of pulse of interest in degrees
    
    :param array B1peak_POI: peak B1 of pulse of interest (POI)
    :param dict POI_parx: contains relevant parameters of POI
    :param str shapefile_path: path where RF shape files located
    :param bool on_resonance:  are spins assumed to be on resonance?
    :param int acq_ndim: number of dimensions in the acquisition (2D or 3D)
    :return:
        flip angle map of pulse of interest
    '''
    
    gamma = 267.513e6 # radians/(s*T)
    integ_ratio = POI_parx['integralratio']
    pulse_width = POI_parx['width']
    if on_resonance:
        '''
        flip angle calc for on_resonant spins: gamma*(time integral of B1)
        '''
        flipangle_POI = gamma*B1peak_POI*integ_ratio*pulse_width
    else:
        '''
        flip angle calc which accounts for off-resonant spins due to slice
        select gradient (assumes that slice profile = FT of gamma*|B1(t)|,
        which is valid at small flip angles:  
        '''
        pulse = \
            sarpy.io.BRUKERIO.readRFshape(os.path.join(shapefile_path, POI_parx['name']))
        B1_amp = numpy.asarray(pulse['amp'])
        # take Fourier transform of B1(t)
        pulseprofile = abs(scipy.fftpack.fftshift(scipy.fftpack.fft(B1_amp)))
        n_profile = pulseprofile.size
        dt = pulse_width/pulse['NPOINTS']
        df = 1/(dt*n_profile)
        dw = 2*math.pi*df
        # interpolate pulse profile to improve FWHM positions and integral across slice        
        interp_factor = 10
        profile_eqn = \
             scipy.interpolate.interp1d(numpy.linspace(0,(n_profile-1)*dw,n_profile),
                                        pulseprofile,kind='quadratic')
        interp_profile = profile_eqn(numpy.linspace(0,(n_profile-1)*dw,interp_factor*n_profile))
        # normalize pulse amplitude:  centre of slice (zero freq offset) has flip angle = gamma*(time integral of B1(t))
        interp_profile = interp_profile/max(interp_profile)*gamma*integ_ratio*pulse_width
        # truncate pulse to its FWHM        
        lo, hi = find_FWHM_limits(interp_profile)
        BW = dw/interp_factor*(hi-lo)
        trunc_interp_profile = interp_profile[lo:hi]
        n_trunc_profile = trunc_interp_profile.shape[0]
  
        if acq_ndim == 2:
            '''
            2D excitation:  
                flip angle for a slice = (integral of flip angle across slice)/BW
            '''
            integrand = trunc_interp_profile
            flipangle_integral = \
                scipy.integrate.trapz(integrand, None, dw/interp_factor)
            flipangle_POI = flipangle_integral*B1peak_POI/BW
        elif acq_ndim == 3:
            '''
            3D excitation:  
                flip angle for a PE2 slice = gamma*(integral of flip angle across PE2 section)/(width of PE2 section)
            '''
            nPE2 = POI_parx['nPE2']
            n_interval = math.floor(n_trunc_profile/nPE2)
            interval = BW/nPE2
            
            for k in xrange(B1peak_POI.shape[2]):
                integrand = \
                    trunc_interp_profile[k*n_interval:(k+1)*n_interval]
                flipangle_integral = \
                    scipy.integrate.trapz(integrand, None, dw/interp_factor)
                flipangle_POI[:,:,k] = flipangle_integral*B1peak_POI[:,:,k]/interval
            
    return flipangle_POI*180/math.pi
 
 
 
def find_FWHM_limits(pulse):
    '''
    return the indices where the FWHM of a pulse is found
    
    :param array pulse: 1D array of pulse shape
    :returns:
        lo and hi indices of FWHM points
    '''
    n_points = pulse.shape[0]
    centre_index = pulse.argmax()
    centre_ampl = pulse.max()
    halfmax_ampl = centre_ampl/2
    index = centre_index
    ampl = centre_ampl

    # start from centre of peak and look forward
    while ampl > halfmax_ampl and index < n_points:
        index += 1
        ampl = pulse[index]

    assert index < n_points, 'half max not found wrt pulse centre'

    hi = index
    halfwidth = index - centre_index
    lo = centre_index - halfwidth
        
    return lo, hi


   
def calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, off_BSminus_phase, off_BSplus_phase):
    '''
    calculates B1peak of the Bloch-Siegert pulse, based on phase differences  

    :param float KBS: Bloch-Siegert calibration constant in radians/T^2
    :param array on_BSminus_phase: imageset with BS pulse on (-ve offset)
    :param array on_BSplus_phase: imageset with BS pulse on (+ve offset)
    :param array off_BSminus_phase: imageset with zero flip angle BS pulse (-ve offset)
    :param array off_BSplus_phase: imageset with zero flip angle BS pulse (+ve offset)
    :return:
        imageset of B1peak values of the Bloch-Siegert pulse
    :rtype: array
    '''

    pi = math.pi
    offset = off_BSplus_phase - off_BSminus_phase
    offset[offset>pi] = offset[offset>pi] - 2*pi

 
    phase_diff = on_BSplus_phase - on_BSminus_phase - offset
    phase_diff[phase_diff>pi] = phase_diff[phase_diff>pi] - 2*pi
    phase_diff[phase_diff<-pi] = phase_diff[phase_diff<-pi] + 2*pi

    if debug:
        i_slice = 12        
        #fig1=pylab.imshow(180/pi*offset[:,:,i_slice],interpolation='None')
        pylab.imshow(180/pi*phase_diff[:,:,i_slice],interpolation='None')
        pylab.colorbar()
        pylab.show()

    B1peak_BS = numpy.sqrt(numpy.absolute(phase_diff)/(2*KBS))

    return B1peak_BS
    
   

def calc_KBS(freqoffset, pulsewidth, pulseshape, shapefile_path, mode='gradientecho'):
    '''
    calculates KBS in the units of radians/T^2 (assumes freq offset is constant in time).  

    :param float freqoffset: BS offset in Hz
    :param float pulsewidth: length of pulse in seconds
    :param string pulseshape: filename of pulse shape
    :param string shapefile_path: directory of pulse_shape
    :param string mode: 'gradientecho' or 'spinecho'
    :return:
        value of KBS in radians/T^2
    :rtype: float
    
    Divide the result by 1e8 if units of radians/Gauss^2.  See literature reference in the module docstring for details.

Example:
    
    >>> shapefile_path = os.path.expanduser(os.path.join('~','wave'))
    >>> KBS = calc_KBS(4000,8e-3,'fermi.exc',shapefile_path)
    >>> KBS
    7109603109.3280497
    '''
    filename = os.path.join(shapefile_path, pulseshape)
    pulse = sarpy.io.BRUKERIO.readRFshape(filename)
    dt = pulsewidth/pulse['NPOINTS']
    amp = numpy.asarray(pulse['amp'])
    amp = amp/max(amp)

    gamma = 267.513e6 # radians/(s*T)
    integrand = (gamma*amp)**2/(2*2*3.141592654*freqoffset)
    KBS = scipy.integrate.trapz(integrand, None, dt)

    assert mode == 'gradientecho' or mode == 'spinecho', (
        'mode is not gradientecho or spinecho')
    
    if mode == 'spinecho':
        return 2*KBS
    else:
        return KBS


    
def BS_test_originaldata():
    '''
    Test function using Bloch-Siegert data contained in 4 scans (assumes phase
    is reconstructed in the same procno across all 4 scans)
    
    Shows the flip angle map of the pulse of interest    
    '''
    shapefile_path = os.path.expanduser(os.path.join('~','wave'))
   
    studypath = 'dBlochSiegert1.gP2'
    BS_off_neg = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'36'))
    BS_off_pos = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'37'))
    BS_on_neg = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'46'))
    BS_on_pos = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'47'))
    
    POI_scan = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'32'))

    B1peak, flipangle = BS_B1_flipangle_wrt_4ScanObjects(BS_on_pos,BS_on_neg,BS_off_pos,BS_off_neg, 2, POI_scan, shapefile_path)
    
    
    i_slice = int(BS_on_pos.acqp.NSLICES/2)
    pylab.imshow(flipangle[:,:,i_slice],interpolation='none')
    pylab.title('flip angle map of pulse inside ' + POI_scan.shortdirname)
    pylab.colorbar()
    pylab.show()            



def BS_B1_flipangle_wrt_4ScanObjects(Scan_pos_on, Scan_neg_on, Scan_pos_off, Scan_neg_off, BS_procno, POI_scan, shapefile_path):
    '''
    Return peak B1 and flipangle for a RF pulse of interest, given some sarpy.io.BRUKER_classes Scan objects (4 scans for BS data, 1 scan for pulse of interest)

    :param Scan Scan_pos_on: Scan object with +ve BS pulse
    :param Scan Scan_neg_on: Scan object with -ve BS pulse
    :param Scan Scan_pos_off: Scan object with zero-power +ve BS pulse
    :param Scan Scan_neg_off: Scan object with zero-power -ve BS pulse
    :param int BS_procno: procno number with Bloch-Siegert phase data
    :param str POI_scan: Scan object containing RF pulse of interest (POI)
    :param str shapefile_path: path of RF shape files
    :return:
        image sets of peak B1 and flip angle for pulse of interest
    '''
    import re
    
    BS_method = Scan_pos_on.method.Method
    assert BS_method == 'bSB1mapFLASH' or BS_method == 'bSB1mapMSME', 'wrong scan type'    
    
    try:
        dBpwr_BS = Scan_pos_on.method.BSPulse[3]
        width_BS = Scan_pos_on.method.BSPulse[0]*1e-3
        freqoffset = Scan_pos_on.method.BSFreqOffset
        BS_shape = Scan_pos_on.method.BSPulse[8]
        BS_shape = re.findall(r'<\s*(\S+)\s*>',BS_shape)[0]  #regex to strip off <>      

        POI_parx =  {}
        POI_parx['dBpwr'] = POI_scan.method.ExcPulse[3]
        POI_parx['integralratio'] = POI_scan.method.ExcPulse[10]
        POI_parx['width'] = POI_scan.method.ExcPulse[0]*1e-3
        POI_parx['shape'] = \
                re.findall(r'<\s*(\S+)\s*>',POI_scan.method.ExcPulse[8])[0]
    except:
        print('something wrong happened with pulse parameter query')
        
    try:
        on_BSplus_phase = Scan_pos_on.pdata[BS_procno-1].data
        on_BSminus_phase = Scan_neg_on.pdata[BS_procno-1].data
        off_BSplus_phase = Scan_pos_off.pdata[BS_procno-1].data
        off_BSminus_phase = Scan_neg_off.pdata[BS_procno-1].data
    except:
        print('error in retrieving phase images')

    if BS_method=='bSB1mapFLASH':
        seqtype = 'gradientecho'
    elif BS_method=='bSB1mapMSME':
        seqtype = 'spinecho'
        
    KBS = calc_KBS(freqoffset, width_BS, BS_shape, shapefile_path, seqtype)

    B1peak_BS = calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, off_BSminus_phase, off_BSplus_phase)
    B1peak_POI = calc_POI_B1peak(B1peak_BS, dBpwr_BS, POI_parx['dBpwr'])
    
    flipangle_POI = calc_POI_flipangle(B1peak_POI, POI_parx, shapefile_path,
                                        on_resonance=True)        
    return B1peak_POI, flipangle_POI


   
if __name__ == '__main__':
    studyname = 'dBlochSiegert2.j41'
    BS_expno = 1
    BS_procno = 2
    BS_test('dBlochSiegert2.j41',1,2,2)
#    BS_test_originaldata()
    import doctest
    doctest.testmod()






