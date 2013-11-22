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
    
1) The user is interested in using the peak B1 from the Bloch-Siegert experiment to determine the peak B1 or flip angle of a RF pulse of interest ("POI" from a different scan).  We must assume that the geometry and acquisition mode (2D or 3D) of the POI scan matches the BS scan.

2) The Bloch-Siegert pulse parameters depend on the choice of transmit coil and its B1 inhomogeneity. The experimenter has set appropriate Bloch-Siegert offsets and powers which produce phase difference maps which have an in-plane variation of at most one phase wrap.  The code can add an integer number of phase wraps in the Bloch-Siegert B1 calculations, to accoun for situations where the phase difference has experienced several phase wraps across the whole volume.  

3) The Bloch-Siegert input data is a phase image reconstructed from Paravision.

4) The pulse shapes are available in ~/wave (use symbolic links if they are somewhere else).
   
5) Sample data used in test functions are available in ~/data/<anyusername>/nmr

=====================

Several function wrappers exist for calculation of flip angle maps of a pulse of interest (POI) which differ by their input type:
    
1) BS_flipangle_wrt_Bruker2dseq: file path of Bloch-Siegert expno and POI expno 

2) BS_flipangle_wrt_ScanObject: sarlabpy Scan object (1 or 4 Scan objects)

3) calc_POI_flipangle:  arrays for the individual Bloch-Siegert acquired scans

Test functions include BS_test and BS_test_originaldata (data acquired as 4 separate expnos)

TODO:
1) write_BS_flipangle_2dseq creates a 2dseq and visu_pars file for the calculated flip angle map, thus enabling direct visualization of the results in Paravision image display. 

2) make consistency checks for POI vs. BS scans 

"""

from __future__ import division

import numpy
import os
import pylab
import sarpy
import scipy.integrate
import scipy.interpolate
import scipy.optimize
import scipy.fftpack
import math
import re
import sarpy.io.BRUKERIO

debug=True

def write_BS_2dseq(B1map,flipanglemap,expno_path,procno):

    FGparx = [{},{}]    

    procno_path =  os.path.join(expno_path,'pdata',str(procno))
#    B1_flipangle_data = numpy.zeros(B1map.shape + (2L,))
#    B1_flipangle_data[:,:,:,0] = B1map
#    B1_flipangle_data[:,:,:,1] = flipanglemap
#    wordtype = 'int32'
#    frame_count = 2*B1map.shape[2]
#
#    
#
#
#    visu_slope, datamin, datamax, dataoffs, binarydata = \
#        sarpy.io.BRUKERIO.write2dseq(B1_flipangle_data,procno_path,
#                                     wordtype,frame_count)
#
#    FGparx[0]['grouplen'] = B1map.shape[2]
#    FGparx[0]['groupid'] = '<FG_SLICE>'
#    FGparx[0]['groupdesc'] = '<>' 
#    FGparx[0]['groupdepval'] = ['VisuCoreOrientation','VisuCorePosition']
#    FGparx[0]['elemid'] = ['']
#    FGparx[0]['elemcomments'] = ['']
#    FGparx[0]['elemunits'] = ['']
#
#    FGparx[1]['grouplen'] = 2
#    FGparx[1]['groupid'] = '<FG_BS>'
#    FGparx[1]['groupdesc'] = '<Bloch-Siegert maps>' 
#    FGparx[1]['groupdepval'] = ['VisuCoreDataUnits', 'VisuFGElemId',
#                                'VisuFGElemComment']
#    FGparx[1]['elemid'] = ['<BS_B1>','<BS_FLIP>']
#    FGparx[1]['elemcomments'] = ['<B1peak of POI>', '<eff flip angle of POI>']
#    FGparx[1]['elemunits'] = ['<T>','<degrees>']
#
#
#    sarpy.io.BRUKERIO.populate_visupars(procno_path,binarydata,
#                                        frame_count,FGparx,datamax,datamin,
#                                        dataoffs,visu_slope,wordtype,is2D=True)
                                        
    # test by loading original data and trying to write it
    BS_scan = sarpy.io.BRUKER_classes.Scan(expno_path, False)
    BS_data=BS_scan.pdata[3].data
    wordtype = 'int16'
    frame_count = BS_data.shape[2]*BS_data.shape[3]*BS_data.shape[4]
    BS_data = -BS_data
    visu_slope, datamin, datamax, dataoffs, binarydata = \
        sarpy.io.BRUKERIO.write2dseq(BS_data,procno_path,
                                     wordtype,frame_count)    



    return 1

def BS_test(studyname,BS_expno,BS_procno,POI_expno,uname,n_phasewraps=0,offConsole=True):
    '''
    Test function for Bloch-Siegert data with visualization of flip angle map
    
    Arguments:
        string studyname:  directory name of study
        int BS_expno:  expno of Bloch-Siegert input data
        int BS_procno:  procno of Bloch-Siegert input data phase
        int POI_expno:  expno of the RF pulse of interest
        int n_phasewraps:  number of phase wraps across whole volume
        bool offConsole:  is the program being run on or off the scanner computer?  (affects our choice of path for RF pulse shapes)
    '''
    uname_path = os.path.join(os.path.expanduser('~'),'data',uname,'nmr')
    BS_path = os.path.join(uname_path,studyname,str(BS_expno))
    POI_path = os.path.join(uname_path,studyname,str(POI_expno))

    if offConsole:
        shapefile_path = os.path.expanduser(os.path.join('~','wave'))
    else:
        shapefile_path = 'opt/PV5.1/exp/stan/nmr/lists/wave'

    B1peakmap, flipanglemap = BS_B1_flipangle_wrt_BrukerPath(BS_path,BS_procno,POI_path,shapefile_path,n_phasewraps)
    
#    BS_parmap_procno = 4
#    write_BS_2dseq(B1peakmap,flipanglemap,BS_path,BS_parmap_procno)    
    
    #visualize the flip angle of POI in the middle slice of the BS slice pack
    i_slice = int(math.floor(flipanglemap.shape[2]/2))
    pylab.imshow(flipanglemap[:,:,i_slice],interpolation='none')
    pylab.title('flip angle map of pulse inside ' + POI_path)
#    pylab.imshow(B1peakmap[:,:,i_slice],interpolation='none')
#    pylab.title('B1peak map of pulse inside ' + POI_path)
    pylab.colorbar()
    pylab.show()


    
def BS_B1_flipangle_wrt_BrukerPath(
 BS_path,BS_procno,POI_path,shapefile_path,n_phasewraps=0,abspath=False):
    '''
    Return peak B1 & flip angle map for pulse of interest (input = pathnames)
    
    Arguments:
        str BS_path: expno path of Bloch-Siegert scan (all data in one expno)
        int BS_procno: procno number with Bloch-Siegert phase data
        str POI_path: expno path containing RF pulse of interest (POI)
        int n_phasewraps:  number of phase wraps across whole volume
        str shapefile_path: path of RF shape files
        bool abspath: 
            if False: paths given as 'studydirname/expno' and root path is ~/data/<anyusername)/nmr
            if True: paths given as absolute paths
    Returns:
       array: images of peak B1 and flip angle for pulse of interest
    '''
    
    BS_scan = sarpy.io.BRUKER_classes.Scan(BS_path, abspath)
    POI_scan = sarpy.io.BRUKER_classes.Scan(POI_path, abspath)
    
    B1peak, flipangle = BS_B1_flipangle_wrt_ScanObject(POI_scan,shapefile_path,
                                                       n_phasewraps,BS_scan,
                                                       BS_procno)

    return B1peak, flipangle


    
def BS_B1_flipangle_wrt_ScanObject(POI_scan, shapefile_path, n_phasewraps,
                                   BS_scan, BS_procno, 
                                   BS_scan_neg_on=0, BS_scan_neg_on_procno=0,
                                   BS_scan_pos_off=0, BS_scan_pos_off_procno=0,
                                   BS_scan_neg_off=0, BS_scan_neg_off_procno=0):
    '''
    Return peak B1 & flipangle of pulse of interest, given some sarpy.io.BRUKER_classes Scan objects

    Arguments:
        Scan POI_scan: Scan object containing RF pulse of interest (POI)
        str shapefile_path: path of RF shape files
        Scan BS_scan: Scan object containing Bloch-Siegert data
            -if only 1 scan specified, this includes all data
            -if 4 scans specified, this is the powered +ve frequency scan
        int BS_procno: procno number with Bloch-Siegert phase data
        Scan BS_scan_neg_on: Scan object with powered -ve frequency scan
        int BS_scan_neg_on_procno: procno of powered -ve frequency scan
        Scan BS_scan_pos_off: Scan object with zero-power +ve frequency scan
        int BS_scan_pos_off_procno: procno of zero-power +ve frequency scan
        Scan BS_scan_neg_off: Scan object with zero-power -ve frequency scan
        int BS_scan_neg_off_procno: procno with zero-power -ve frequency scan
        int n_phasewraps:  number of phase wraps across whole volume
        
    Returns:
        image sets of peak B1 and flip angle for pulse of interest
    '''
    from sarpy.io.BRUKERIO import readRFshape
    import numpy as np    
    
    BS_method = BS_scan.method.Method
    assert BS_method == 'bSB1mapFLASH' or BS_method == 'bSB1mapMSME', \
                        'wrong scan type'    
     
    try:
        freqoffset = BS_scan.method.BSFreqOffset
        dBpwr_BS = BS_scan.method.BSPulse[3]
        width_BS = BS_scan.method.BSPulse[0]*1e-3
        #regex to strip off <>      
        BS_shapename = \
            re.findall(r'<\s*(\S+)\s*>',BS_scan.method.BSPulse[8])[0]  
        BS_pulseshape = \
            readRFshape(os.path.join(shapefile_path,BS_shapename))['amp']
        BS_pulseshape = numpy.asarray(BS_pulseshape)

        POI_nom_flipangle = POI_scan.method.ExcPulse[2]
        POI_dBpwr = POI_scan.method.ExcPulse[3]
        POI_integ_ratio = POI_scan.method.ExcPulse[10]
        POI_pulse_width = POI_scan.method.ExcPulse[0]*1e-3
        POI_shapename = \
                re.findall(r'<\s*(\S+)\s*>',POI_scan.method.ExcPulse[8])[0]
        POI_pulseshape = \
            readRFshape(os.path.join(shapefile_path,POI_shapename))['amp']
        POI_pulseshape = numpy.asarray(POI_pulseshape)
        acq_ndim = POI_scan.acqp.ACQ_dim
        if acq_ndim == 3:
            POI_nPE2 = POI_scan.acqp.ACQ_size[2]
        else:
            POI_nPE2 = 1
        
    except:
        print('something wrong happened with pulse parameter query')
        
    # load array data
    if BS_scan_neg_on == 0 and BS_scan_pos_off == 0 and BS_scan_neg_off == 0:
        # one expno contains all the BS data
        assert BS_scan.pdata[BS_procno-1].reco.RECO_image_type == \
            'PHASE_IMAGE', 'BS data is not a phase image'
        BS_data=BS_scan.pdata[BS_procno-1].data
        BS_dim=list(BS_data.shape)
        if len(BS_dim)==3:
            BS_data = numpy.reshape(BS_data,[BS_dim[0],BS_dim[1],1,BS_dim[2]])

        on_BSplus_phase = BS_data[:,:,:,0]
        on_BSminus_phase = BS_data[:,:,:,1]
        off_BSplus_phase = BS_data[:,:,:,2]
        off_BSminus_phase = BS_data[:,:,:,3]
        
    else:
        #four expnos contain the BS data
        BS_scan_pos_on_pdata = BS_scan.pdata[BS_procno-1]
        BS_scan_neg_on_pdata = BS_scan_neg_on.pdata[BS_scan_neg_on_procno-1]
        BS_scan_pos_off_pdata = BS_scan_pos_off.pdata[BS_scan_pos_off_procno-1]
        BS_scan_neg_off_pdata = BS_scan_neg_off.pdata[BS_scan_neg_off_procno-1]

        assert \
            BS_scan_pos_on_pdata.reco.RECO_image_type == \
            BS_scan_neg_on_pdata.reco.RECO_image_type == \
            BS_scan_pos_off_pdata.reco.RECO_image_type == \
            BS_scan_neg_off_pdata.reco.RECO_image_type == \
            'PHASE_IMAGE', \
            'some or all of BS data not a phase image'
        
        assert \
            BS_scan_pos_on_pdata.data.shape == \
            BS_scan_neg_on_pdata.data.shape == \
            BS_scan_pos_off_pdata.data.shape == \
            BS_scan_neg_off_pdata.data.shape, \
            'unequal data dimensions in the 4 BS expnos'
            
        BS_dim=BS_scan_pos_on_pdata.data.shape

        on_BSplus_phase = BS_scan_pos_on_pdata.data
        on_BSminus_phase = BS_scan_neg_on_pdata.data
        off_BSplus_phase = BS_scan_pos_off_pdata.data
        off_BSminus_phase = BS_scan_neg_off_pdata.data
        
        BS_dim = list(on_BSplus_phase.shape)
        if len(BS_dim)==2:
            on_BSplus_phase = np.reshape(on_BSplus_phase,BS_dim.append(1))
            on_BSminus_phase = np.reshape(on_BSminus_phase,BS_dim.append(1))
            off_BSplus_phase = np.reshape(off_BSplus_phase,BS_dim.append(1))
            off_BSminus_phase = np.reshape(off_BSminus_phase,BS_dim.append(1))
        

    if BS_method=='bSB1mapFLASH':
        seqtype = 'gradientecho'
    elif BS_method=='bSB1mapMSME':
        seqtype = 'spinecho'
        
    KBS = calc_KBS(freqoffset, width_BS, BS_pulseshape, seqtype)

    B1peak_BS = calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, 
                               off_BSminus_phase, off_BSplus_phase,
                               n_phasewraps)

    B1peak_POI = calc_POI_B1peak(B1peak_BS, dBpwr_BS, POI_dBpwr)
    
    flipangle_POI = calc_POI_flipangle(B1peak_POI, POI_integ_ratio,
                                       POI_pulse_width, POI_pulseshape,
                                       acq_ndim, POI_nPE2, 
                                       on_resonance = True )
        
    return B1peak_POI, flipangle_POI



def calc_POI_B1peak(B1peak_BS, dBpwr_BS, dBpwr_POI):
    '''
    Return peak B1 of a pulse of interest in Tesla
    
    Arguments:
        array B1peak_BS: peak B1 of Block-Siegert pulse
        float dBpwr_BS:  power of Bloch-Siegert pulse in dB
        float dBpwr_POI:  power of pulse of interest in dB
    Returns:for j=[1 3 4 8 9 10]

        array B1peak_POI
    '''    
    B1peak_POI = B1peak_BS * (math.pow(10,(dBpwr_BS-dBpwr_POI)/20))
    return B1peak_POI
 
def calc_norm_flipangle_profile(pulseshape, integ_ratio, pulse_width, interp_factor=100):
    '''
    returns the flip angle profile across a selective slice, divided by peak B1
    
    Arguments:
        array pulseshape: amplitude of pulse shape (arb. units)
        float integratio: integral ratio of pulse relative to block pulse
        float pulsewidth: length of pulse in seconds
        int interp_factor: interpolation factor
    Returns:
        array trunc_interp_profile: interpolated flip angle profile/peak B1
            (only the points within the bandwidth are included)
        float BW:  bandwidth of pulse in radians
        float interp_dw:  stepsize of interpolated pulse in radians

    Given the B1(t) pulse shape, the flip angle profile can be approximated
    as the fourier transform of B1(t) multiplied by gamma.  The width of the
    slice is equivalent to the bandwidth of the pulse.
    
    '''
    import scipy as sp
    import numpy as np
 
    gamma = 267.513e6 # radians/(s*T)
    
    # Fourier transform of B1(t) shape = shape of flip angle profile
    pulseprofile = abs(sp.fftpack.fftshift(sp.fftpack.fft(pulseshape)))
    n_profile = pulseprofile.size
    dt = pulse_width/pulseprofile.shape[0]
    df = 1/(dt*n_profile)
    dw = 2*math.pi*df
    
    # interpolate pulse to improve FWHM calc & later integration across profile
    profile_eqn = \
         sp.interpolate.UnivariateSpline(
             np.linspace(0,(n_profile-1)*dw,n_profile),pulseprofile)
    interp_profile = \
        profile_eqn(np.linspace(0,(n_profile-1)*dw,interp_factor*n_profile))

    # normalize pulse amplitude:  centre of slice (zero freq offset) has
    # flip angle = gamma*(time integral of normalized B1(t))*B1peak
    interp_profile = \
        interp_profile/max(interp_profile)*gamma*integ_ratio*pulse_width
        
    interp_dw = dw/interp_factor

    # truncate pulse to its FWHM        
    lo, hi = find_FWHM_limits(interp_profile)
    BW = dw/interp_factor*(hi-lo)
    trunc_interp_profile = interp_profile[lo:hi]
    
    if debug:
        pylab.plot(trunc_interp_profile)
        pylab.show()
    
    return trunc_interp_profile, BW, interp_dw
  
def calc_POI_flipangle(B1peak_POI, integ_ratio, pulse_width, pulseshape,
                        acq_ndim=2, nPE2=1, on_resonance=False):
    '''
    return flip angle of pulse of interest in degrees
    
    Arguments:
        array B1peak_POI: peak B1 of pulse of interest (POI)
        dict POI_parx: contains relevant parameters of POI
        array pulseshape: shape of pulse (maybe not normalized)
        bool on_resonance:  are spins assumed to be on resonance?
        int acq_ndim: number of dimensions in the acquisition (2D or 3D)
    Returns
        flip angle map of pulse of interest
    '''

    import scipy as sp
   
    gamma = 267.513e6 # radians/(s*T)

    if on_resonance:
        '''
        flip angle calc for on_resonant spins: gamma*(time integral of B1)
        '''
#        flipangle_POI = gamma*B1peak_POI*integ_ratio*pulse_width
        flipangle_POI = gamma*B1peak_POI*1*0.2e-3
    else:
        '''
        flip angle calc which accounts for off-resonant spins due to slice
        select gradient (assumes that slice profile = FT of gamma*|B1(t)|,
        which is valid at small flip angles:  
        '''
        sliceprofile, BW, dw = \
            calc_norm_flipangle_profile(pulseshape,integ_ratio,pulse_width)

        n_sliceprofile = sliceprofile.shape[0]
  
        if acq_ndim == 2:
            '''
            effective flip angle for 2D excitation = 
                (integral of flip angle across slice)/BW
            '''
            integrand = sliceprofile
            flipangle_integral = \
                sp.integrate.trapz(integrand, None, dw)
            flipangle_POI = flipangle_integral*B1peak_POI/BW
        elif acq_ndim == 3:
            '''
            effective flip angle for a PE2 slice in 3D excitation = 
                (integral of flip angle across PE2 section)/(BW/nPE2)
            '''
            n_interval = math.floor(n_sliceprofile/nPE2)
            interval = BW/nPE2
            flipangle_POI = numpy.zeros(B1peak_POI.shape)            
            
            for k in xrange(B1peak_POI.shape[2]):
                integrand = \
                    sliceprofile[k*n_interval:(k+1)*n_interval]
                flipangle_integral = \
                    sp.integrate.trapz(integrand, None, dw)
                flipangle_POI[:,:,k] = \
                    flipangle_integral*B1peak_POI[:,:,k]/interval
            
    return flipangle_POI*180/math.pi

def find_FWHM_limits(pulse):
    '''
    return the indices where the FWHM of a pulse is found
    
    :param array pulse: 1D array of pulse shape
    :returns:
        lo and hi indices of FWHM points
        
    example:
        >>> x,dx=numpy.linspace(-99,99,1000,retstep=True)
        >>> BW=50; c=BW/2.35482
        >>> gauss= 1*numpy.exp((-(x-0)**2/(2*c**2)))
        >>> lo,hi=find_FWHM_limits(gauss)
        >>> calc_BW = (hi-lo)*dx
        >>> numpy.testing.assert_approx_equal(BW,calc_BW,significant=2)
        
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


   
def calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, off_BSminus_phase,
                   off_BSplus_phase, n_phasewraps):
    '''
    calculates B1peak of the Bloch-Siegert pulse, based on phase differences  

    :param float KBS: Bloch-Siegert calibration constant in radians/T^2
    :param array on_BSminus_phase: imageset with BS pulse on (-ve offset)
    :param array on_BSplus_phase: imageset with BS pulse on (+ve offset)
    :param array off_BSminus_phase: imageset with zero flip angle BS pulse (-ve offset)
    :param array off_BSplus_phase: imageset with zero flip angle BS pulse (+ve offset)
    :param int n_phasewraps:  number of phase wraps across whole volume
    :
    :return:
        imageset of B1peak values of the Bloch-Siegert pulse
    :rtype: array
    '''

    pi = math.pi
    offset = off_BSplus_phase - off_BSminus_phase
    offset[offset>pi] = offset[offset>pi] - 2*pi
    offset[offset<-pi] = offset[offset<-pi] + 2*pi

 
    phase_diff = on_BSplus_phase - on_BSminus_phase - offset
    phase_diff[phase_diff>pi] = phase_diff[phase_diff>pi] - 2*pi
    phase_diff[phase_diff<-pi] = phase_diff[phase_diff<-pi] + 2*pi
#    phase_diff[phase_diff>0.75*pi] = phase_diff[phase_diff>0.75*pi] - 2*pi

    phase_diff = phase_diff + n_phasewraps*2*pi

    if debug:
        i_slice = math.floor(phase_diff.shape[2]/2)
#        pylab.imshow(180/pi*offset[:,:,i_slice],interpolation='None')
        pylab.imshow(180/pi*phase_diff[:,:,i_slice],interpolation='None')
        pylab.colorbar()
        pylab.show()

    B1peak_BS = numpy.sqrt(numpy.absolute(phase_diff)/(2*KBS))

    return B1peak_BS
    
   

def calc_KBS(freqoffset, pulsewidth, pulseshape, mode='gradientecho'):
    '''
    calculates KBS in the units of radians/T^2 (assumes freq offset is constant in time).  

    :param float freqoffset: BS offset in Hz
    :param float pulsewidth: length of pulse in seconds
    :param array pulseshape: amplitude of pulse shape
    :param string shapefile_path: directory of pulse_shape
    :param string mode: 'gradientecho' or 'spinecho'
    :return:
        value of KBS in radians/T^2
    :rtype: float
    
    Divide the result by 1e8 if units of radians/Gauss^2.  See literature reference in the module docstring for details.

Example:
    
    >>> shapefile_path = os.path.expanduser(os.path.join('~','wave'))
    >>> pulseshape = sarpy.io.BRUKERIO.readRFshape(os.path.join(shapefile_path,'fermi.exc'))['amp']
    >>> pulseshape = numpy.asarray(pulseshape)    
    >>> KBS = calc_KBS(4000,8e-3,pulseshape)
    >>> KBS
    7109603109.3280497
    '''
    amp = pulseshape/max(pulseshape)
    dt = pulsewidth/pulseshape.shape[0]

    gamma = 267.513e6 # radians/(s*T)
    integrand = (gamma*amp)**2/(2*2*3.141592654*freqoffset)
    KBS = scipy.integrate.trapz(integrand, None, dt)

    assert mode == 'gradientecho' or mode == 'spinecho', (
        'mode is not gradientecho or spinecho')
    
    if mode == 'spinecho':
        return 2*KBS
    else:
        return KBS


    
def BS_test_originaldata(BS_on_pos, BS_on_pos_procno,
                         BS_on_neg, BS_on_neg_procno,
                         BS_off_pos, BS_off_pos_procno,
                         BS_off_neg, BS_off_neg_procno,
                         POI_scan):
    '''
    Test function using Bloch-Siegert data contained in 4 scans (assumes phase
    is reconstructed in the same procno across all 4 scans)
    
    Shows the flip angle map of the pulse of interest    
    '''
    shapefile_path = os.path.expanduser(os.path.join('~','wave'))

    B1peak, flipangle = \
    BS_B1_flipangle_wrt_ScanObject(POI_scan, shapefile_path,
                                   BS_on_pos, BS_on_pos_procno,
                                   BS_on_neg, BS_on_neg_procno,
                                   BS_off_pos, BS_off_pos_procno,
                                   BS_off_neg, BS_off_neg_procno)
    
    
    i_slice = int(BS_on_pos.acqp.NSLICES/2)
    pylab.imshow(flipangle[:,:,i_slice],interpolation='none')
    pylab.title('flip angle map of pulse inside ' + POI_scan.shortdirname)
    pylab.colorbar()
    pylab.show()            
  
if __name__ == '__main__':
#    studyname = 'dBlochSiegert2.j41'
#    BS_expno = 1
#    BS_procno = 2
#    POI_expno = 2
#    uname = 'mriuser'
#    BS_test(studyname,BS_expno,BS_procno,POI_expno,uname)

    studyname = 'dCPMG_edgetest.kL1'
    BS_expno = 10
    BS_procno = 2
    POI_expno = 9
    uname = 'mriuser'
    n_phasewraps = 3
    BS_test(studyname,BS_expno,BS_procno,POI_expno,uname,n_phasewraps)
 
#    studypath = 'dBlochSiegert1.gP2'
#    BS_off_neg = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'36'))
#    BS_off_pos = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'37'))
#    BS_on_neg = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'46'))
#    BS_on_pos = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'47'))
#    POI_scan = sarpy.io.BRUKER_classes.Scan(os.path.join(studypath,'32'))
#    BS_test_originaldata(BS_on_pos,2,BS_on_neg,2,BS_off_pos,2,BS_off_neg,2,POI_scan)
    
#    import doctest
#    doctest.testmod()






