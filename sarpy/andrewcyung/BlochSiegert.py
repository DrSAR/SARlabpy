# -*- coding: utf-8 -*-
# Copyright (C) 2012-2013 Andrew Yung and UBC 7T MRI lab
"""
Analysis routines for Bloch-Siegert B1 and flip angle mapping on Bruker data. 

Ref:  B1 mapping by Bloch-Siegert shift. Sacolick LI, Wiesinger F, Hancu I, 
Vogel MW.  Magn Reson Med. 2010 May;63(5):1315-22. 

Brief description of method:  off-resonance pulse applied after spin excitation, far away from Larmor frequency that it does not perturb magnetization, but strong enough to perturb direction of magnetic field.  Resultant phase shift is related to the peak B1 of the off-resonance pulse.  This B1 measurement can be used to map B1 and flip angle of other pulses (with knowledge of relative transmit power levels).

The technique requires at least two scans, repeated at positive and negative offset centered around the Larmor frequency.  The phase difference is proportional to B1, which is independent of B0 inhomogeneity up to the first order.  As of version 1.0, the method requires an additonal two scans with zero-power off-resonant pulses to determine the phase shift originating from switching the frequency synthesizer for the Bloch-Siegert pulse.

The analysis routines assume that the user is interested in using the peak B1 from the Bloch-Siegert experiment to determine the B1 or flip angle of a RF pulse of interest ("POI" from a different scan).  The user must specify a Bruker pulse shape directory to allow calculation of KBS constant.

Several function wrappers exist for calculation of flip angle maps of a pulse of interest (POI) which differ by their input type:
    
1) BS_flipangle_wrt_Bruker2dseq: file path of Bloch-Siegert expno and POI expno 

2) BS_flipangle_wrt_ScanObject: sarlabpy Scan objects

3) calc_BS_flipangle:  arrays for the individual Bloch-Siegert acquired scans

write_BS_flipangle_2dseq creates a 2dseq and visu_pars file for the calculated flip angle map, thus enabling direct visualization of the results in Paravision image display. 

Test functions assume that sample data is available in the ~/data, and the RF shape files are in ~/wave.  Use symbolic links if they are somewhere else.
"""

from __future__ import division

import numpy
import os
import pylab
import sarpy
import sarpy.utils.plt_utils
import pdb
import scipy.integrate
import scipy.optimize
import scipy.fftpack
import math

def BS_test(studyname,BS_expno,BS_procno,POI_expno,offConsole=True):
    
    BS_path = os.path.join(studyname,str(BS_expno))
    POI_path = os.path.join(studyname,str(POI_expno))

    if offConsole:
        shapefile_path = os.path.expanduser(os.path.join('~','wave'))
    else:
        shapefile_path = 'opt/PV5.1/exp/stan/nmr/lists/wave'

    flipanglemap = BS_flipangle_wrt_BrukerPath(BS_path,BS_procno,POI_path,shapefile_path)
    
    i_slice = 4
    pylab.imshow(flipanglemap[:,:,i_slice],interpolation='none')
    pylab.show()


def BS_flipangle_wrt_BrukerPath(BS_path,BS_procno,POI_path,shapefile_path):

    BS_scan = sarpy.io.BRUKER_classes.Scan(BS_path)
    POI_scan = sarpy.io.BRUKER_classes.Scan(POI_path)
    
    flipanglemap = BS_flipangle_wrt_ScanObject(BS_scan, BS_procno, POI_scan, shapefile_path)

    return flipanglemap
    
def BS_flipangle_wrt_ScanObject(BS_scan, BS_procno, POI_scan, shapefile_path):

    import re
    
    try:
        BS_method = BS_scan.method.Method
        dBpwr_BS = BS_scan.method.BSPulse[3]
        dBpwr_POI = POI_scan.method.ExcPulse[3]
        integralratio_POI = POI_scan.method.ExcPulse[10]
        width_BS = BS_scan.method.BSPulse[0]*1e-3
        width_POI = POI_scan.method.ExcPulse[0]*1e-3
        freqoffset = BS_scan.method.BSFreqOffset
        BS_shape = BS_scan.method.BSPulse[8]
        BS_shape = re.findall(r'<\s*(\S+)\s*>',BS_shape)[0]  #regex to strip off <>      
    except:
        print('something wrong happened with pulse parameter query')
        
    try:
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

    B1peak_BS = calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, off_BSminus_phase, off_BSplus_phase)
    POI_maps = calc_POI_B1_flipangle(B1peak_BS, dBpwr_BS, dBpwr_POI, integralratio_POI, width_POI)
        
    return POI_maps['flipangle']
    
def calc_POI_B1_flipangle(B1peak_BS, dBpwr_BS, dBpwr_POI, integralratio_POI, width_POI):
    
    gamma = 267.513e6
    B1peak_POI = B1peak_BS * (math.pow(10,(dBpwr_BS-dBpwr_POI)/20))
    flipangle_POI = gamma*B1peak_POI*integralratio_POI*width_POI
         
    return {'flipangle':flipangle_POI, 'B1peak':B1peak_POI}
    
def calc_BS_B1peak(KBS, on_BSminus_phase, on_BSplus_phase, off_BSminus_phase, off_BSplus_phase):
    '''
    calculates B1peak of an off-resonant pulse, based on phase differences  

    :param float KBS: Bloch-Siegert calibration constant in radians/T^2
    :param array on_BSminus_phase: imageset with BS pulse on (-ve offset)
    :param array on_BSplus_phase: imageset with BS pulse on (+ve offset)
    :param array off_BSminus_phase: imageset with zero flip angle BS pulse (-ve offset)
    :param array off_BSplus_phase: imageset with zero flip angle BS pulse (+ve offset)
    :return:
        imageset of B1peak values of the Bloch-Siegert pulse
    :rtype: array
    
Example:
    
    >>> shapefile_path = os.path.expanduser(os.path.join('~','wave'))
    >>> KBS = calc_KBS(4000,8e-3,'fermi.exc',shapefile_path)
    >>> KBS
    7109603109.3280497
    '''
    debug=True
    
    offset = off_BSplus_phase - off_BSminus_phase
    
    if debug:
        pylab.imshow(offset[:,:,4],interpolation='None')
        pylab.show()
        
    phase_diff = on_BSplus_phase - on_BSminus_phase + offset
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
    
    Divide the result by 1e8 if units of radians/Gauss^2.  See literature reference in      the module docstring for details.

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
    
            
   
if __name__ == '__main__':
#    shapefile_path = 'opt/PV5.1/exp/stan/nmr/lists/wave'
#    shapefile_path = os.path.expanduser(os.path.join('~','wave'))
#    KBS = calc_KBS(4000,8e-3,'fermi.exc',shapefile_path)
#    print KBS/1e8
    BS_test('dBlochSiegert2.j41',3,2,1)
    import doctest
    doctest.testmod()






