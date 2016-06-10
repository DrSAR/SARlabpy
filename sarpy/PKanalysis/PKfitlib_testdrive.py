# -*- coding: utf-8 -*-
"""
Attempts to test-drive the PK fit library while investigating the findings
by Sourbron and friends.
"""

import numpy
import sarpy.PKanalysis.PKfitlib as PKfitlib
import matplotlib.pyplot as plt

def parameter_stability(N=100):

    #create some data
    t=numpy.arange(400.)/60
    
    # set up AIFs
    aif_Gd = PKfitlib.aif(t,aifchoice='Lyng')
    ca_Gd = numpy.concatenate((numpy.zeros(len(t)-1), aif_Gd))

    aif_HPG = PKfitlib.aif(t,aifchoice='HPG')
    ca_HPG = numpy.concatenate((numpy.zeros(len(t)-1), aif_HPG)) 

    # create and fit intermediate case
    print('Gd fit with intermediate vascularization')
    p0 = [25./100, 50./100, 5./100, 5./100]
    print(p0)
    param1_i=[]
    param2_i=[]
    for i in numpy.arange(N):
        c_ivasc = PKfitlib.conc_2CXM(t, p0, ca_Gd) + numpy.random.randn(len(t))/4
        p1, success = PKfitlib.fit_2CXM(t, ca_Gd, c_ivasc, p0)
        param1_i.append(p1[1])
        param2_i.append(p1[2])

    
    # create and fit high vasc case
    print('Gd fit with high vascularization')
    p0 = [25./100, 50./100, 0.01, 0.09]
    print(p0)
    param1_h=[]
    param2_h=[]
    for i in numpy.arange(N):
        c_hvasc = PKfitlib.conc_2CXM(t, p0, ca_Gd) + numpy.random.randn(len(t))/4
        p1, success = PKfitlib.fit_2CXM(t, ca_Gd, c_hvasc, p0)
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
        c_hvasc = PKfitlib.conc_2CXM(t, p0[0:4], ca_Gd) + numpy.random.randn(len(t))/4
        cHPG_hvasc = PKfitlib.conc_2CXM(t, numpy.r_[p0[4],p0[1:4]], ca_HPG
                          ) + numpy.random.randn(len(t))/4
        p1, success = PKfitlib.fit_double_2CXM(t, ca_Gd, ca_HPG, 
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
    p0 = [25./100, 50./100, .05, 0.05]
    
    # set up AIFs
    aif_Gd = PKfitlib.aif(t,aifchoice='Lyng')
    ca_Gd = numpy.concatenate((numpy.zeros(len(t)-1), aif_Gd))

    aif_HPG = PKfitlib.aif(t,aifchoice='HPG')
    ca_HPG = numpy.concatenate((numpy.zeros(len(t)-1), aif_HPG)) 

    # create and fit intermediate case
    c_ivasc = PKfitlib.conc_2CXM(t, p0, ca_Gd) + numpy.random.randn(len(t))/4
    
    print('Gd fit with intermediate vascularization')
    print(p0)
    p1, success = PKfitlib.fit_2CXM(t, ca_Gd, c_ivasc, p0)
    print(p1)
    print((p1-p0)/p0*100)
    fit_ivasc = PKfitlib.conc_2CXM(t, p1, ca_Gd)    

    # create and fit highly vascular case    
    p0 = [1./100, 9./100, 25./100, 50./100, 0.01/100]
    p0 = [25./100, 50./100, .01, .09, 0.01/100]
    c_hvasc = PKfitlib.conc_2CXM(t, p0[0:4], ca_Gd) + numpy.random.randn(len(t))/4
    print('Gd fit with high vascularization')
    print(p0)
    p1, success = PKfitlib.fit_2CXM(t, ca_Gd, c_hvasc, p0[0:4])
    print(p1)
    print((p1-p0[0:4])/p0[0:4]*100)
    fit_hvasc = PKfitlib.conc_2CXM(t, p1[0:4], ca_Gd)
    
    # create HPG case (PS=0)
    cHPG_hvasc = PKfitlib.conc_2CXM(t, numpy.r_[p0[4], p0[1:4]], ca_HPG
                     ) + numpy.random.randn(len(t))/4
    print('Gd fit with high vascularization and PS=0')
    print(p0)
    p1, success = PKfitlib.fit_2CXM(t, ca_HPG, cHPG_hvasc, numpy.r_[p0[4],p0[1:4]])
    print(p1)
    print((p1-p0[0:4])/p0[0:4]*100)
    fitHPG_hvasc = PKfitlib.conc_2CXM(t, p1, ca_HPG)
    
    #double fit
    print('double fit')
    print(p0)
    p1, success = PKfitlib.fit_double_2CXM(t, ca_Gd, ca_HPG, c_hvasc, cHPG_hvasc, p0)
    print(p1)
    print((p1-p0)/p0*100)
    fitd_hvasc = PKfitlib.conc_2CXM(t, p1[0:4], ca_Gd)
    fitHPGd_hvasc = PKfitlib.conc_2CXM(t, numpy.r_[p0[0:2],p0[4],p0[3]], ca_HPG)

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
    aif_Gd = PKfitlib.aif(t,aifchoice='Lyng')
    ca_Gd = PKfitlib.aif(t,aifchoice='Lyng', zeropad=True)
    
    aif_HPG = PKfitlib.aif(t,aifchoice='HPG')
    ca_HPG = PKfitlib.aif(t,aifchoice='HPG', zeropad=True)
    
    # intermediately vascualarized
    ve = 20./100
    vp = 10./100
    PS = 25./100
    Fpl = 50./100
    print(PKfitlib.paramconv_modeltosoln((ve, vp, PS, Fpl)))
    c_ivasc = PKfitlib.conc_2CXM(t, (ve, vp, PS, Fpl), ca_Gd)
        
    # intermediately vascularized and slow tracer exchange regime
    ve =  20./100
    vp =  10./100
    PS =   0.01/100
    Fpl = 50./100
    cHPG_ivasc = PKfitlib.conc_2CXM(t, (ve, vp, PS, Fpl), ca_HPG)
        
    # weakly vascularized
    ve = 29./100
    vp = 1./100
    PS = 25./100
    Fpl = 50./100
    c_wvasc = PKfitlib.conc_2CXM(t, (ve, vp, PS, Fpl), ca_Gd)
        
    # weakly vascularized and slow tracer exchange regime
    ve =  29./100
    vp =   1./100
    PS =   0.01/100
    Fpl = 50./100
    cHPG_wvasc = PKfitlib.conc_2CXM(t, (ve, vp, PS, Fpl), ca_HPG)

    # highly vascularized
    ve = 1./100
    vp = 29./100
    PS = 25./100
    Fpl = 50./100
    c_hvasc = PKfitlib.conc_2CXM(t, (ve, vp, PS, Fpl), ca_Gd)

    # highly vascularized and low tracer exchange
    ve = 1./100
    vp = 29./100
    PS = 0.01/100
    Fpl = 50./100
    cHPG_hvasc = PKfitlib.conc_2CXM(t, (ve, vp, PS, Fpl), ca_HPG)
        
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
    ca = PKfitlib.aif(t, aifchoice='Parker',zeropad=True)
     # intermediately vascularized
    vp = 10./100 # ml/ml
    ve = 20./100 # ml/ml
    PS = 25./100 # ml/min/ml
    Fpl = 50./100 # ml/min/ml
    ''' consider equation [32] and [33] to obtain a Ktrans,
    according to [32] it is Ktrans = 17  # ml/min/100ml '''
    Ktrans = PS*Fpl / (PS+Fpl) # ml/min/ml - eqn [32]
  

    c_2CXM = PKfitlib.conc_2CXM( t, (PS, Fpl, ve, vp), ca, dt)
    c_XTofts =  PKfitlib.conc_XTofts(t,(Ktrans, ve, vp), ca, dt)
    c_Tofts =  PKfitlib.conc_Tofts(t,(Ktrans, ve), ca, dt)
    plt.subplot(331)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,c_XTofts,'b:')
    plt.plot(t,c_Tofts,'g--')
    plt.axis([0,6,0,.4])

    '''now fit the data fit_generic(t, ca, data, p0, model):'''
    resTofts  = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_Tofts, (.17, 0.2), ca, dt)
    resXTofts = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_XTofts, (.17, 0.2, .1), ca, dt)
    res2CXM   = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_2CXM, (.25, .5,
                                                                     .2, .1), ca, dt)

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
     
    c_2CXM = PKfitlib.conc_2CXM( t, (PS, Fpl, ve, vp), ca, dt)
    c_XTofts =  PKfitlib.conc_XTofts(t,(Ktrans, ve, vp), ca, dt)
    c_Tofts =  PKfitlib.conc_Tofts(t,(Ktrans, ve), ca, dt)
    plt.subplot(334)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,c_XTofts,'b:')
    plt.plot(t,c_Tofts,'g--')
    plt.axis([0,6,0,.25])

    '''now fit the data fit_generic(t, ca, data, p0, model):'''
    resTofts  = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_Tofts, (.17, 0.2), ca, dt)
    resXTofts = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_XTofts, (.17, 0.2,
                                                                      .1), ca, dt)
    res2CXM   = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_2CXM, (.25, .5,
                                                                     .2, .1), ca, dt)

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

    c_2CXM = PKfitlib.conc_2CXM(t,(PS, Fpl, ve, vp), ca, dt)
    c_XTofts =  PKfitlib.conc_XTofts(t,(Ktrans, ve, vp), ca, dt)
    c_Tofts =  PKfitlib.conc_Tofts(t,(Ktrans, ve), ca, dt)
    plt.subplot(337)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,c_XTofts,'b:')
    plt.plot(t,c_Tofts,'g--')
    plt.axis([0,6,0,.4])

    '''now fit the data fit_generic(t, ca, data, p0, model):'''
    resTofts  = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_Tofts, (.17, 0.2), ca, dt)
    resXTofts = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_XTofts, (.17, 0.2,
                                                                      .1), ca, dt)
    res2CXM   = PKfitlib.fit_generic(t, c_2CXM, PKfitlib.conc_2CXM, (.25, .5,
                                                                     .2, .1), ca, dt)

    print resTofts[0]
    print resXTofts[0]
    print res2CXM[0]
    
    plt.subplot(338)
    plt.plot(t,c_2CXM,'r-')
    plt.plot(t,resTofts[2],'g--')
    plt.plot(t,resXTofts[2],'b:')


    return (resTofts, resXTofts, res2CXM)
