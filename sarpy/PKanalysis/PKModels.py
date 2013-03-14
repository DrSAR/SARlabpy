#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Collection of pharmaco-kinetic model classes. Currently provides:
   - 2CXM
   - Tofts
   - extended Tofts
   
Created on Tue Apr 24 22:36:10 2012
@author: *Stefan A Reinsberg*, UBC, Vancouver
"""

#traits import
from enthought.traits.api import (HasTraits, Array, Range, Float,
                                  Enum, on_trait_change, Property,
                                  Tuple)
from traitsui.api import View, Item, Group, TupleEditor, RangeEditor
                                 
from numpy import sqrt                                 
    
    
def paramconv_modeltosoln(params):
    ''' convert parameters used to set up the model (ve, vp, PS, Fpl) into 
    parameters more useful in the solution (Fp, Fm, Kp, Km)
    
    function paramconv_solntomodel() is the inversion'''
    
    (ve, vp, PS, Fpl) = params
    
    e = ve / (vp + ve)  # extravascular fraction of the extracellular volume
    
    E = PS / (PS + Fpl)  # extraction fraction
    
    root = sqrt(1-4*E*e*(1-E)*(1-e)/((E-E*e+e)**2))
    taup = (E - E*e +e)/(2*E) * (1 + root)
    taum = (E - E*e +e)/(2*E) * (1 - root)
    
    Kp = Fpl / ((vp + ve)*taum)
    Km = Fpl / ((vp + ve)*taup)
    
    
    Fp =  Fpl * (taup - 1)/(taup - taum)
    Fm = -Fpl * (taum - 1)/(taup - taum)
    
    return (Fp, Fm, Kp, Km)

def paramconv_solntomodel(params):
    ''' convert parameters more useful in the solution (Fp, Fm, Kp, Km) into
    parameters used to set up the model (ve, vp, PS, Fpl) 
    
    
    this is the inverse of paramconf_modeltosoln()'''
    
    (Fp, Fm, Kp, Km) = params
    
    taum = (Fp/Fm+1) / (Fp/Fm + Kp/Km)
    taup = taum * Kp/Km
    
    Fpl = Fp * (taup- taum) / (taup - 1)
    
    E = (taup + taum -1 - taup*taum) / ((taup+taum)*(taup+taum-1)-taup*taum)
    
    e = (taup + taum -1 - taup*taum)  / (taup+taum-1)    
    
    vpplve = Fpl/Kp/taum
    ve = e * vpplve
    
    vp =  vpplve - ve
    
    PS = Fpl * E / (1-E)

    return (ve, vp, PS, Fpl)
                             
class PharmocoKineticModel(HasTraits):
    time = Array
    AIFlabel = Enum("Lyng", 
                    "Parker", 
                    "HPG")
    #concentration = Property(Array)
    
    
    
    
#here we should use a decorator to register this class with the visualizer
class TwoCXModel(PharmocoKineticModel):
    
    # set of parameters that makes sense in the physiological model    
    param_model = Tuple((.05, .05, .29, .17),
                        editor=TupleEditor(labels=['Ve', 'Vp', 'PS', 'Flow']))
    # parameters that are useful in the solution and map to the model parameter                    
    param_solution = Property(Tuple(0.,0.,0.,0.),
                              editor = TupleEditor(labels=['Fp','Fm','Kp','Km']),
                              depends_on=['param_model'])

    def __init__(self):
        #this group could be used if enclosed in a different viewer
        super(TwoCXModel, self).__init__()
        self.ViewGroupinst = Group(Item(name="param_model", 
                                    label="Model parameter",
                                    style='custom'),
                               Item(name="param_solution", 
                                    label="Solution parameter"))

    ViewGroup = Group(Item(name="param_model", 
                            label="Model parameter",
                            style='custom'),
                       Item(name="param_solution", 
                            label="Solution parameter"))
    #would be nice to mix sliders as in Range and Tuple Editor
    traits_view = View(ViewGroup,
                       buttons=['OK', 'Cancel'], resizable=True)

    def _get_param_solution(self):
        param_solution =  paramconv_modeltosoln(self.param_model)
        return param_solution
                       
    def _set_param_solution(self, value):
        self.param_model=paramconv_solntomodel(value)
        

class ToftsModel(PharmocoKineticModel):
    param_model = Tuple((0.05, 0.17),
                        editor=TupleEditor(labels=['ve','Ktrans']))
                        
    ViewGroup = Group(Item(name="param_model", label="Tofts Model parameter",
                           style='custom'))
    
    traits_view = View(ViewGroup,
                       buttons=['OK', 'Cancel'], resizable=True)
    
if __name__ == '__main__':
    #should run some standard or unit tests
    
    model1 = TwoCXModel()
    model1.configure_traits()
    
    model2 = ToftsModel()                          
    model2.configure_traits()
    
    view = View(model1.trait("ViewGroup"), model2.trait("ViewGroup"),
                buttons=['OK', 'Cancel'], resizable=True)
                