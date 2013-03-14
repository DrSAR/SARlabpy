# -*- coding: utf-8 -*-
from enthought.traits.api import (HasTraits, Instance, Array, Range, Float,
                                  Enum, on_trait_change, Property)
from enthought.traits.ui.api import View, Item, Group
from enthought.chaco.api import Plot, ArrayPlotData, OverlayPlotContainer
from enable.component_editor import ComponentEditor
from numpy import arange, max
import TwoCXM

class vis2CXM(HasTraits):
    plot = Instance(OverlayPlotContainer)
    time = Array
    dt=1./60                  # sample at 1s time resolution
    
    # declaration of model parameters, could this be encapsulated in 
    # some classes that are imported and registered with this tool?
    # some of the model traits are dependent on some of theit parameters and
    # there might be circular dependecies in cases of multiple representations
    # for one model
    
    AIFlabel = Enum("Lyng", 
                    "Parker", 
                    "HPG")
    aif = Property(Array, depends_on=['AIFlabel'])

    # std Tofts
    concentration1 = Property(Array, depends_on=['ve', 'Ktrans','AIFlabel'])
    ve = Range(low=0.0,high=2.0,value=.05)
    Ktrans = Range(low=.0,high=2.0,value=0.17)
    Grp_Tofts = Group(Item(name='ve'),
                      Item(name='Ktrans'),
                      Item(name='AIFlabel', label='Arterial Input Function'),
                      label="Tofts", dock='tab')

    # addition for XTofts
    concentration1X = Property(Array, depends_on=['ve', 'vp','Ktrans',
                                                  'AIFlabel'])
    vp = Range(low=0.0,high=1.0,value=.05)

    Grp_XTofts = Group(Item(name='ve'),
                       Item(name='vp'),
                       Item(name='Ktrans'),
                       Item(name='AIFlabel', 
                            label='Arterial Input Function'),
                       label="XTofts", dock='tab')    
    

    # 2CXM (see Sourbron)
    concentration2 = Property(Array, depends_on=['Kp', 'Km', 'Fp','Fm',
                                                 'AIFlabel'])
    Kp = Range(low=.0,high=10.0,value=0.17)
    Km = Range(low=.0,high=10.0,value=0.17)
    Fp = Range(low=.0,high=10.0,value=0.17)
    Fm = Range(low=.0,high=10.0,value=0.17)
    
    Grp_2CXM = Group(Item(name='Kp'),
                     Item(name='Km'),
                     Item(name='Fp'),
                     Item(name='Fm'),
                     Item(name='AIFlabel', 
                          label='Arterial Input Function'),
                     label="2CXM", dock='tab')
    
    # some plot characteristics, we should also have caller and line type...
    plot_type = Enum("line", "scatter")

                        
    # this is an internal view but maybe some of the alternatives should be 
    # tried. Question, can you add Items to an existing View
    traits_view = View( Item('plot',editor = ComponentEditor(),
                             resizable=True,
                             show_label=False),
                        Group(Grp_Tofts,
                              Grp_XTofts,
                              Grp_2CXM, layout="tabbed"),
                        #Item(name='plot_type'),
                        resizable = True,
                        buttons = ["OK"],
                        title='concentration time curves',
                        width=900, height=800)
                        
    #traits_view.add_trait(Item(name='plot_type'))

    def __init__(self):
        super(vis2CXM, self).__init__()
        self.plotdata = ArrayPlotData (x=self.time, 
                                       y1=self.concentration1,
                                       y1X=self.concentration1X,
                                       y2=self.concentration2, aif=self.aif)
        # tissue curves and aif have to be on different panels but overlaid 
        # since their max is quite different                                   
        tissuecurve = Plot(self.plotdata)
        tissuecurve.plot(("x", "y1"), type="line", color="blue")
        tissuecurve.plot(("x", "y1X"), type="line", color="black")
        tissuecurve.plot(("x", "y2"), type="line", color="green")
        aifcurve = Plot(self.plotdata)
        aifcurve.plot(("x", "aif"), type="line", color="red")
        
        container = OverlayPlotContainer(tissuecurve, aifcurve)
        tissuecurve.y_axis.title = "tissue concentration"
        tissuecurve.legend.visible = True 
        print tissuecurve.legend.plots.keys
        print tissuecurve.legend.plots.values
#        tissuecurve.legend.labels = ["Tofts", "Extended Tofts", "2CXM"]
#        
#        
#        legend_plot_dict = {}
        for i in tissuecurve.legend.plots.keys():
            for j in tissuecurve.legend.labels:
                print i,j
#                label_set = set(j.split(' '))
#                plot_set = set(i.split('|'))
#                if label_set.issubset(plot_set):
#                    legend_plot_dict[j] = self.legend.plots[i]
#        tissuecurve.legend.plots = legend_plot_dict
                    
                    
        aifcurve.y_axis.orientation = "right"
        aifcurve.y_axis.title = "AIF - concentration"

        self.plot = container
        

    def _time_default(self):
        """ Default handler for time Trait Array. """
        t=arange(0,6, self.dt)   # time running from 0 to 6min
        return t
    
    def _get_aif(self):
       return (TwoCXM.aif(self.time,aifchoice=self.AIFlabel))

    def _get_concentration1(self):
        """Recalculate when any trait the property depends on changes."""
        modelparams = (self.Ktrans, self.ve)
        ca = TwoCXM.aif(self.time,aifchoice=self.AIFlabel, zeropad=True)
        return TwoCXM.conc_Tofts(self.time, modelparams, ca, self.dt)

    def _get_concentration1X(self):
        """Recalculate when any trait the property depends on changes."""
        modelparams = (self.vp, self.Ktrans, self.ve)
        ca = TwoCXM.aif(self.time,aifchoice=self.AIFlabel, zeropad=True)
        return TwoCXM.conc_XTofts(self.time, modelparams, ca, self.dt)

    def _get_concentration2(self):
        """Recalculate when any trait the property depends on changes."""
        modelparams = (self.Kp, self.Km, self.Fp, self.Fm)
        ca = TwoCXM.aif(self.time,aifchoice=self.AIFlabel, zeropad=True)
        return TwoCXM.conc_2CXM(self.time, modelparams, ca, self.dt)

    def _AIFlabel_changed(self):
        self.plotdata.set_data("aif", self.aif)
        
    def _concentration1_changed(self):
        self.plotdata.set_data("y1", self.concentration1)
        
    def _concentration1X_changed(self):
        self.plotdata.set_data("y1X", self.concentration1X)
        
    def _concentration2_changed(self):
        self.plotdata.set_data("y2", self.concentration2)

if __name__ == '__main__':
    viewer = vis2CXM()
    viewer.configure_traits()