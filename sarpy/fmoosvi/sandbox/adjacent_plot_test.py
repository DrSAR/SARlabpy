# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 12:14:45 2013

@author: firas
"""


import numpy
import pylab

a=numpy.empty([4,50])
a[0,:] = numpy.array(numpy.linspace(0,100,num=50))
a[1,:] = numpy.linspace(0,200,num=50)
a[2,:] = numpy.linspace(0,300,num=50)
a[3,:] = numpy.linspace(0,50,num=50)

fig = pylab.figure()
G = pylab.matplotlib.gridspec.GridSpec(1, 4, wspace=0.0, hspace=0.0)   

for col in xrange(4):

    fig.add_subplot(G[0, col])
    pylab.plot(a[col,:])
        
    xlim(0,80)
    ylim(0,350)
    pylab.locator_params(axis='both',which='both',nbins=3,labelsize=3)
    
    if col != 0:
        pylab.gca().axes.get_yaxis().set_visible([])
        pylab.gca().axes.get_xaxis().set_visible([])

    ax = pylab.gca()
    
    for label in ax.get_xticklabels():
        label.set_fontsize(8)
        pylab.text(pylab.xlim()[1]*0.5,pylab.ylim()[1]*0.5,'Test',horizontalalignment='center')
        
#    ax.tick_params(direction='in', pad=5)


#fig.tight_layout()
    

