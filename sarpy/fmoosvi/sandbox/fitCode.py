# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 14:19:16 2013

@author: firas
"""

# Fitting code

import pylab
import scipy
import scipy.optimize

num_points = 150
Tx = linspace(5., 8., num_points)
Ty = Tx

tX = 11.86*cos(2*pi/0.81*Tx-1.32) + 0.64*Tx+4*((0.5-rand(num_points))*exp(2*rand(num_points)**2))
tY = -32.14*cos(2*pi/0.8*Ty-1.94) + 0.15*Ty+7*((0.5-rand(num_points))*exp(2*rand(num_points)**2))

# Create a functiont to data to 
fitfunc = lambda p, x: p[0]*cos(2*pi/p[1]*x+p[2]) + p[3]*x # Target function
errfunc = lambda p, x, y: fitfunc(p, x) - y # Distance to the target function

# Initial guesses for parameters
p0 = [-15., 0.8, 0., -1.] # Initial guess for the parameters
p1, success = scipy.optimize.leastsq(errfunc, p0[:], args=(Tx, tX))

time = linspace(Tx.min(), Tx.max(), 100)

plot(Tx, tX, "ro", time, fitfunc(p1, time), "r-")