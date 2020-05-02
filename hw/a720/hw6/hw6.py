# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 00:48:58 2020

@author: Ryan
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks as fp


lc = np.genfromtxt("lightcurve_data.txt")
rv = np.genfromtxt("RV_data.txt")

jd = lc[:,0] - lc[0,0]
peaks = (lc[:,1]-1)*-1

# find peaks and subtract their x value to create folded rotation curve.
t = fp(peaks,height = 0.007, width = 5)
light = []
phase = []
for ind, i in enumerate(lc[:,1]):
    for ind2, j in enumerate(t[0]):
        if j == ind:
            for k in lc[ind-10:ind+10,1]:
                light.append(k)
                
            for l in jd[ind-10:ind+10] - jd[ind]:
                phase.append(l)
            
        else: 
            continue
        
# plot the resulting light curve
phase = np.array(phase)
light = np.array(light)

plt.figure()
plt.scatter(phase, light)
plt.xlabel('Phase')
plt.ylabel('Amplitude')


def box(x, center, width, amplitude, offset = 1):
    '''
    Box model
    Parameters
    ----------
    x : array or list
    c : float
        center.
    w : float
        half width.
    amp : float
        amplitude.
    '''
    
    a = []
    
    for i in x:
        
        if i <= (center-width) or i >= (center+width):
            a.append(float(offset))
        else:
            a.append(float(offset + amplitude))
    return a


# by eye fit of the light curve
ybox = box(phase[0:19], 0, 0.083,-0.00717)
plt.figure()
plt.scatter(phase, light)
plt.plot(phase[0:19], ybox, 'r')
plt.xlabel('Phase')
plt.ylabel('Amplitude')

def Like(y, signal, sigma = 1):
    '''
    calculate likelihood
    '''

    return (1/np.sqrt(2*np.pi*sigma**2))*(np.exp(-0.5*((y - signal)/sigma)**2))


def MCMC(it, x, y, init, goffset=1):
    '''
    x: array of x values
    y: array of y values
    
    given initial guesses for the parameters

    '''
    
    iteration = [0]
    
    c_fit = [init[0]]
    w_fit = [init[1]]
    amp_fit = [init[2]]
    
    count = 0
    
    while count <= it:
        iteration.append(count)
        count += 1
        
        init_model = box(x, init[0], init[1], init[2])
        Likelihood_prior = np.prod(Like(y, init_model))
        
        # list of propsal parameters
        proposal = [float(np.random.normal(init[0], abs(0.05*init[0]), 1)),float(np.random.normal(init[1], abs(0.05*init[1]), 1)), float(np.random.normal(init[2], abs(0.05*init[2]), 1))]
        
        # each box model uses one proposal parameter and current parameters
        p1 = box(x, proposal[0], init[1], init[2])
        p2 = box(x, init[0], proposal[1], init[2])
        p3 = box(x, init[0], init[1], proposal[2])
        
        # Likelihood
        L1 = np.prod(Like(y,p1))
        L2 = np.prod(Like(y,p2))
        L3 = np.prod(Like(y,p3))
        
        # metropolis ratios
        if L1/Likelihood_prior > np.random.rand():
            c_fit.append(proposal[0])
            init[0] = proposal[0]
        else:
            c_fit.append(init[0])
            init[0] = init[0]
            
        if L2/Likelihood_prior > np.random.rand():
            w_fit.append(proposal[1])
            init[1] = proposal[1]
        else:
            c_fit.append(init[1])
            init[1] = init[1]
            
        if L3/Likelihood_prior > np.random.rand():
            amp_fit.append(proposal[2])
            init[2] = proposal[2]
        else:
            c_fit.append(init[2])
            init[2] = init[2]
                    
    return iteration, c_fit, w_fit, amp_fit


init = [0, 0.08,-0.007]
t = MCMC(10000, phase, light, init)
  
n = np.linspace(0,10001,10002 )
plt.plot(n, t[3])
plt.ylabel('Amplitude')
plt.xlabel('Iteration')

center = np.mean(t[1])
width = np.mean(t[2])
amplitude = np.mean(t[3])


