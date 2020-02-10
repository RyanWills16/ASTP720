# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 22:59:56 2020

@author: Ryan
"""

import matplotlib.pyplot as plt
import numpy as np
import numcalc as nm
from astropy import units as u

def velfunc(x, c, v_200):
    '''
    returns function for v_c^2
    '''
    return np.sqrt((v_200**2/x)*((np.log(1 + c*x) - (c*x)/(1 + c*x))/(np.log(1 + c) - (c)/(1 + c))))

def massfunc(r, velocity, G):
    return r*velocity**2/G

G = 4.299e-6*u.kpc*u.km**2/u.solMass
c = 15
v_200 = 200*u.km/u.second
r_c = 250*u.kpc


x = np.linspace(0.0001,2,1000)
radius = x*r_c


velocity = velfunc(x, c, v_200)
mass = massfunc(radius, velocity, G)

plt.figure()
plt.plot(radius, velocity)
plt.xlabel('r')
plt.ylabel('velocity (km/s)')
                
plt.figure()
plt.plot(radius, mass)
plt.xlabel('radius')
plt.ylabel('Mass enclosed (solar masses)')    

def rotcurve(r_c, v_200, c, velfunc, massfunc):
    for 
    