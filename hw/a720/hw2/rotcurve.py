# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 22:59:56 2020

@author: Ryan
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u

def velfunc(x, c, v_200):
    '''
    returns function for v_c^2
    '''
    return np.sqrt((v_200**2/x)*((np.log(1 + c*x) - (c*x)/(1 + c*x))/(np.log(1 + c) - (c)/(1 + c))))

def massfunc(r, velocity, G):
    return r*velocity/G


radius = np.linspace(0.0001,1,1000)
radius2 = np.linspace(0,300, 1000)


velocity = velfunc(radius, 15, 200)
mass = massfunc(radius2, velocity, 4.299e-6)
                
    