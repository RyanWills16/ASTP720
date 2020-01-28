#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:13:24 2020

@author: ryan
"""

import rootfind as rt
import interpolation as ip
import numpy as np
import matplotlib.pyplot as plt


# problem #1 

def func(x):
    return(x**2 - 45)
    
def dfunc(x):
    return(2*x)


rt.bisection(func,4,11, 0.0000000001, verbose = True, numiter = True)

rt.Newton(func, dfunc, 1, 0.00000000000001, verbose = True, numiter = True)

rt.secant(func, -2, -1, verbose = True, numiter = True)


# Problem #2

x = np.arange(-4.5,5, 0.5)

y =  2*x + 11 -x**2

x2 = []
y2 = []

for i in np.arange(-4.5,x[-1],0.25):
    newpoint = ip.linterp(x,y)
    y_new = newpoint(i)
    
    x2.append(i)
    y2.append(y_new)
    
plt.scatter(x2,y2, marker = 'o', edgecolors = 'r', facecolors = 'none')

plt.scatter(x,y, marker = 'o', edgecolors = 'b', facecolors = 'none')




