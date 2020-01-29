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

x = np.arange(-4.5,5, 0.5) # randomly chosen range of numbers

y =  2*x + 11 -x**2 # randomly chosen quadratic function to get a set of y values


p = ip.linterp(x,y) # call the linear interpolation function to prime it with sets of points x and y

p(-4.5, verbose = True) # interpolate a random point to check that the function works


x2 = [] # defining a new empty set of x values
y2 = [] # defining a new empty set of y values

for i in np.arange(-4.5,x[-1],0.25):    # iterating over a range of new x points at a higher resolution
    newpoint = ip.linterp(x,y)  # priming the function with the original x, y data points
    y_new = newpoint(i)         # using linear interpolation to find new x data points
    
    x2.append(i)        # appending new points to empty array
    y2.append(y_new)    # appending new y points to empty array
    
plt.scatter(x2,y2, marker = 'o', edgecolors = 'r', facecolors = 'none', label = 'interpolated points') # plot new interpolated points

plt.scatter(x,y, marker = 'o', edgecolors = 'b', facecolors = 'none', label = 'Original points') # plot original points in a different color

plt.legend()

# Problem 3






