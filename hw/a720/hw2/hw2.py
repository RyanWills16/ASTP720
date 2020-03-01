#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:54:41 2020

@author: ryan
"""

import numcalc as nd
import numpy as np

def func(x):
    return -x**2+5

x = list(np.linspace(-4,4, 129))

y = []
for i in x:
    y_val = func(i)
    y.append(y_val)

# testing numerical derivative
y_prime = []
for ind, i in enumerate(x):
    if ind > 0:
        y_new = divtest = nd.numdiv(x,y, i)
        y_prime.append(y_new)


# Integration Testing
    
# midpoint
test1 = nd.midpoint(x,y)

test1(-3.5,3.5) 
# takes more than 2049 points to get up to get 
#5 decimals of precision, checked with integral calculator
# requuires a good deal of resolution to be accurate


# Trapezoid Rule

test2 = nd.traprule(x,y)

test2(-3.5, 3.5)

# Simpson's Rule

test3 = nd.simpson(x,y)

test3(-3.5, 3.5)

