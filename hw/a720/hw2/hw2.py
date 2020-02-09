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

x = list(np.linspace(-4,4, 17))

y_val = []
for i in x:
    y = func(i)
    y_val.append(y)


# Integration Testing
x2 = nd.midpoint(x,y_val)

x2(-3.4,3.4) 
# takes more than 2049 points to get up to get 
#5 decimals of precision, checked with integral calculator
# requuires a good deal of precision to be accurate


# Trapezoid Rule

x3 = nd.traprule(x,y)

x3(-3.5, 3.5)