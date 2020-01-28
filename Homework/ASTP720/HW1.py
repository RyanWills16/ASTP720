#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:13:24 2020

@author: ryan
"""

import rootfind as rt


def func(x):
    return(x**2 - 45)
    
def dfunc(x):
    return(2*x)


rt.bisection(func,4,11, 0.0000000001, verbose = True, numiter = True)



rt.Newton(func, dfunc, 1, 0.00000000000001, verbose = True, numiter = True)

rt.secant(func, -2, -1, verbose = True, numiter = True)
