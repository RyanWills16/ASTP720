#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:13:24 2020

@author: ryan
"""

import rootfind as rt


def func(x):
    return(x**2 - 45)


rt.bisection(func,-11,-4, 0.00000000000001, numiter = True)
