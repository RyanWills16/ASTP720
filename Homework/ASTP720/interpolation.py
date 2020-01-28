#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 19:44:31 2020

@author: ryan
"""


x = [1,2,5,8,10,11]
y = [1,3,7,9,13,14]

def linpiecewise(x, y):
    
    for ind, i in enumerate(x):
        if x[ind-1] < x[ind]:
            
            def calc_point(x_p):
                if x[ind-1] < x_p < x[ind]:
                    return(float(y[ind-1] * (1 - ((x_p - x[ind-1])/(x[ind] - x[ind-1]))) + y[ind] * ((x_p - x[ind-1])/(x[ind] - x[ind-1]))))
            return calc_point
        
numbers = linpiecewise(x, y)
x2 = numbers(10.5)         

x_p = 10.5       
for ind, i in enumerate(x):
    if x[ind-1] < x[ind]:
        if x[ind-1] < x_p < x[ind]:
            x_new = y[ind-1] * (1 - ((x_p - x[ind-1])/(x[ind] - x[ind-1]))) + y[ind] * ((x_p - x[ind-1])/(x[ind] - x[ind-1]))
        
        
