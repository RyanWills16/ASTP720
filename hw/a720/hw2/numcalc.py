#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:27:04 2020

@author: ryan
"""



def numdiv(x, y, x_p):
    
    if x_p == x[0] or x_p == x[-1]:
        return(print("cannot find derivative at the first or last point"))
            
    for ind, i in enumerate(x):
            
        if x[ind] == x_p:
            
            h = x[ind+1] - x_p
            y = (y[ind +1] - y[ind -1])/(2*h)
            return(y)
            
def midpoint(x,y):
    import interpolation as ip
    p = ip.linterp(x,y)
    def calculate(a, b, verbose = False):
        
        area = 0
        
        for ind, i in enumerate(x):
            
            if a <= i < b:
                
                h = abs(x[ind+1] - i)
                x_point = ((i + x[ind+1])/2)
                y_p = p(x_point)
                
                A = y_p*h
                area = area + A
                
                if verbose == True:
                    print(width = f'{h}', x_point = f'{x_point}', y_p = f'{y_p}', Area = f'{A}')
        
            
        return area
    return calculate

def traprule(x,y):
    def calculate(a, b, verbose = False):
        
        area = 0
        
        for ind, i in enumerate(x):
            
            if a <= i < b:
                
                h = abs(x[ind+1] - i)
                y_p = (y[ind] + y[ind+1])/2
                A = y_p*h
                area = area + A
                
                if verbose == True:
                    print(width = f'{h}', y_p = f'{y_p}', Area = f'{A}')
        
            
        return area
    return calculate