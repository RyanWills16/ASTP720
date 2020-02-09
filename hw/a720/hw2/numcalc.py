#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:27:04 2020

@author: ryan
"""



def numdiv(x, y, x_p):
    
    '''
    Input:
        x: list of x data points
        y: list of y data points
        x_p: the point at which the derivative is to be found
        
    Returns:
        The derivative 
    '''
    
    if x_p == x[0] or x_p == x[-1]:
        return(print("cannot find derivative at the first or last point"))
            
    for ind, i in enumerate(x):
            
        if x[ind] == x_p:
            
            h = x[ind+1] - x_p
            y = (y[ind +1] - y[ind -1])/(2*h)
            return(y)
            
def midpoint(x,y):
    '''
    '''
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
                    print(f"width = {h}, x_point = {x_point}, y_p = {y_p}, Area = {A}")
        
            
        return area
    return calculate

def traprule(x,y):
    '''
    '''
    def calculate(a, b, verbose = False):
        
        area = 0
        
        for ind, i in enumerate(x):
            
            if a <= i < b:
                
                h = abs(x[ind+1] - i)
                y_p = (y[ind] + y[ind+1])/2
                
                A = y_p*h
                area = area + A
                
                if verbose == True:
                    print(f"indicy = {ind}, x = {i}, x[ind+1] = {x[ind+1]}")
                    print(f"y[ind] = {y[ind]}, y[ind + 1] = {y[ind+1]}")
                    print(f"width = {h}, y_p = {y_p}, Area = {A}")
        
            
        return area
    return calculate

def simpson(x,y):
    '''
    '''
    def calculate(a, b, verbose = False):
        
        area = 0
            
        x_new = []
        y_new = []
        for ind, j in enumerate(x):
            if a <= j <= b:
                x_new.append(j)
                y_new.append(y[ind])
                
        length = len(x_new)
        
        h = abs((b - a)/(length - 1))

        
        for ind, i in enumerate(x_new):
            
            if 0 <= ind <= (length-1)/2 -1:
                
                
                
                A = (h/3) * (y_new[2*ind] + 4*y_new[2*ind +1] + y_new[2*ind + 2])
                
                area = area + A
                
                if verbose == True:
                    print('first',ind,'//', h,'//', x_new[ind+1],'//', i ,'//', area)
                    print('second','//',y_new[2*ind],'//', y_new[2*ind +1],'//', y_new[2*ind + 2])
                
        return area
    
    return calculate
            