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
        The derivative at x_p
    '''
    
    if x_p == x[0] or x_p == x[-1]:
        return(print("cannot find derivative at the first or last point"))
            
    for ind, i in enumerate(x):
            
        if x[ind] == x_p:
            
            h = x[ind+1] - x_p # calculate width
            y = (y[ind +1] - y[ind -1])/(2*h)   # calculate derivative
            return(y)
            
def midpoint(x,y):
    '''
    Input:
        x: list of x points
        y: list of corresponding y points
        
    Output:
        helper function calculate is output
    '''
    import interpolation as ip
    # using linear interpolation to find point between x_1 and x_2
    p = ip.linterp(x,y) 
    def calculate(a, b, verbose = False):
        '''
        Input:
            a: lower bound for integration
            b: upper bound for integration
        Output:
            The area under the curve over the range of (a,b) calculated via 
            midpoint rule
        '''
        
        area = 0
        
        for ind, i in enumerate(x):
            
            if a <= i < b:
                
                h = abs(x[ind+1] - i)
                x_point = ((i + x[ind+1])/2) # finding midpoint
                y_p = p(x_point)    # finding y_p at x_p
                
                A = y_p*h # calculate area
                area = area + A # add areas up
                
                if verbose == True:
                    print(f"width = {h}, x_point = {x_point}, y_p = {y_p}, Area = {A}")
        
            
        return area
    return calculate

def traprule(x,y):
    '''
    Input:
        x: list of x points
        y: list of corresponding y points
        
    Output:
        helper function calculate is output
    '''
    def calculate(a, b, verbose = False):
        '''
        Input:
            a: lower bound for integration
            b: upper bound for integration
        Output:
            The area under the curve over the range of (a,b) calculated via 
            trapezoid rule
        '''
        area = 0
        
        for ind, i in enumerate(x):
            
            if a <= i < b:
                
                h = abs(x[ind+1] - i)   # calculate width
                y_p = (y[ind] + y[ind+1])/2 # calculate y
                
                A = y_p*h  # calculate area
                area = area + A # add areas up
                
                if verbose == True:
                    print(f"indicy = {ind}, x = {i}, x[ind+1] = {x[ind+1]}")
                    print(f"y[ind] = {y[ind]}, y[ind + 1] = {y[ind+1]}")
                    print(f"width = {h}, y_p = {y_p}, Area = {A}")
        
            
        return area
    return calculate

def simpson(x,y):
    '''
    Input:
        x: list of x points
        y: list of corresponding y points
        
    Output:
        helper function calculate is output
    '''
    def calculate(a, b, verbose = False):
        '''
        Input:
            a: lower bound for integration
            b: upper bound for integration
        Output:
            The area under the curve over the range of (a,b) calculated via 
            Simpson's rule
        '''
        
        area = 0
            
        x_new = []
        y_new = []
        
        # appending the numbers in the range a:b to new lists
        for ind, j in enumerate(x):
            if a <= j <= b:
                x_new.append(j)
                y_new.append(y[ind])
                
        length = len(x_new) # find length of new x list
        
        h = abs((b - a)/(length - 1))   # calculate width based on lenght

        
        for ind, i in enumerate(x_new):
            
            # condition for summation
            if 0 <= ind <= (length-1)/2 -1:
                
                # calculate area
                A = (h/3) * (y_new[2*ind] + 4*y_new[2*ind +1] + y_new[2*ind + 2])
                
                # add areas up
                area = area + A
                
                # option to get more information about the workings
                if verbose == True:
                    print('first',ind,'//', h,'//', x_new[ind+1],'//', i ,'//', area)
                    print('second','//',y_new[2*ind],'//', y_new[2*ind +1],'//', y_new[2*ind + 2])
                
        return area
    
    return calculate
            