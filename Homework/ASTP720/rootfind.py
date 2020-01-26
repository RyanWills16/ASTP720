#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 14:20:36 2020

@author: ryan
"""

def bisection(function, a, b, threshold = 0.00000001, verbose = False):
    '''
    Simple bisection method for determining roots of a function. The variables a and b represent a lower bound and upper bound
    for the range in which one root is expected to fall. 
    
    INPUTS:
        function: the function for which to find the roots.
        a: one of the guesses for lower or upper bound.
        b: one of the guesses for lower of upper bound.
        threshold: defines the range plus or minus zero that c and fall in to be considered the root. defaults to 0.00000001 unless specified
        verbose: If True, this option will print out the values for a, b, and c for each iteration.
    
    OUTPUT:
        The root of a function, within a specified uncertainy, inside a given range. 
    
    '''
    c = (a+b)/2     # Calculates midpoint between guesses

    if function(c) == 0:    # Checks to see if you made a lucky guess of bounds
        return(c)
        
    if function(a) * function(b) > 0:   # Exit condition if you choose incorrect bounds
        return(print('Error: root not within bounds'))
        
    elif function(a) * function(b) < 0: # if none of the previous conditions are met, proceeds to looped calculations
        
        # added condition to account for any situation in which func(a) is positive and func(b) negative 
        # or func(b) is positive and func(a) is negative. Might help for retrieving negative roots
        if function(a) < 0 and function(b) > 0: 
            
            while function(a) < 0 and function(b) > 0:
                
                if function(a)*function(c) < 0:     # Checks if a * c is negative, if so a new upper bound b is chosen 
                    b = float(c)
                    c = (a + b)/2                   # Calculates new midpoint with new upper bound b
                    
                    # If verbose = True, prints values for each iteration
                    if verbose == True:
                        print('a * b < 0')
                        print((f"lower = {format(a,'.6f')}, upper = {format(b,'.6f')}, mid = {format(c,'.6f')}")) 
                        
                elif function(b)*function(c) < 0:   # Checks if B * c is negative, if so a new lower bound a is chosen
                    a = float(c)
                    c = (a+b)/2
                    
                    # If verbose = True, prints values for each iteration
                    if verbose == True:
                        print((f"lower = {format(a,'.6f')}, upper = {format(b,'.6f')}, mid = {format(c,'.6f')}"))
                        
                if abs(function(c)) < threshold:    # Exit condition for when the function converges to a value less than the error threshold
                    return(c)
                        
        # condition for when func(a) is positive and func(b) is negative
        elif function(a) > 0 and function(b) < 0:
                
            while function(a) > 0 and function(b) < 0: # revers condition from above while loop
                
                if function(a)*function(c) < 0:    
                    b = float(c)
                    c = (a + b)/2                   
    
                    if verbose == True:
                        print('a * b < 0')
                        print((f"lower = {format(a,'.6f')}, upper = {format(b,'.6f')}, mid = {format(c,'.6f')}")) 
                        
                elif function(b)*function(c) < 0:   
                    a = float(c)
                    c = (a+b)/2
                    
                    
                    if verbose == True:
                        print((f"lower = {format(a,'.6f')}, upper = {format(b,'.6f')}, mid = {format(c,'.6f')}"))
                
                if abs(function(c)) < threshold:    # Exit condition for when the function converges to a value less than the error threshold
                    return(c)
                    
    
    
    
    
    
    
    
    
#def Newton():
#    #not putting anything here yet
#    
#
#    
#    
#def secant():
#    nothing his here yet
#    
#    