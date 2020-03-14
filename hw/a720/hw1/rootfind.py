#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 14:20:36 2020

@author: ryan
"""

def bisection(function, a, b, threshold = 0.00000001, verbose = False, numiter = False):
    '''
    Simple bisection method for determining roots of a function. The variables a and b represent a lower bound and upper bound
    for the range in which one root is expected to fall. The must be only one root in the given range.
    
    INPUTS:
        function: the function for which to find the roots.
        a: one of the guesses for lower or upper bound.
        b: one of the guesses for lower of upper bound.
        threshold: defines the range plus or minus zero that c and fall in to be considered the root. defaults to 0.00000001 unless specified
        verbose: If True, this option will print out the values for a, b, and c for each iteration.
        numiter: If True, this option will return the number of iterations needed for convergence.
    
    OUTPUT:
        The root of a function, within a specified uncertainy, inside a given range. 
    
    '''
    c = (a+b)/2     # Calculates midpoint between guesses as primer for function
    
    count = 0   # used for counting iterations

    if function(c) == 0:    # Checks to see if you made a lucky guess of bounds
        return(c)
        
    if function(a) * function(b) > 0:   # Exit condition if you choose incorrect bounds
        return(print('Error: root not within bounds'))
        
    elif function(a) * function(b) < 0: # if none of the previous conditions are met, proceeds to looped calculations
        
        while abs(function(c)) > threshold:
            
            count = count + 1   # adds one to count in each iteration
                
            if function(a)*function(c) < 0:     # Checks if a * c is negative, if so a new upper bound b is chosen 
                b = float(c)
                c = (a + b)/2                   # Calculates new midpoint with new upper bound b
                
                # If verbose = True, prints values for each iteration
                    
            elif function(b)*function(c) < 0:   # Checks if b * c is negative, if so a new lower bound a is chosen
                a = float(c)
                c = (a+b)/2
                
            
            if verbose == True:
                print((f"lower = {format(a,'.10f')}, upper = {format(b,'.10f')}, mid = {format(c,'.10f')}"))
                    
            if abs(function(c)) < threshold:    # Exit condition for when the function converges to a value less than the error threshold
                if numiter == True:
                    print(f"NumIter = {count}")
                    return(c, count)
                return(c)
        
        
    
    
def Newton(function, dfunction, x, threshold = 0.00000000000001, verbose = False, numiter = False):
    '''
    Newton's method for finding the nearest root of a given function and its derivative. 
    
    Input
        function: the function for which to find a root
        dfunction: the previous function's derivative
        x: An initial guess for a point near the desired root
        threshold: The amount of acceptable error for finding the root
        verbose: Prints the new position guess for each iteration if True.
        numiter: Prints the number of iterations if True.
        
    Output: 
        The nearest root of the given function.
    '''
    
    x_1 = x - function(x)/dfunction(x) # calculating new coordinate closer to the root
    count = 0                          # Priming the count for number of iterations
    
    while abs(function(x_1)) > threshold: # condition for iteration, while the function evaluates to greater than threshold
        
        count = count + 1
        
        x_1 = float(x_1) - function(float(x_1))/dfunction(float(x_1)) # calculating new coordinate each iteration 
        
        if verbose == True:                     #Prints the new coordinate each iteration if verbose is True
            print(f" Pos {count} = {x_1}")
        
        # Exit condition for when the function at the new 
        #coordinate evaluates to less than the threshold
        if abs(function(x_1)) < threshold:  
            
            if numiter == True:             # print the number of iterations
                print(f"NumIter = {count}")
                return(x_1, count)
                
            return(x_1)

    
    
def secant(function, x_0, x_1, threshold = 0.00000000000001, verbose = False, numiter = False):
    '''
    Secant method for finding nearest root to given coordinates. The must be only one root in the given range.
    
    Input:
        function: the function for which to find the nearest root.
        x_0: the first coordinate near the desired root.
        x_1: the second coordinate near the desired root.
        threshold: The allowed uncertainty on the final answer for the root
        verbose: An option to print out the values for x_0, x_1, and x_2 for each iteration
        numiter: An option to print out the number iterations taken by the function.
        
    Output: 
        The root nearest to the given guesses
    '''
    
    x_2 = x_1 - function(x_1)*((x_1 - x_0)/(function(x_1) - function(x_0))) # equation to find new coordinate
    
    if abs(function(x_2)) == 0:
        return(x_2)
    
    count = 0   # Primer for counting the iterations
    
    while abs(function(x_2)) > threshold: # condition for looping while the function evaluated at x_2 is greater than threshold
        
        count = count + 1 # counting the iterations
        
        x_0 = float(x_1)    # replacing previous values with those from the last iteration
        x_1 = float(x_2)
        
        x_2 =  x_1 - function(x_1)*((x_1 - x_0)/(function(x_1) - function(x_0))) # calculating new x_2 with new X_0 and x_1
        
        if verbose == True:
            print(f" Pos {count} ~ ~ x_0 = {x_0}, x_1 = {x_1}, x_2 = {x_2}, Func(x_2) = {function(x_2)}")
        
        if abs(function(x_2)) < threshold:
            
            if numiter == True:
                print(f"NumIter = {count}")
                return(x_2, count)
                
            return(x_2)
    
    