#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 19:44:31 2020

@author: ryan
"""


def linterp(x, y):
    
    '''
    INPUTS:
        x: list or array of x coordinates
        y: list or array of y coordinates cooresponding to x
        
    OUTPUT:
        The function calcpoint, which outputs a coordinate.
    '''
         
    def calcpoint(x_p, verbose = False): # Function that takes x_p as input, which is the x coordinate of the desired y coordinate
        
        '''
        INPUTS:
            x_p: the desired x point at which to find the coresponding y point using linear interpolation. 
            
        Output:
            
            The y coordinate corresponding to the input x coordinate
             
            Also prints extra information like the points used, the slope, and intercept for the interpolated line if the option verbose is True. 
        '''
        
        for ind, i in enumerate(x): # iterating through the list, keeping trak of indicy
            
            
            if x[ind-1] < x_p < x[ind] and y[ind-1] == y[ind]: # accounting for the possibility of a horizontal line between two points
                return(y[ind])
                
                
            if x[ind-1] <= x_p <= x[ind]: # if the given point falls between two x coordinates, use those in the calculation
                
                # equation for finding the y coordinate of the given x coordinate, x_d is subbed into the equation already.
                y_p = y[ind-1] * (1 - ((x_p - x[ind-1])/(x[ind] - x[ind-1]))) + y[ind] * ((x_p - x[ind-1])/(x[ind] - x[ind-1]))
                
                if verbose == True: # do some extra calculations for slope and intercept and print out a bunch of info.
                    
                    slope = (y[ind] - y[ind-1]) / (x[ind] - x[ind-1]) # calculate slope of line between two points
                
                    b = y[ind] - slope * x[ind] # calculate 
                
                    return(f"point = {y_p}, x_1 = {x[ind]}, x_0 = {x[ind-1]}, y_1 = {y[ind]}, y_0 = {y[ind-1]}, slope = {slope}, intercept = {b}")
                    
                elif verbose == False:
                    return(y_p)
                
    return calcpoint # returns the callable function



def splineterp(x, y):
    
    '''
    '''
    