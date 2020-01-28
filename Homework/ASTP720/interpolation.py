#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 19:44:31 2020

@author: ryan
"""


def linterp(x, y):
    
    '''
    '''
         
    def calcpoint(x_p): # Function that takes x_p as input, which is the x coordinate of the desired y coordinate
        
        '''
        '''
        
        for ind, i in enumerate(x): # iterating through the list, keeping trak of indicy
            
            if x[ind-1] < x_p < x[ind] and y[ind-1] == y[ind]: # accounting for the possibility of a horizontal line between two points
                
                return(y[ind])
                
                
            if x[ind-1] < x_p < x[ind]: # if the given point falls between two x coordinates, use those in the calculation
                
                # equation for finding the y coordinate of the given x coordinate, x_d is subbed into the equation already.
                y_p = y[ind-1] * (1 - ((x_p - x[ind-1])/(x[ind] - x[ind-1]))) + y[ind] * ((x_p - x[ind-1])/(x[ind] - x[ind-1]))
                
                return(y_p)
                
    return calcpoint # returns the callable function


        
        
