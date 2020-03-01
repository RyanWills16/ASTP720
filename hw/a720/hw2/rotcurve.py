# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 22:59:56 2020

@author: Ryan
"""

import matplotlib.pyplot as plt
import numpy as np
import numcalc as nm
from astropy import units as u

def velfunc(x, c, v_200):
    '''
    returns function for v_c^2
    '''
    return np.sqrt((v_200**2/x)*((np.log(1 + c*x) - (c*x)/(1 + c*x))/(np.log(1 + c) - (c)/(1 + c))))

def massfunc(r, velocity, G):
    '''
    returns function for mass enclosed
    '''
    return r*velocity**2/G


def rotcurve(r_200,x, v_200, c, velfunc, savefig = False):
    '''
    Input:
        r_200: radius at which density is 200 times the critical density
        x: a list of the parameterized r/r_200 values
        v_200: the velocity at 200 kpc
        c: the coefficient in the velocity equation
        velfunc: the function for calculating velocity
    Output:
        plots of the rotation curve at different c and v_200 values
    '''
    
    for ind, i in enumerate(c): 
        
        plt.figure()
        plt.xlabel('R (kpc)')
        plt.ylabel('Circular Velocity (km/s)')
        
        for bind, j in enumerate(v_200):
           radius = x*r_200 # calculate radii from parameter x and r_200
           velocity = velfunc(x,i,j) # calculate the velocity distribution
           
           plt.plot(radius, velocity, label = f'V_200 = {j}')
           plt.title(f'c = {i}')
           plt.legend()
           
        if savefig == True:
            plt.savefig(f'rotcurve_{i}.pdf')
    
    
def massenc(r_200,x, v_200, c, velfunc, massfunc, savefig = False):
    '''
    Input:
        r_200: radius at which density is 200 times the critical density
        x: a list of the parameterized r/r_200 values
        v_200: the velocity at 200 kpc
        c: the coefficient in the velocity equation
        velfunc: the function for calculating velocity
        massfunc: the function for calculating mass enclosed
        
    Output:
        plots of the rotation curve at different c and v_200 values
    '''
    
    # gravitation constant in kpc/s^2 km^2/solarMass
    G = 4.299e-6*u.kpc/u.second**2*u.km**2/u.solMass

    for ind, i in enumerate(c): # iterate through c values
        plt.figure()    # prime the figure for plotting
        plt.xlabel('R (kpc)')
        plt.ylabel('Mass (solar mass)')

        for bind, j in enumerate(v_200): # iterate through v_200 values

            radius = x*r_200 # calculate radii from parameter x and r_200
            velocity = velfunc(x,i,j) # calculate the velocity distribution
            mass = massfunc(radius, velocity, G) # calculate mass enclosed
   
   
            plt.plot(radius, mass, label = f'V_200 = {j}')
            plt.title(f'c = {i}')
            plt.legend()
               
        if savefig == True:
            plt.savefig(f'massenc_{i}.pdf')
                
def calc_mass(r_200,x, v_200, c, velfunc, massfunc, savefig = False):
    G = 4.299e-6*u.kpc/u.second**2*u.km**2/u.solMass

    for ind, i in enumerate(c):
        plt.figure()
        plt.xlabel('R (kpc)')
        plt.ylabel('M(r)')

        for bind, j in enumerate(v_200):

            radius = x*r_200 # calculate radii from parameter x and r_200
            velocity = velfunc(x,i,j) # calculate the velocity distribution
            mass = massfunc(radius, velocity, G) # calculate mass enclosed
            
            radius2 = []    # empty list for new radii
            mass2 = []  # empty list for new masses
    
            for ind, mass_enc in enumerate(mass):
                
                # calculate mass at each radius
                # define new radius as averages of two radii
                if ind > 0:
                    m = mass[ind]/(u.solMass)-mass[ind -1]/(u.solMass)
                    r = (radius[ind]/(u.kpc)+radius[ind-1]/(u.kpc))/2
                       
                    mass2.append(m)
                    radius2.append(r)
                    
            plt.plot(radius2, mass2, label = f'V_200 = {j}')
            plt.title(f'c = {i}')
            plt.legend()
               
        if savefig == True:
            plt.savefig(f'mass_profile_{i}.pdf')
                
def mass_grad(r_200,x, v_200, c, velfunc, massfunc, savefig = False):
    G = 4.299e-6*u.kpc/u.second**2*u.km**2/u.solMass

    for ind, i in enumerate(c):
        plt.figure()
        plt.xlabel('R (kpc)')
        plt.ylabel('dM(r)/dr')

        for bind, j in enumerate(v_200):

            radius = x*r_200 # calculate radii from parameter x and r_200
            velocity = velfunc(x,i,j) # calculate the velocity distribution
            mass = massfunc(radius, velocity, G) # calculate mass enclosed
            
            radius2 = []    # empty list for new radii
            mass2 = []  # empty list for new masses
    
            for ind, mass_enc in enumerate(mass):
                
                # calculate mass at each radius
                # define new radius as averages of two radii
                if ind > 0:
                    m = mass[ind]/(u.solMass)-mass[ind -1]/(u.solMass)
                    r = (radius[ind]/(u.kpc)+radius[ind-1]/(u.kpc))/2
                       
                    mass2.append(m)
                    radius2.append(r)
                    
            dM_r = []   # empty list for derivative of mass dist
            radius3 = []    # empty list for radii
            
            # iterate over values in radius2
            for ind, k in enumerate(radius2):
                
                # condition so we don't use the first or last point
                if radius2[0] < k < radius2[-1]:
                    y = nm.numdiv(radius2,mass2, k) # calculating derivative
                    radius3.append(k)
                    dM_r.append(y)
                    
            plt.plot(radius3, dM_r, label = f'V_200 = {j}')
            plt.title(f'c = {i}')
            plt.legend()
               
        if savefig == True:
            plt.savefig(f'mass_grad_{i}.pdf')


# defining constants for use in the above functions
v_200 = [150, 250, 300]*u.km/u.second
c = [5, 20]
r_200 = 200*u.kpc
x = np.linspace(0.0001, 2, 1000)

# creating plot of mass enclosed
massenc(r_200, x, v_200, c, velfunc, massfunc, savefig = True)

# creating plot of mass distribution
calc_mass(r_200, x, v_200, c, velfunc, massfunc, savefig = True)

# creating plot of rotation curve
rotcurve(r_200, x, v_200, c, velfunc, savefig = True)

# creating plot of mass gradient
#mass_grad(r_200, x, v_200, c, velfunc, massfunc, savefig = True)
