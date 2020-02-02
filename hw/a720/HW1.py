#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 15:13:24 2020

@author: ryan
"""

import rootfind as rt
import interpolation as ip
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def func(x):
    return(x**2 - 45)
    
def dfunc(x):
    return(2*x)


rt.bisection(func,4,11, 0.0000000001, numiter = True)

rt.Newton(func, dfunc, 1, 0.00000000000001, numiter = True)

rt.secant(func, -2, -1, numiter = True)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

x1 = np.arange(-4.5,5, 0.5) # randomly chosen range of numbers

y1 =  2*x1 + 11 -x1**2 # randomly chosen quadratic function to get a set of y values


p = ip.linterp(x1,y1) # call the linear interpolation function to prime it with sets of points x and y

p(-4.5, verbose = True) # interpolate a random point to check that the function works


x1_2 = [] # defining a new empty set of x values
y1_2 = [] # defining a new empty set of y values

for i in np.arange(-4.5,x1[-1],0.25):    # iterating over a range of new x points at a higher resolution
    newpoint = ip.linterp(x1,y1)  # priming the function with the original x, y data points
    y_new = newpoint(i)         # using linear interpolation to find new x data points
    
    x1_2.append(i)        # appending new points to empty array
    y1_2.append(y_new)    # appending new y points to empty array
    
plt.figure()
    
plt.scatter(x1_2,y1_2, marker = 'o', edgecolors = 'r', facecolors = 'none', label = 'interpolated points') # plot new interpolated points

plt.scatter(x1,y1, marker = 'o', edgecolors = 'b', facecolors = 'none', label = 'Original points') # plot original points in a different color

plt.title('11 + 2x - x^2')
plt.legend()

plt.savefig('interp.pdf')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def func2(x):
    return(((1+x**2)**(-0.5) - 0.5))
    
def dfunc2(x):
    return(-1*x*(1+x**2)**-1.5)
    
thresholds = np.e**np.linspace(-30,-6, 15)

# defining empty arrays for the points generated below
points_bi =[]
points_newton = []
points_secant = []

iterns_bi = []
iterns_newton = []
iterns_secant = []

# generating points for the number of iterations based on in input list of threshold values.
for ind, i in enumerate(thresholds):
    fwhm_bi, itern_bi = rt.bisection(func2, 1, 4, threshold = i, numiter = True) # gettin point and number of iteration for bisect
    fwhm_newton, itern_newton = rt.Newton(func2, dfunc2, 1.5, threshold = i, numiter = True)    # gettin point and number of iteration for bisect
    fwhm_secant, itern_secant = rt.secant(func2, 1,3, threshold = i, numiter = True)    # gettin point and number of iteration for bisect
    
    # appending the points from each method to the empty arrays
    points_bi.append(fwhm_bi)
    iterns_bi.append(itern_bi)
    
    points_newton.append(fwhm_newton)
    iterns_newton.append(itern_newton)
    
    points_secant.append(fwhm_secant)
    iterns_secant.append(itern_secant)

# plotting and saving the figure for number of iteration vs threshold
plt.figure()

# plotting thresholds on x-axis and number of iterations on y-axis
plt.plot(np.log(thresholds), iterns_bi, color = 'r', label = "Bisection method")
plt.plot(np.log(thresholds), iterns_newton, color = 'b',  label = "Newton's method")
plt.plot(np.log(thresholds), iterns_secant, color = 'g', label = "Secant Method")

ticks = np.linspace(-30, -5, 6)    # generating x values for tick marks
plt.xticks(ticks, fontsize = 11)    

plt.xlabel('ln(Threshold)', fontname = 'Times New Roman', fontsize = 16)
plt.ylabel('Number of Iterations', fontname = 'Times New Roman', fontsize = 16)

plt.yticks(fontsize = 11)
plt.tight_layout()
plt.legend(fontsize = 12)


plt.savefig('threshvsiter.pdf')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# defining constants
N_0 = 0.01*206265*u.au*u.cm**-3
wv = 21*u.cm 
r_e = 2.817 * 10**-13*u.cm
D = 1*u.kpc*2.063*10**8*(u.au/u.kpc)
a = 1*u.au

def gausstrace(x):
        return x*(1 + Q*np.exp((-(x)**2))) - x_prime # output will be in AU since Q is unitless and x_prime is in AU
    
    
Q = wv**2*r_e*N_0*D*np.pi**-1*a**-2 # the coefficient in front of the guassian part of the coordinate transformation

xprime = list(np.arange(0.1,2.1,0.2))   # create list of xprime values from 0 to 2 AU

xpoints = []    # create an empty list for xvalues at the lense

for i in xprime:    # iterate through xprime values
    print('x_prime =', i)
    x_prime = i
    x_lens = rt.secant(gausstrace,-1, 2.5, threshold = 0.000000000001, verbose = True) # calculate root
    xpoints.append(x_lens)  # append x values to list

# create figure of ray tracing by drawing lines between (x_prime,0) and (x,D)
plt.figure(figsize = (8,5.5))
plt.tight_layout()
plt.xlabel('x_prime (AU)', fontsize = 14)
plt.ylabel('Distance (AU)', fontsize = 14)
plt.title('Gaussian', fontsize = 14)
for ind, j in enumerate(xpoints):   # iterating over points to add one line each iteration
    point1 = [xprime[ind], 0]
    point2 = [j, D/u.au]
    
    x_val = [point1[0], point2[0]]
    y_val = [point1[1], point2[1]]
    
    
    plt.plot(x_val, y_val, color = 'r')
plt.savefig('question4.pdf')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBLEM 5
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
N_0 = 0.01*206265*u.au*u.cm**-3
wv = 21*u.cm 
r_e = 2.817 * 10**-13*u.cm
D = 1*u.kpc*2.063*10**8*(u.au/u.kpc)
r_c = 1*u.au


def psuedoiso(x):
        return x*(1 + R*(1 + x**2)**-(3/2)) - x_prime2 # output will be in AU since R is 
                                                       #unitless and x_prime2 is in AU
    
    
R = 0.5*wv**2*r_e*N_0*D*np.pi**-1*r_c**-2 # the coefficient calculation

xprime2 = list(np.arange(0,2.1,0.15))   # create list of xprime values from 0 to 2 AU

xpoints2 = []    # create an empty list for xvalues at the lense

for k in xprime2:   # iterate through given xprime2 values
    print('x_prime =', k)
    x_prime2 = k
    x_lens2 = rt.secant(psuedoiso,-1, 2.5, threshold = 0.000000000001, verbose = True)  # find root
    xpoints2.append(x_lens2) # append to empty list


# plot the ray tracing 
plt.figure(figsize = (8,5.5))
plt.tight_layout()
plt.xlabel('x_prime (AU)', fontsize = 14)
plt.ylabel('Distance (AU)', fontsize = 14)
plt.title('Psuedo - Isothermal Spherical', fontsize = 14)
for ind, l in enumerate(xpoints2):  # iterate though points to add one line each iteration
    point3 = [xprime2[ind], 0]
    point4 = [l, D/u.au]
    
    x_val2 = [point3[0], point4[0]]
    y_val2 = [point3[1], point4[1]]
    
    
    plt.plot(x_val2, y_val2, color = 'r')
plt.savefig('question5.pdf')

    







