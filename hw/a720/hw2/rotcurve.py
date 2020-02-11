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
    return r*velocity**2/G

def densityfunc(ro_0,r_200,c )
G = 4.299e-6*u.kpc/u.second**2*u.km**2/u.solMass
c = 15
v_200 = 200*u.km/u.second
r_200 = 250*u.kpc


x = np.linspace(0.0001,2,1000)
radius = x*r_200
y = np.linspace(5,995,199)


velocity = velfunc(x, c, v_200)
mass = massfunc(radius, velocity, G)

plt.figure()
plt.plot(radius, velocity)
plt.xlabel('r')
plt.ylabel('velocity (km/s)')
                
plt.figure()
plt.plot(radius, mass)
plt.xlabel('radius')
plt.ylabel('Mass enclosed (solar masses)')    

def rotcurve(r_200,x, v_200, c, velfunc, savefig = False):
    
    G = 4.299e-6*u.kpc*u.km**2/u.solMass
    
    for ind, i in enumerate(c): 
        plt.figure()
        
        plt.xlabel('R (kpc)')
        plt.ylabel('Mass (solar mass)')
        for bind, j in enumerate(v_200):
           radius = x*r_200
           velocity = velfunc(x,i,j)
           
           plt.plot(radius, velocity, label = f'V_200 = {j}')
           plt.title(f'c = {i}')
           plt.legend()
           if savefig == True:
               plt.savefig(f'rotcurve_{ind}_{bind}.pdf')
           
def massenc(r_200,x, v_200, c, velfunc, massfunc, savefig = False):
    G = 4.299e-6*u.kpc*u.km**2/u.solMass
    
    for ind, i in enumerate(c): 
        plt.figure()
        
        plt.xlabel('R (kpc)')
        plt.ylabel('Mass (solar mass)')
        for bind, j in enumerate(v_200):
           radius = x*r_200
           velocity = velfunc(x,i,j)
           mass = massfunc(radius, velocity, G)
           
           plt.plot(radius, mass, label = f'V_200 = {j}')
           plt.title(f'c = {i}')
           plt.legend()
           if savefig == True:
               plt.savefig(f'massenc_{ind}_{bind}.pdf')
           
    
v_200 = [150,200,250,300]*u.km/u.second
c = [0.5, 10, 50]
r_200 = 200*u.kpc
x = np.linspace(0.0001, 2, 1000)

massenc(r_200, x, v_200, c, velfunc, massfunc)

rotcurve(r_200, x, v_200, c, velfunc)
