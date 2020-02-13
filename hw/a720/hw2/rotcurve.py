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

G = 4.299e-6*u.kpc/u.second**2*u.km**2/u.solMass

def rotcurve(r_200,x, v_200, c, velfunc, savefig = False):
    
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
        plt.figure(f'{i}')
        plt.xlabel('R (kpc)')
        plt.ylabel('Mass (solar mass)')
        
        for bind, j in enumerate(v_200):
           radius = x*r_200
           velocity = velfunc(x,i,j)
           mass = massfunc(radius, velocity, G)
           
           plt.plot(radius, mass, label = f'V_200 = {j}')
           plt.title(f'c = {i}')
           plt.legend()
           
           m_r = []
           radius2 = []
           
           for ind, mass_enc in enumerate(mass):
              
               if ind > 0:
                   m = (mass[ind]-mass[ind -1])/(u.solMass)
                   m_r.append(m)
                   r = ((radius[ind]+radius[ind])/2)/(u.kpc)
                   radius2.append(r)
                   
           plt.figure()
           plt.xlabel('R (kpc)')
           plt.ylabel('M(r) (solar mass)')
           plt.plot(radius2, m_r/s, label = f'V_200 = {j}')
           plt.title(f'c = {i}')
           plt.legend()
           
           if savefig == True:
               plt.savefig(f'massenc_{ind}_{bind}.pdf')
        return(m_r, radius2)


v_200 = [150,200,250,300]*u.km/u.second
c = [5, 20]
r_200 = 200*u.kpc
x = np.linspace(0.0001, 2, 1000)

massenc(r_200, x, v_200, c, velfunc, massfunc)

rotcurve(r_200, x, v_200, c, velfunc)
