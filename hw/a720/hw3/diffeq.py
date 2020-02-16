# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 15:31:35 2020

@author: Ryan
"""

def euler(func, initial, time, step = None):
    
    if not step == None:
        import numpy as np
        h = step
        n = ((time[-1] - time[0])/h + 1)
        time = np.linspace(time[0], time[-1],n)
        
    else:
        h = (time[-1] - time[0])/(len(time)-1)
    
    
    f = []
    y_prime = []
    f.append(initial[0])
    y_prime.append(initial[1])
    
    for ind, t in enumerate(time):
        if ind > 0:
            y_prime_new = y_prime[ind-1] + h*func(time[ind-1],f[ind-1],y_prime[ind-1])[1]
            
            y_prime.append(y_prime_new)
            
            y_new = f[ind-1] + h*func(time[ind-1],f[ind-1],y_prime[ind])[0]
            
            
            f.append(y_new)
    return f,y_prime, time

def heun(func, initial, time, step = None):
    
    if not step == None:
        import numpy as np
        h = step
        n = (time[-1] - time[0]/h + 1)
        time = np.linspace(time[0], time[-1],n)
        
    else:
        h = (time[-1] - time[0])/(len(time)-1)
        
    f = []
    y_prime = []
    f.append(initial[0])
    y_prime.append(initial[1])
    for ind, t in enumerate(time):
        if ind > 0:
            
            f_dzdt = func(time[ind -1], f[ind-1], y_prime[ind-1])[1]
            f_z = func(time[ind-1], f[ind-1],  y_prime[ind-1])[0]
            
            f_dzdt2 = func(t, f[ind-1] + f_z*h, y_prime[ind-1]+f_dzdt*h)[1]
            f_z2 = func(t,f[ind-1] + f_z*h,y_prime[ind-1]+f_dzdt*h)[0]
            
            dzdt = y_prime[ind-1] + 0.5*h*(f_dzdt + f_dzdt2)
            z = f[ind -1] + 0.5*h*(f_z + f_z2)
            
            y_prime.append(dzdt)
            f.append(z)
            
            
    return f,y_prime, time