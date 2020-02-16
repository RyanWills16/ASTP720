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
            y_prime_new = y_prime[ind-1] + h*func(time[ind-1],f[ind-1],y_prime[ind-1])
            
            y_new = f[ind-1] + h*y_prime_new
            
            y_prime.append(y_prime_new)
            f.append(y_new)
    return f,y_prime, time

def heun(func, initial, time):
    
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
            
            step = h*func(time[ind-1],f[ind-1],y_prime[ind-1])
            yp1 = y_prime[ind-1] + step
            print('yp1 = ',yp1)
            
            step2 = 0.5*h*(func(time[ind-1],f[ind-1],y_prime[ind-1]) + func(t,f[ind-1],yp1))
            yp2 = y_prime[ind-1] + step2
            print('yp2 = ', yp2)
            
            y1 = f[ind-1] + h*yp1
            print('y1 = ', y1)
            y2 = f[ind-1] + 0.5*h*(y1 + yp2)
            print('y2 = ', y2)
            y_prime.append(yp1)
            f.append(y1)
    return f,y_prime