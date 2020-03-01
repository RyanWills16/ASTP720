# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 15:31:35 2020

@author: Ryan
"""

def euler(func, initial, time, step = None):
    '''
    Takes a list of one or two funciton and a list of one or two intial conditions
    
    Input:
        func: a function of two variables used for the calculation, obtained by setting 
            z = dydt (func[0]) and its derivative dzdt = d2ydt2 (func[1])
        initial: the initial conditions y = [y(0),dydt(0)]
        time: either the full series, or stop and start times if step is defined
        step: the amount by which to step the calculation between t[0] and t[1]
    Output
        list of values pertaining to the analytical solution for the second
        order differential equation input
    '''
    
    # if step is defined, build time series based on that step
    # and the given time range
    if not step == None:
        import numpy as np
        h = step
        n = ((time[-1] - time[0])/h + 1)
        time = np.linspace(time[0], time[-1],n)
        
    # otherwise, just use the step derived from the time range
    else:
        h = (time[-1] - time[0])/(len(time)-1)
    
    # create empty lists for vales of y and y'
    f = []
    y_prime = []
    
    if len(initial) == 1:
        f.append(initial[0])
        
        for ind, t in enumerate(time):
            
            if ind > 0:
                
                # Using the Rung-Katta 2 method
                # The first calculation: predictor
                f_z = func(time[ind-1], f[ind-1])[0]
                z = f[ind-1] + f_z*h
                # appending to empty lists
                f.append(z)
        return f, time
    
    else: 
        f.append(initial[0])
        y_prime.append(initial[1])
        # iterate through the time series
        for ind, t in enumerate(time):
            
            # start with the second value in the list
            if ind > 0:
                
                # calculate new values of y' from the input function, which
                # is a list of functions obtained from the subs z = dy/dt
                y_prime_new = y_prime[ind-1] + h*func(time[ind-1],f[ind-1],y_prime[ind-1])[1]
                
                # append the new y' value to the empty list
                y_prime.append(y_prime_new)
                
                # find the new y values using 
                y_new = f[ind-1] + h*func(time[ind-1],f[ind-1],y_prime[ind])[0]
                
                # append new y value
                f.append(y_new)
                
        return f,y_prime, time

def heun(func, initial, time, step = None):
    
    '''
    Takes a list of one or two funciton and a list of one or two intial conditions
    
    Input:
        func: a function of two variables used for the calculation, obtained by setting 
            z = dydt (func[0]) and its derivative dzdt = d2ydt2 (func[1])
        initial: the initial conditions y = [y(0),dydt(0)]
        time: either the full series, or stop and start times if step is defined
        step: the amount by which to step the calculation between t[0] and t[1]
    Output
        list of values pertaining to the analytical solution for the second
        order differential equation input
        '''
    
    if not step == None:
        import numpy as np
        h = step
        n = (time[-1] - time[0]/h + 1)
        time = np.linspace(time[0], time[-1],n)
        
    else:
        h = (time[-1] - time[0])/(len(time)-1)
        
    f = []
    y_prime = []
    
    if len(initial) == 1:
        f.append(initial[0])
        
        for ind, t in enumerate(time):
            
            if ind > 0:
                
                # Using the Rung-Katta 2 method
                # The first calculation: predictor
                f_z = func(time[ind-1], f[ind-1])[0]
                
                # secondary calculation: corrector
                f_z2 = func(t,f[ind-1] + f_z*h)[0]
                
                # putting calculation together
                z = f[ind -1] + 0.5*h*(f_z + f_z2)
                
                # appending to empty lists
                f.append(z)
        return f, time
                
    else:
        f.append(initial[0])
        y_prime.append(initial[1])
        
        for ind, t in enumerate(time):
            
            if ind > 0:
                
                # Using the Rung-Katta 2 method
                # The first calculation: predictor
                f_dzdt = func(time[ind -1], f[ind-1], y_prime[ind-1])[1]
                f_z = func(time[ind-1], f[ind-1],  y_prime[ind-1])[0]
                
                # secondary calculation: corrector
                f_dzdt2 = func(t, f[ind-1] + f_z*h, y_prime[ind-1]+f_dzdt*h)[1]
                f_z2 = func(t,f[ind-1] + f_z*h,y_prime[ind-1]+f_dzdt*h)[0]
                
                # putting calculation together
                dzdt = y_prime[ind-1] + 0.5*h*(f_dzdt + f_dzdt2)
                z = f[ind -1] + 0.5*h*(f_z + f_z2)
                
                # appending to empty lists
                y_prime.append(dzdt)
                f.append(z)
            
            
        return f,y_prime, time

def rk4(func, initial, time, step = None):
    from astropy import units as u
    '''
    Takes a list of one or two funciton and a list of one or two intial conditions
    
    Input:
        func: a function of two variables used for the calculation, obtained by setting 
            z = dydt (func[0]) and its derivative dzdt = d2ydt2 (func[1])
        initial: the initial conditions y = [y(0),dydt(0)]
        time: either the full series, or stop and start times if step is defined
        step: the amount by which to step the calculation between t[0] and t[1]
    Output
        list of values pertaining to the analytical solution for the second
        order differential equation input
        '''
        
    if not step == None:
        import numpy as np
        h = step
        n = ((time[-1] - time[0])/h + 1)
        time = np.linspace(time[0], time[-1],n)
        
    else:
        h = (time[-1] - time[0])/(len(time)-1)
    print(h)
    f = []
    y_prime = []
        
    if len(initial) == 2:
        f.append(initial[0])
        y_prime.append(initial[1])
        
        for i, t in enumerate(time):
            
            
            if i > 0:
                
                # Rung-Katta 4 method
                
#                if time[i-1] == 0:
#                    k1_dzdt = initial[1]
##                    print('k1 =',k1_dzdt)
#                    k1_z = initial[0]
##                    print('k1z =',k1_dzdt)
                    
                    
#                else:
                # k1 for each variable
                k1_dzdt = func(time[i -1], f[i-1], y_prime[i-1])[1]
#                    print('k1 =',k1_dzdt)
                k1_z = func(time[i-1], f[i-1],  y_prime[i-1])[0]
#                    print('k1z =',k1_dzdt)
                
                
                # k2 for each variable
                k2_dzdt = func(time[i -1] + h/2, f[i-1] + k1_z*h/2, y_prime[i-1]+k1_dzdt*h/2)[1]
#                print('k2 =',k2_dzdt)
                k2_z = func(time[i -1] + h/2,f[i-1] + k1_z*h/2,y_prime[i-1]+k1_dzdt*h/2)[0]
#                print('k2z =',k2_z)
                
                # k3 for each variable
                k3_dzdt = func(time[i-1] + h/2, f[i-1] + k2_z*h/2, y_prime[i-1] + k2_dzdt*h/2)[1]
#                print('k3 =',k3_dzdt)
                k3_z = func(time[i-1] + h/2, f[i-1] + k2_z*h/2, y_prime[i-1] + k2_dzdt*h/2)[0]
#                print('k3z =',k3_z)
                
                # k4 for each variable
                k4_dzdt = func(t, f[i-1] + k3_z*h, y_prime[i-1] + k3_dzdt*h)[1]
#                print('k4 =',k4_dzdt)
                k4_z = func(t, f[i-1] + k3_z*h, y_prime[i-1] + k3_dzdt*h)[0]
#                print('k4z =',k4_z)
                
                # putting calculations together
                dzdt = y_prime[i-1] + (1/6)*h*(k1_dzdt + 2*k2_dzdt + 2*k3_dzdt + k4_dzdt)
#                print('dzdt =',dzdt)
                z = f[i-1] + (1/6)*h*(k1_z + 2*k2_z + 2*k3_z + k4_z)
#                print('z =',z,'\n')
                
                # appending to lists
                y_prime.append(dzdt)
                f.append(z)
                
        return f, y_prime, time
    
    else:
        f.append(initial[0])
        
        for i, t in enumerate(time):
            
            if i > 0:
                
                # Rung-Katta 4 method, multiple corrector steps
                # k1 for each variable
                k1_z = func(time[i-1], f[i-1])[0]
                
                # k2 for each variable
                k2_z = func(time[i -1] + h/2,f[i-1] + k1_z*h/2)[0]
                
                # k3 for each variable
                k3_z = func(time[i-1] + h/2, f[i-1] + k2_z*h/2)[0]
                
                # k4 for each variable
                k4_z = func(t, f[i-1] + k3_z*h)[0]
                
                # putting calculations together
                z = f[i-1] + (1/6)*h*(k1_z + 2*k2_z + 2*k3_z + k4_z)
                
                # appending to lists
                f.append(z)
                
        return f, time