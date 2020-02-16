# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 18:51:28 2020

@author: Ryan
"""

import numpy as np
import matplotlib.pyplot as plt
import diffeq as df
from scipy.integrate import odeint

def examplepend(y,t,b,c):
    theta,omega = y
    dydt = [omega,-b*omega - c*np.sin(theta)]
    return dydt


def mypend(t,y,y1,b=0.25,c=5.0):
    theta = y 
    omega = y1
    dydt = [omega,-b*omega - c*np.sin(theta)]
    return dydt

b = 0.25
c = 5.0

y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0,10,101)



sol1 = odeint(examplepend, y0, t, args = (b,c))

# Testing Euler method
theta1, omega1, time = df.euler(mypend, y0,t, step = 0.01)    

error_theta1 = sol1[:,0] - theta1 # comparing two methods by finding difference
error_omega1 = sol1[:,1] - omega1

plt.figure()
plt.plot(time,theta1, label = 'Theta Difference')
plt.plot(time,omega1, label = 'Omega Difference')
plt.plot(t, error_theta1, label = 'Theta Difference')
plt.plot(t,error_omega1, label = 'Omega Difference')
plt.legend()

# Testing Heun's Method
theta2, omega2, time = df.heun(mypend, y0, t)

error_theta2 = sol1[:,0] - theta2 # comparing two methods by finding difference
error_omega2 = sol1[:,1] - omega2

plt.figure()
plt.plot(t, theta2, label = 'Theta')
plt.plot(t, omega2, label = 'Omega')
plt.plot(t, error_theta2, label = 'eTheta')
plt.plot(t, error_omega2, label = 'eOmega')
plt.legend()
