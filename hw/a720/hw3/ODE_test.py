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
    theta, omega = y 
    dydt = [omega,-b*omega - c*np.sin(theta)]
    return dydt

def pend(t,y,y1,b=0.25,c=5.0):
    theta = y 
    omega = y1
    dydt = [omega,-b*omega - c*np.sin(theta)]
    return dydt

b = 0.25
c = 5.0

y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0,10,101)



sol1 = odeint(examplepend, y0, t, args = (b,c))

# Testing Euler method, plotting
theta1, omega1, time = df.euler(pend, y0,t)    
r1_theta = np.corrcoef(sol1[:,0], theta1)[0,1]
r1_omega = np.corrcoef(sol1[:,1], omega1)[0,1]

plt.figure()
plt.title('Euler Method')
plt.plot(time,theta1, label = 'Theta')
plt.plot([],[], label = f"R^2 = {format(r1_theta,'.5')}")
plt.plot(time,omega1,'--', label = 'Omega')
plt.plot([],[], label = f"R^2 = {format(r1_omega,'.5')}")
plt.legend()

# Testing Heun's Method, plotting
theta2, omega2, time2 = df.heun(pend, y0, t)
r2_theta = np.corrcoef(sol1[:,0], theta2)[0,1]
r2_omega = np.corrcoef(sol1[:,1], omega2)[0,1]

plt.figure()
plt.title('Huen Method (rk2)')
plt.plot(t, theta2, label = 'Theta')
plt.plot([],[], label = f"R^2 = {format(r2_theta,'.5')}")
plt.plot(t, omega2,'--', label = 'Omega')
plt.plot([],[], label = f"R^2 = {format(r2_omega,'.5')}")
plt.legend()

# Testing RK4 method, plotting
theta3, omega3, time3 = df.rk4(pend, y0, t)
r3_theta = np.corrcoef(sol1[:,0], theta3)[0,1]
r3_omega = np.corrcoef(sol1[:,1], omega3)[0,1]

plt.figure()
plt.title('RK4 method')
plt.plot(t, theta3, label = 'Theta')
plt.plot([],[], label = f"R^2 = {format(r3_theta,'.5')}")
plt.plot(t, omega3, '--', label = 'Omega')
plt.plot([],[], label = f"R^2 = {format(r3_omega,'.5')}")
plt.legend()


# Fun little plot
L = 0.5 # pendulum length in meters
x = L*np.sin(theta3)
y = L - L*np.cos(theta3)

plt.figure()
plt.scatter(x,y)
plt.scatter(0,L, color = 'r')

# Homework Question 3
def stiff(t,y,lam = 12):
    dydt = [-lam*(y - np.cos(t))]
    return dydt

y1 = [0]

# get solutions from each method
sol4_euler, time_euler = df.euler(stiff, y1, t)
sol4_huen, time_huen = df.heun(stiff, y1, t)
sol4_rk4, time_rk4 = df.rk4(stiff, y1, t)

# create data points for the provided solution to the stif ode
y_sol = []
lam = 12
for i in t:
    a = lam/(1+lam**2)*np.exp(-lam*i)
    b = lam/(1+lam**2)*np.sin(i)
    c = lam**2/(1+lam**2)*np.cos(i)
    f = a + b + c
    
    y_sol.append(f)

# plot huen results
plt.figure()
plt.plot(t, y_sol, label = 'provided solution')
plt.plot(t,sol4_huen, label = 'huen')
plt.legend()

# plot euler results
plt.figure()
plt.plot(t,sol4_euler, label = 'euler')
plt.plot(t, y_sol, label = 'provided solution')
plt.legend()

# plot rk4 results
plt.figure()
plt.plot(t,sol4_rk4, label = 'rk4')
plt.plot(t, y_sol, label = 'provided solution')
plt.legend()


