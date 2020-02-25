# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 18:51:28 2020

@author: Ryan
"""

import numpy as np
import matplotlib.pyplot as plt
import diffeq as df
from scipy.integrate import odeint
from astropy import units as u

# defining functions for question 1
# the ode_int example
def examplepend(y,t,b,c):
    theta, omega = y 
    dydt = [omega,-b*omega - c*np.sin(theta)]
    return dydt

# modified it to my specifications
def pend(t,y,y1,b=0.25,c=5.0):
    theta = y 
    omega = y1
    
    # output list of functions
    dydt = [omega,-b*omega - c*np.sin(theta)]
    return dydt

# define coefficients
b = 0.25
c = 5.0

# define initial conditions, time series
y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0,10,101)

# generate ode_int solution
sol1 = odeint(examplepend, y0, t, args = (b,c))

# Testing my Euler method, plotting
theta1, omega1, time = df.euler(pend, y0,t)

    # calculating correlation coefficient 
r1_theta = np.corrcoef(sol1[:,0], theta1)[0,1]
r1_omega = np.corrcoef(sol1[:,1], omega1)[0,1]

    # plotting
plt.figure(figsize = (6.75,5))
plt.xlabel('Time', fontsize = 14)
plt.title('Euler Method')
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(time,theta1, label = 'Theta')
plt.plot([],[], label = f"R^2 = {format(r1_theta,'.5')}")
plt.plot(time,omega1,'--', label = 'Omega')
plt.plot([],[], label = f"R^2 = {format(r1_omega,'.5')}")
plt.tight_layout()
plt.legend(fontsize = 12)
plt.savefig('euler1.pdf')

# Testing my Heun's Method, plotting
theta2, omega2, time2 = df.heun(pend, y0, t)
r2_theta = np.corrcoef(sol1[:,0], theta2)[0,1]
r2_omega = np.corrcoef(sol1[:,1], omega2)[0,1]

plt.figure(figsize = (6.75,5))
plt.xlabel('Time', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.title('Huen Method (rk2)')
plt.plot(t, theta2, label = 'Theta')
plt.plot([],[], label = f"R^2 = {format(r2_theta,'.5')}")
plt.plot(t, omega2,'--', label = 'Omega')
plt.plot([],[], label = f"R^2 = {format(r2_omega,'.5')}")
plt.tight_layout()
plt.legend(fontsize = 12)
plt.savefig('huen1.pdf')

# Testing RK4 method, plotting
theta3, omega3, time3 = df.rk4(pend, y0, t)
r3_theta = np.corrcoef(sol1[:,0], theta3)[0,1]
r3_omega = np.corrcoef(sol1[:,1], omega3)[0,1]

plt.figure(figsize = (6.75,5))
plt.xlabel('Time', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.title('RK4 method')
plt.plot(t, theta3, label = 'Theta')
plt.plot([],[], label = f"R^2 = {format(r3_theta,'.5')}")
plt.plot(t, omega3, '--', label = 'Omega')
plt.plot([],[], label = f"R^2 = {format(r3_omega,'.5')}")
plt.tight_layout()
plt.legend(fontsize = 12)
plt.savefig('rk41.pdf')

# Fun little plot of x,y coordinates of pendulum
# not very useful
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

y1 = [1]

# get solutions from each method
sol4_euler, time_euler = df.euler(stiff, y1, t)
sol4_huen, time_huen = df.heun(stiff, y1, t)
sol4_rk4, time_rk4 = df.rk4(stiff, y1, t)

# create data points for the provided solution to the stif ode
y_sol = []
lam = 25
for i in t:
    a = lam/(1+lam**2)*np.exp(-lam*i)
    b = lam/(1+lam**2)*np.sin(i)
    c = lam**2/(1+lam**2)*np.cos(i)
    f = a + b + c
    
    y_sol.append(f)

# plot huen results
plt.figure(figsize = (6.75,5))
plt.xlabel('Radius (R_sun)', fontsize = 14)
plt.ylabel('Menc(r) (g)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.scatter(t, y_sol, label = 'solution', color = 'black')
plt.plot(t,sol4_huen, label = 'huen', color = 'green')
plt.tight_layout()
plt.legend(fontsize = 12)
plt.savefig('huen2.pdf')

# plot euler results
plt.figure(figsize = (6.75,5))
plt.xlabel('Radius (R_sun)', fontsize = 14)
plt.ylabel('Menc(r) (g)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.scatter(t,sol4_euler, label = 'euler', color = 'black')
plt.plot(t, y_sol, label = 'solution',color = 'green')
plt.tight_layout()
plt.legend(fontsize = 12)
plt.savefig('euler2.pdf')

# plot rk4 results
plt.figure(figsize = (6.75,5))
plt.xlabel('Radius (R_sun)', fontsize = 14)
plt.ylabel('Menc(r) (g)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.scatter(t,sol4_rk4, label = 'rk4', color = 'black')
plt.plot(t, y_sol, label = 'solution', color = 'green')
plt.tight_layout()
plt.legend(fontsize = 12)
plt.savefig('rk42.pdf')

# White Dwarf Hydrostatic equilibrium
def hydroequil(r, M, P, mu=2):
    G = 6.667*10**-8
    r_sun = 69.634*10**9    # radius of sun in cm
    m_sun =  1.989*10**33 # mass of sun in grams
    
    # defining some constants for the pressure and mass equations
    C_mass = 4*np.pi*mu*(1e13)**(-3/5)*r_sun**2
    C_pres = G*mu*(1e13)**(-3/5 )*r_sun**-2
    
    # return list of functions
    function = [C_mass*r**2*P**(3/5), -C_pres*M*P**(3/5)*r**-2]
    return function

# define functions for ODE_int
def hydroequil2(init,r, mu=2):
    G = 6.667*10**-8
    r_sun = 69.634*10**9    # radius of sun in cm
    m_sun =  1.989*10**33 # mass of sun in grams
    
    m,p = init # unpacking initial conditions
    
    # defining some constants for the pressure and mass equations
    C_mass = 4*np.pi*mu*(1e13)**(-3/5)*r_sun**2
    C_pres = G*mu*(1e13)**(-3/5 )*r_sun**-2
    
    # return list of functions
    function = [C_mass*r**2*P**(3/5), -C_pres*M*P**(3/5)*r**-2]
    return function

# define radius in terms of stellar radii
radius = np.linspace(1e-10, 3, 1000)

# initial conditions
M = 0
P = 1*10**13*((10**4)/2)**(5/3)
initial = [M,P]

# get rk4 solution
mass, pressure, radius = df.rk4(hydroequil, initial, radius)

# plot radius and mass
plt.figure(figsize = (6.75,5))
plt.xlabel('Radius (R_sun)', fontsize = 14)
plt.ylabel('Menc(r) (g)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(radius, mass, color = 'black')
plt.tight_layout()
plt.savefig('massenc.pdf')

# plot radius and pressure
plt.figure(figsize = (6.75,5))
plt.xlabel('Radius (R_sun)', fontsize = 14)
plt.ylabel('Menc(r) (g)', fontsize = 14)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.plot(radius, pressure, color = 'black')
plt.tight_layout()
plt.savefig('massenc.pdf')

# ode_int solution
# returns the same thing as my solution
hyrdo = odeint(hydroequil2, initial, radius)

