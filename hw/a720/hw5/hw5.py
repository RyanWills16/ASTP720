# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 22:04:14 2020

@author: Ryan
"""

import numpy as np
import matplotlib.pyplot as plt

# loading in relevant columns from the data file
data = np.genfromtxt('cepheid_data.txt', delimiter = ',', skip_header = 1, usecols = (1,2,3,7,8))

# calculating absolute magnitude using apparent mag, d in pc, Rv = 3.1, and E(B-V) for each obs
M = []
for i in range(len(data)):
    M.append(data[i,2] + 5 - 5*np.log10(data[i,1]*10**3) - 3.1*data[i,3])

# defining data vector for P and data vecotr for metallicity
P = np.log10(data[:,0])
z = data[:,4]
def cmodel(x, z, a0, a1, a2):
    return a0 + a1*x + a2*z

X = np.zeros([len(M), 3])
for i in range(len(M)):
    X[i,0] = 1
    X[i,1] = P[i]
    X[i,2] = data[i,4]
    
# matrix math for parameter estimates
par = (np.linalg.inv((np.transpose(X) @ X))) @ (np.transpose(X) @ M)

# matrix math for uncertainty in parameter estimates
unc = np.sqrt(np.linalg.inv((np.transpose(X) @ X)))

# calculating points for fit line
x_mod = np.array([np.min(P),np.max(P)])
z_mod = np.array([z[301],z[205]])
y_mod = cmodel(x_mod, z_mod, par[0], par[1], par[2])


plt.figure(figsize = (5.5,4))
plt.tight_layout()
plt.scatter(P,M, s = 4, c = 'b')
plt.plot(x_mod,y_mod,c='r',label = 'Best Fit:')
# plt.errorbar(x,y_n,yerr = unc, fmt = 'o', markersize = 3, lw = 0.75)
plt.plot([],[],' ', label = rf"$\alpha$ = {format(par[0],'0.3f')} $\pm$ {format(unc[0,0],'0.3f')}")
plt.plot([],[], ' ', label = rf"$\beta$ = {format(par[1],'0.3f')} $\pm$ {format(unc[1,1],'0.3f')}")
plt.plot([],[], ' ', label = rf"$\gamma$ = {format(par[2],'0.3f')} $\pm$ {format(unc[2,2],'0.3f')}")
plt.ylabel('Absolute Magnitude')
plt.xlabel('$log_{10}(Period) (log_{10}(days))$')
plt.title("Without Weighting")
plt.tight_layout()
plt.legend()

#uncertainty
u = 0.1 # magnitudes
# noise matrix
N = np.zeros([len(P),len(P)])
# fill the noise matrix with 1/uncertainty^2, which is the inverse
np.fill_diagonal(N,1/u**2)
# matrix math to solve for parameters with weighting
par2 = (np.linalg.inv((np.transpose(X) @ N) @ X)) @ (((np.transpose(X) @ N) @ M))

unc2 = np.sqrt(np.linalg.inv((np.transpose(X)@ N @ X)))

# calculating points for fit line
x_mod2 = np.array([np.min(P),np.max(P)])
z_mod2 = np.array([z[301],z[205]])
y_mod2 = cmodel(x_mod2, z_mod2, par2[0], par2[1], par2[2])

plt.figure(figsize = (5.5,4))
plt.tight_layout()
plt.scatter(P,M, s = 4, c = 'b')
plt.plot(x_mod2,y_mod2,c='r',label = 'Best Fit:')
# plt.errorbar(x,y_n,yerr = unc, fmt = 'o', markersize = 3, lw = 0.75)
plt.plot([],[],' ', label = rf"$\alpha$ = {format(par2[0],'0.3f')} $\pm$ {format(unc2[0,0],'0.3f')}")
plt.plot([],[], ' ', label = rf"$\beta$ = {format(par2[1],'0.3f')} $\pm$ {format(unc2[1,1],'0.3f')}")
plt.plot([],[], ' ', label = rf"$\gamma$ = {format(par2[2],'0.3f')} $\pm$ {format(unc2[2,2],'0.3f')}")
plt.ylabel('Absolute Magnitude')
plt.xlabel('$log_{10}(Period) (log_{10}(days))$')
plt.title("With Weighting")
plt.tight_layout()
plt.legend()