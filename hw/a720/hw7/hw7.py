# -*- coding: utf-8 -*-
"""
Created on Fri May  1 17:16:16 2020

@author: Ryan
"""
import numpy as np
import matplotlib.pyplot as plt

data = np.load('strain.npy')
time = np.linspace(0,1048575,1048576)
time_seconds = time * 60
time_day = time * 0.000694444

plt.figure()
plt.plot(time_day, data, lw = 0.75)
plt.xlabel('Time (Days)')
plt.ylabel('Strain h')

fnyq = 1/(2*(time_seconds[1] - time_seconds[0]))
freq = np.fft.rfftfreq(len(data), time_seconds[1] - time_seconds[0])/fnyq

fourier = np.fft.rfft(data-np.average(data))
fourier = 2*fourier/len(data)
fourier[0] = fourier[0]/2
fourier[-1] = fourier[-1]/2
fourier = abs(fourier * np.conj(fourier))

plt.figure()
plt.plot(np.log10(freq), np.log10(fourier))
plt.xlabel('$log_{10}(f/f_{nyq})$')
plt.ylabel('log_{10} (Power)')

f = 10**(-1.85732)*fnyq
amp = np.sqrt(10**(-40.06))

D = 12
c1 = 2.6e-21
c2 = 10e-4

M = ((D**3*amp**3*c2**2)/(c1**3*f*2))**(1/5)

print(M, 'solar masses')
