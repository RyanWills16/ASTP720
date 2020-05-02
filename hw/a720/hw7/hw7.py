# -*- coding: utf-8 -*-
"""
Created on Fri May  1 17:16:16 2020

@author: Ryan
"""
import numpy as np
import matplotlib.pyplot as plt

data = np.load('strain.npy')

fourier = np.fft.rfft(data)
fabs = 2*abs(fourier * np.conj(fourier))**2 / len(data)

time = np.linspace(0,1048575,1048576)
time_seconds = time * 60
time_day = time * 0.000694444

fnyq = 1/(2*(time_seconds[1] - time_seconds[0]))
freq = np.fft.rfftfreq(len(data), time_seconds[1] - time_seconds[0])/fnyq



plt.plot(time_day, data, lw = 0.75)
plt.xlabel('Time (Days)')
plt.ylabel('Strain h')

plt.plot(np.log10(freq), np.log10(fabs))
