#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:33:16 2020

@author: ryan
"""

from matrix import matrix
import numpy as np
import matplotlib as plt

file = '/home/ryan/hw/a720/hw2/A_coefficients.dat'
# Testing add    
matrixa = matrix(9,9,filename = file)
matrixb = matrix(9,9, filename = file)

ansa = matrixa+matrixb

#print(ansa.values)

# Testing Transpose
ansb = matrixa.transpose()
#print(ansb.values)

# Testing multiplication
matrixc = matrix(1,3,filename = 'test1.dat')
matrixd = matrix(3,2, filename = 'test2.dat')

mult = matrixa*matrixb
#print(mult.values)
#print(matrixc.value, matrixd.values)



