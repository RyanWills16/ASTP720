# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 14:25:56 2020

@author: Ryan
"""
import numpy as np


class matrix:
    
    def __init__(self, m, n, filename = 'None'):
        self.m = m
        self.n = n
        
        self.values = [[0]*n for x in range(m)]
        
        if not filename == 'None':
             
            valList = np.loadtxt(filename, delimiter = ',')
            
            for row in valList:
                self.values[int(row[0])-1][int(row[1])-1] = float(row[2])
                
        
            
        
    def __add__(self, matrixB):
        
        assert self.m == matrixB.m
        assert self.n == matrixB.n
        
        answer = matrix(self.m,self.n)
        
        for r, row in enumerate(self.values):
            for c, val in enumerate(row):
                answer.values[r][c] = val + matrixB.values[r][c]
                
        return answer
    
    def __sub__(self, matrixB):
        
        assert self.m == matrixB.m
        assert self.n == matrixB.n
        
        answer = matrix(self.m,self.n)
        
        for r, row in enumerate(self.values):
            for c, val in enumerate(row):
                answer.values[r][c] = val - matrixB.values[r][c]
                
        return answer
    
        
    def transpose(self):
        answer = matrix(self.n,self.m)
        
        for r, row in enumerate(self.values):
            for c, val in enumerate(row):
                
                answer.values[c][r] = val 
                
        return answer
    
    def __mul__(self,matrixB):
        assert self.n == matrixB.m
        answer = matrix(self.m,matrixB.n)
        
        for r, row in enumerate(self.values):
            for c, col in enumerate(matrixB.transpose().values):
                
                total = 0
                
                for i in range(len(row)):
                    x = row[i]
                    y = col[i]
                    total = total + y*x
                answer.values[r][c] = total
                
        return answer       
                
#    def trace(self):
#        for r, row in enumerate(self.values):
#            for c, col in enumerate(matrixB.transpose().values):

