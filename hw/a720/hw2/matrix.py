# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 14:25:56 2020

@author: Ryan
"""
import numpy as np


class matrix:
    
    def __init__(self, m, n, filename = 'None'):
        '''
        Input:
            m: the number of rows
            n: the number of columns
            filename: the file from which to generate the matrix
        Output:
            Either a blank matrix full of zeroes with m rows and n columns
            or a matrix generated from a text file using lists
        '''
        self.m = m
        self.n = n
        
        self.values = [[0]*n for x in range(m)] # create matrix full of zeroes
        
        if not filename == 'None':
             
            valList = np.loadtxt(filename, delimiter = ',')
            
            # appending values to matrix based on given text file
            for row in valList:
                self.values[int(row[0])-1][int(row[1])-1] = float(row[2])
                
        
            
        
    def __add__(self, matrixB):
        '''
        Input:
            self: first matrix to add
            matrixB: second matrix to add
            
        Output:
            The matrix result of an addition between self and matrixB
        '''
        
        # cheking that the matrices have the same sizes
        assert self.m == matrixB.m
        assert self.n == matrixB.n
        
        answer = matrix(self.m,self.n)  # creating blank answer matrix
        
        # iterate through rows, then columns
        for r, row in enumerate(self.values):
            for c, val in enumerate(row):
                
                # add the values there together
                answer.values[r][c] = val + matrixB.values[r][c]
                
        return answer
    
    def __sub__(self, matrixB):
        '''
        Input:
            self: first matrix
            matrixB: matrix to subtract from self
            
        Output: 
            resulting matrix
        '''
        
        # again, check they are the same sizes
        assert self.m == matrixB.m
        assert self.n == matrixB.n
        
        answer = matrix(self.m,self.n)
        
        for r, row in enumerate(self.values):
            for c, val in enumerate(row):
                answer.values[r][c] = val - matrixB.values[r][c]
                
        return answer
    
        
    def transpose(self):
        '''
        Input:
            self: the matrix to find the transpose of
        
        Output:
            The transpose of the given matrix
        '''
        
        # generate answer matrix
        answer = matrix(self.n,self.m)
        
        # iterate through rows, then columns
        for r, row in enumerate(self.values):
            for c, val in enumerate(row):
                
                # reverse order of the column and row values
                answer.values[c][r] = val 
                
        return answer
    
    def __mul__(self,matrixB):
        '''
        Input:
            self: the first matrix in the multiplication
            matrixB: the second matrix in the multiplication
        Output:
            the product of the two matrices A*B = C
        '''
        
        # assert that the columns of matrix 1 have the same size as the
        # rows of matrix 2
        assert self.n == matrixB.m
        
        # generate blank matrix
        answer = matrix(self.m,matrixB.n)
        
        # iterate through rows, then columns rows of transpose of matrix 2
        for r, row in enumerate(self.values):
            for c, col in enumerate(matrixB.transpose().values):
                
                total = 0
                
                # iterate through values and multiply and add them together
                for i in range(len(row)):
                    x = row[i]
                    y = col[i]
                    total = total + y*x
                    
                # append the calculate value to its appropriate row and col
                answer.values[r][c] = total
                
        return answer       
                

