# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 15:55:13 2020

@author: Ryan
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u


class galaxy:
    
    def __init__(self, x, y, z, M = 1e12):
        
        '''
        Initialize galaxy object with its coordinates and position vector
        
        '''
        
        self.x = x
        self.y = y
        self.z = z
        self.M = M
        self.r = (self.x, self.y, self.z)

class node:
    
    def __init__(self, x_b, y_b, z_b, parent=None, numGal = None):
        
        '''
        Initialize node with is min and max values in each dimension, children,
        galaxies, p
        
        '''
        
        self.xmin, self.xmax = x_b
        self.ymin, self.ymax = y_b
        self.zmin, self.zmax = z_b
        
        if parent is not None:
            self.parent = parent
        
        self.children = [] # list of child node objects
        
        self.gals = [] # list of galaxy objects
        
        self.parent_gal = 0 # number of galaxies under the node if a parent
        
        self.COM = 0
        
        if numGal is not None:
            self.numGal = numGal # same as parent_gal, unless root or leaf
        
    def addchild(self, child):
        
        '''
        Method to add child nodes to parent node
        
        '''
        self.children.append(child)
        return
    
    def addgal(self, gal):
        
        '''
        Method to add galaxy objects to node
        
        '''
        
        self.gals.append(gal)
        return
        
class Tree:
    
    def __init__(self, X_mb, Y_mb, Z_mb, file):
        
        '''
        Initialize the tree with the root node and the list of galaxy positions
        
        '''
        
        self.galList = np.load(file)
        
        self.root = node(X_mb, Y_mb, Z_mb)
        
        
    def createChild(self, parent):
        
        '''
        
        Method for creating child nodes from input parents based on galaxy positions
        
        '''
        xmid = (parent.xmin + parent.xmax)/2
        ymid = (parent.ymin + parent.ymax)/2
        zmid = (parent.zmin + parent.zmax)/2
        
        x0 = parent.xmin
        x1 = parent.xmax
        
        y0 = parent.ymin
        y1 = parent.ymax
        
        z0 = parent.zmin
        z1 = parent.zmax
        
        g = self.galList
        
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        n5 = 0
        n6 = 0
        n7 = 0
        n8 = 0 
        
        # iterate through the galaxy list
        # if there's 1 or more galaxies in that region, add 1 to n#
        # node 1
        for i in range(len(g)):
            if g[i,0] >= x0 and g[i,0] < xmid and g[i,1] >= y0 and g[i,1] < ymid and g[i,2] >= z0 and g[i,2] < zmid:
                n1 += 1
                
        # if n# is greater than zero, create a child
        if n1 > 0:
            parent.addchild(node((x0,xmid),(y0,ymid),(z0,zmid), parent = parent, numGal = n1))
            parent.parent_gal += n1
            
        # node 2
        for i in range(len(g)):
            if g[i,0] >= x0 and g[i,0] < xmid and g[i,1] >= y0 and g[i,1] < ymid and g[i,2] >= zmid and g[i,2] < z1:
                n2 += 1
        if n2 > 0:
            parent.addchild(node((x0,xmid),(y0,ymid),(zmid,z1), parent = parent, numGal = n2))
            parent.parent_gal += n2
        
        # node 3
        for i in range(len(g)):
            if g[i,0] >= x0 and g[i,0] < xmid and g[i,1] >= ymid and g[i,1] < y1 and g[i,2] >= zmid and g[i,2] < z1:
                n3 += 1
        if n3 > 0:
            parent.addchild(node((x0,xmid),(ymid,y1),(zmid,z1), parent = parent, numGal = n3))
            parent.parent_gal += n3
        
        # node 4
        for i in range(len(g)):
            if g[i,0] >= x0 and g[i,0] < xmid and g[i,1] >= ymid and g[i,1] < y1 and g[i,2] >= z0 and g[i,2] < zmid:
                n4 += 1
        if n4 > 0:
            parent.addchild(node((x0,xmid),(ymid,y1),(z0,zmid), parent = parent, numGal = n4))
            parent.parent_gal += n4
            
        # node 5
        for i in range(len(g)):
            if g[i,0] >= xmid and g[i,0] < x1 and g[i,1] >= y0 and g[i,1] < ymid and g[i,2] >= z0 and g[i,2] < zmid:
                n5 += 1
        if n5 > 0:
            parent.addchild(node((xmid,x1),(y0,ymid),(z0,zmid), parent = parent, numGal = n5))
            parent.parent_gal += n5
        
        # node 6
        for i in range(len(g)):
            if g[i,0] >= xmid and g[i,0] < x1 and g[i,1] >= y0 and g[i,1] < ymid and g[i,2] >= zmid and g[i,2] < z1:
                n6 += 1
        if n6 > 0:
            parent.addchild(node((xmid,x1),(y0,ymid),(zmid,z1), parent = parent, numGal = n6))
            parent.parent_gal += n6
        
        # node 7
        for i in range(len(g)):
            if g[i,0] >= xmid and g[i,0] < x1 and g[i,1] >= ymid and g[i,1] < y1 and g[i,2] >= zmid and g[i,2] < z1:
                n7 += 1
        if n7 > 0:
            parent.addchild(node((xmid,x1),(ymid,y1),(zmid,z1), parent = parent, numGal = n7))
            parent.parent_gal += n7
        
        # node 8
        for i in range(len(g)):
            if g[i,0] >= xmid and g[i,0] < x1 and g[i,1] >= ymid and g[i,1] < y1 and g[i,2] >= z0 and g[i,2] < zmid:
                n8 += 1
        if n8 > 0:
            parent.addchild(node((xmid,x1),(ymid,y1),(z0,zmid), parent = parent, numGal = n8))
            parent.parent_gal += n8
        
    def constructTree(self,parent):
        
        '''
        Method to construct the entire tree and add galaxies to nodes
        
        Works recursively
        
        '''
        # for the input node, which should be the root node initially, create
        # children
        self.createChild(parent)
        
        # iterate through the children of the parent node
        for node in parent.children:
            
            # if there are multiple galaxies in the node, call the method again
            if node.numGal > 1:
                self.constructTree(node)
                
            # if there is one galaxy, create galaxy object and add to node
            else:
                x0 = node.xmin
                x1 = node.xmax
                y0 = node.ymin
                y1 = node.ymax
                z0 = node.zmin
                z1 = node.zmax
                
                g = self.galList
                
                for i in range(len(g)):
                    if g[i,0] >= x0 and g[i,0] < x1 and g[i,1] >= y0 and g[i,1] < y1 and g[i,2] >= z0 and g[i,2] < z1:
                        node.addgal(galaxy(g[i,0],g[i,1], g[i,2]))
                    
                