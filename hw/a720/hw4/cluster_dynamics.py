# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 15:55:13 2020

@author: Ryan
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u


galaxies = np.load('galaxies0.npy')
class galaxy:
    
    def __init__(self, x, y, z, M = 1e12):
        
        self.x = x
        self.y = y
        self.z = z
        self.M = M

class node:
    
    def __init__(self, x_b, y_b, z_b, parent=None):
        
        self.xmin, self.xmax = x_b
        self.ymin, self.ymax = y_b
        self.zmin, self.zmax = z_b
        
        if parent is not None:
            self.parent = parent
        
        self.children = []
        
        self.galaxies = []
        
        self.numGal = 0
        
    def addchild(self, child):
        self.children.append(child)
        return
    
    def addgal(self, gal):
        self.galaxies.append(gal)
        return
        
class Tree:
    
    def __init__(self, X_mb, Y_mb, Z_mb, file):
        
        self.galList = np.load(file)
        
        self.root = node(X_mb, Y_mb, Z_mb)
        
        
    def createChild(self, parent):
        xmid = (parent.xmin + parent.xmax)/2
        ymid = (parent.ymin + parent.ymax)/2
        zmid = (parent.zmin + parent.zmax)/2
        
        x0 = parent.xmin
        x1 = parent.xmax
        
        y0 = parent.ymin
        y1 = parent.ymax
        
        z0 = parent.zmin
        z1 = parent.zmax
        
        n1 = 0
        n2 = 0
        n3 = 0
        n4 = 0
        n5 = 0
        n6 = 0
        n7 = 0
        n8 = 0 
        
        # node 1
        for i in range(len(self.galList)):
            if self.galList[i,0] >= x0 and self.galList[i,0] < xmid and self.galList[i,1] >= y0 and self.galList[i,1] < ymid and self.galList[i,2] >= z0 and self.galList[i,2] < zmid:
                n1 += 1
        if n1 > 0:
            parent.addchild(node((x0,xmid),(y0,ymid),(z0,zmid), parent = parent))
            parent.numGal = n1
        
        # node 2
        for i in range(len(self.galList)):
            if self.galList[i,0] >= x0 and self.galList[i,0] < xmid and self.galList[i,1] >= y0 and self.galList[i,1] < ymid and self.galList[i,2] >= zmid and self.galList[i,2] < z1:
                n2 += 1
        if n2 > 0:
            parent.addchild(node((x0,xmid),(y0,ymid),(zmid,z1), parent = parent))
            parent.numGal = n2
        
        # node 3
        for i in range(len(self.galList)):
            if self.galList[i,0] >= x0 and self.galList[i,0] < xmid and self.galList[i,1] >= ymid and self.galList[i,1] < y1 and self.galList[i,2] >= zmid and self.galList[i,2] < z1:
                n3 += 1
        if n3 > 0:
            parent.addchild(node((x0,xmid),(ymid,y1),(zmid,z1), parent = parent))
            parent.numGal = n3
        
        # node 4
        for i in range(len(self.galList)):
            if self.galList[i,0] >= x0 and self.galList[i,0] < xmid and self.galList[i,1] >= ymid and self.galList[i,1] < y1 and self.galList[i,2] >= z0 and self.galList[i,2] < zmid:
                n4 += 1
        if n4 > 0:
            parent.addchild(node((x0,xmid),(ymid,y1),(z0,zmid), parent = parent))
            parent.numGal = n4
            
        # node 5
        for i in range(len(self.galList)):
            if self.galList[i,0] >= xmid and self.galList[i,0] < x1 and self.galList[i,1] >= y0 and self.galList[i,1] < ymid and self.galList[i,2] >= z0 and self.galList[i,2] < zmid:
                n5 += 1
        if n5 > 0:
            parent.addchild(node((xmid,x1),(y0,ymid),(z0,zmid), parent = parent))
            parent.numGal = n5
        
        # node 6
        for i in range(len(self.galList)):
            if self.galList[i,0] >= xmid and self.galList[i,0] < x1 and self.galList[i,1] >= y0 and self.galList[i,1] < ymid and self.galList[i,2] >= zmid and self.galList[i,2] < z1:
                n6 += 1
        if n6 > 0:
            parent.addchild(node((xmid,x1),(y0,ymid),(zmid,z1), parent = parent))
            parent.numGal = n6
        
        # node 7
        for i in range(len(self.galList)):
            if self.galList[i,0] >= xmid and self.galList[i,0] < x1 and self.galList[i,1] >= ymid and self.galList[i,1] < y1 and self.galList[i,2] >= zmid and self.galList[i,2] < z1:
                n7 += 1
        if n7 > 0:
            parent.addchild(node((xmid,x1),(ymid,y1),(zmid,z1), parent = parent))
            parent.numGal = n7
        
        # node 8
        for i in range(len(self.galList)):
            if self.galList[i,0] >= xmid and self.galList[i,0] < x1 and self.galList[i,1] >= ymid and self.galList[i,1] < y1 and self.galList[i,2] >= z0 and self.galList[i,2] < zmid:
                n8 += 1
        if n8 > 0:
            parent.addchild(node((xmid,x1),(ymid,y1),(z0,zmid), parent = parent))
            parent.numGal = n8
        
        