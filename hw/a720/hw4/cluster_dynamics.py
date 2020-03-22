# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 15:55:13 2020

@author: Ryan
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u


class galaxy:
    
    def __init__(self, x, y, z, x_p = 0, y_p = 0, z_p = 0, M = 1e12):
        
        '''
        Initialize galaxy object with its coordinates and position vector
        
        '''
        
        self.x = x
        self.y = y
        self.z = z
        self.M = M
        self.r = [x, y, z]
        
        self.accel = [0,0,0]
        
        prev_r = [x_p, y_p, z_p]

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
        
        self.CoM_num = [0,0,0] # numerator of the center of mass equation
        
        self.totM = 0 # total mass within a node
        
        self.rCoM = [0,0,0] # CoM radius vector [x,y,z]
        
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
    
    def calcCOM(self):
        '''
        Calculate the center of mass of each node
        
        '''
        if len(self.children) == 0:
            self.CoM_num[0] = self.gals[0].x
            self.CoM_num[1] = self.gals[0].y
            self.CoM_num[2] = self.gals[0].z
            
            return self.gals[0].x, self.gals[0].y, self.gals[0].z, self.gals[0].M
        
        else:
            
            for child in self.children:
                
                com = child.calcCOM()
                
                self.CoM_num[0] += com[0]*com[3]
                self.CoM_num[1] += com[1]*com[3]
                self.CoM_num[2] += com[2]*com[3]
                self.totM += com[3]
            
            self.rCoM[0] = self.CoM_num[0] / self.totM
            self.rCoM[1] = self.CoM_num[1] / self.totM
            self.rCoM[2] = self.CoM_num[2] / self.totM
            
            return self.rCoM[0], self.rCoM[1], self.rCoM[2], self.totM
            
        
            
        
        
        
class Tree:
    
    def __init__(self, X_mb, Y_mb, Z_mb, galList):
        
        '''
        Initialize the tree with the root node and the list of galaxy positions
        
        '''
        
        # array of galaxy positions
        self.galList = galList
        
        self.numNode = 0
        
        self.root = node(X_mb, Y_mb, Z_mb)
        
        # list of galaxy objects
        self.galaxies = []
        
        self.constructTree(self.root)
        
        self.getGals(self.root)
        
        self.root.calcCOM()
        
        
        
    def createChild(self, parent):
        
        '''
        
        Helper Method for constructTree method which actually creats child nodes 
        from input parents based on galaxy positions
        
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
            self.numNode += 1
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
                        
                        node.CoM_num[0] = node.gals[0].x*node.gals[0].M
                        node.CoM_num[1] = node.gals[0].y*node.gals[0].M
                        node.CoM_num[2] = node.gals[0].z*node.gals[0].M
                        
                        # make the center of mass of leaves the galaxy position
                        node.rCoM[0] = node.gals[0].x
                        node.rCoM[1] = node.gals[0].y
                        node.rCoM[2] = node.gals[0].z
                        
                        # make total mass of leaf the mass of the galaxy
                        node.totM = node.gals[0].M
    
    def getGals(self, node):
        if len(node.children) == 0:
            self.galaxies.append(node.gals[0])
            
        else: 
            for child in node.children:
                self.getGals(child)
        
    def calcAccel(self, node, galaxy):
        
        G = 4.5172e-48 # Mpc^3/ s^2 / solar mass
        
        
        accel = [0,0,0]
                
        if len(node.children) > 0:
            for child in node.children:
                
                r = np.sqrt((child.rCoM[0] - galaxy.x)**2 + (child.rCoM[1] - galaxy.y)**2 + (child.rCoM[2] - galaxy.z)**2)
                node_size = child.xmax - child.xmin
                
                if r == 0:
                    continue
                
                if node_size/r > 1:
                    self.calcAccel(child, galaxy)
                    
                else:
                    M = child.totM
                    accel[0] += G * M * (child.rCoM[0] - galaxy.x)/r**3
                    accel[1] += G * M * (child.rCoM[1] - galaxy.y)/r**3
                    accel[2] += G * M * (child.rCoM[2] - galaxy.z)/r**3
                    
        else:
            r = np.sqrt((node.rCoM[0] - galaxy.x)**2 + (node.rCoM[1] - galaxy.y)**2 + (node.rCoM[2] - galaxy.z)**2)
            if r == 0:
                accel[0] += 0
                accel[1] += 0
                accel[2] += 0
                
            else: 
                accel[0] += G * M * (node.rCoM[0] - galaxy.x)/r**3
                accel[1] += G * M * (node.rCoM[1] - galaxy.y)/r**3
                accel[2] += G * M * (node.rCoM[2] - galaxy.z)/r**3
                
        return accel
                    
            
            
    def calcPos(self, step):
        galaxies = self.galaxies
        new_galaxies = []
        step = step * 3.154e7 # convert step in years to step in seconds
        
        for i in galaxies:
            accel = self.calcAccel(self.root, i, step )
            
            new_x = 2*i.x - i.x_p + (step**2)*(accel[0])
            new_y = 2*i.y - i.y_p + (step**2)*(accel[1])
            new_z = 2*i.z - i.z_p + (step**2)*(accel[2])
            
            new_galaxies.append(galaxy(new_x, new_y, new_z, x_p = i.x, y_p = i.y, z_p = i.z))
        
        return new_galaxies
            
                    
        
#    def calcAccel(self, node, x, y, z, node2):
#        
#        '''
#        Calculate acceleration on node2 due to node
#        '''
#        
#        G = 4.5172e-48
#        
#        if node == node2:
#            return
#        
#        r = np.sqrt((node.rCoM[0] - x)**2 + (node.rCoM[1] - y)**2 + (node.rCoM[3] - x)**2)
#        node_size = node.xmax - node.xmin
#        
#        if node_size/r > 1:
#            for child in node.children:
#                self.calcAccel(child, x, y, z, node2)
#        else:
#            M = node.totM
#            x_a = G * M * (node.rCoM[0] - x)/r**3
#            y_a = G * M * (node.rCoM[0] - y)/r**3
#            z_a = G * M * (node.rCoM[0] - x)/r**3
#            
#        return x_a, y_a, z_a
#    
#    def calcPos(self, nd):
#        
#        '''
#        Calculate new position using help of function calcAccel
#        '''
#        
#        galaxies = []
#        accel = [0,0,0]
#        if nd.children > 0:
#            
#            for node in nd.children:
#                self.calcPos(node)
#                
#        else:
#            x = nd.gals[0].x
#            y = nd.gals[0].y
#            z = nd.gals[0].z
#            
#            a = self.calcAccel(self.root, x, y, z, nd)
#            
#            if a is not None:
#                accel[0] += a[0]
#                accel[1] += a[1]
#                accel[2] += a[2]
            