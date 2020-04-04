# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 15:55:13 2020

@author: Ryan
"""

class galaxy:
    
    def __init__(self, x, y, z, x_p = 0, y_p = 0, z_p = 0, M = 1e12):
        
        '''
        Initialize galaxy object with its coordinates and position vector and 
        acceleration vector
        '''
        
        # galaxy's current coordinates for a particular iteration
        self.x = x
        self.y = y
        self.z = z
        self.M = M
        self.r = [x, y, z]
        
        # galaxy's previous coordinates for a particular iteration, used in 
        # the verlet algorithm
        self.x_p = x_p
        self.y_p = y_p
        self.z_p = z_p
        self.prev_r = [x_p, y_p, z_p]
        
        self.accel = [0,0,0]



class node:
    
    def __init__(self, x_b, y_b, z_b, parent=None, numGal = None):
        
        '''
        Initialize node with is min and max values in each dimension, children,
        galaxies, parent, center of mass, and total mass
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
        Calculate the center of mass of each node, takes the tree root node
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
    
    def __init__(self, X_mb, Y_mb, Z_mb, galaxies):
        
        '''
        Initialize the tree with the root node and the list of galaxy positions.
        Construct the tree and calculate the CoM of all nodes
        '''
        
        self.xmin, self.xmax = X_mb
        self.ymin, self.ymax = Y_mb
        self.zmin, self.zmax = Z_mb
        
        # list of galaxy objects, used in 
        self.galaxies = galaxies
        
        self.numNode = 0
        
        self.root = node(X_mb, Y_mb, Z_mb)
        
        # array of galaxy positions, used in createChild
        self.galList = self.galArray()
        
        # construct the Tree
        self.constructTree(self.root)
        
        # calc all node COMs
        self.root.calcCOM()
    
    def galArray(self):
        import numpy as np
        
        a = np.array([[self.galaxies[0].x,self.galaxies[0].y,self.galaxies[0].z]])
        
        for ind, i in enumerate(self.galaxies):
            if ind > 0:
                a = np.append(a, [[i.x,i.y,i.z]], axis = 0)
                
        return a
                
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
                
                for i in self.galaxies:
                    if i.x >= x0 and i.x < x1 and i.y >= y0 and i.y < y1 and i.z >= z0 and i.z < z1:
                        node.addgal(i)
                        
                        node.CoM_num[0] = node.gals[0].x*node.gals[0].M
                        node.CoM_num[1] = node.gals[0].y*node.gals[0].M
                        node.CoM_num[2] = node.gals[0].z*node.gals[0].M
                        
                        # make the center of mass of leaves the galaxy position
                        node.rCoM[0] = node.gals[0].x
                        node.rCoM[1] = node.gals[0].y
                        node.rCoM[2] = node.gals[0].z
                        
                        # make total mass of leaf the mass of the galaxy
                        node.totM = node.gals[0].M
        
    def calcAccel(self, node, galaxy):
        import numpy as np
        
        '''
        method calculates the acceleration on a given galaxy due to all other 
        galaxies or nodes.
        '''
        
        G = 4.5172e-48 # Mpc^3/ s^2 / solar mass
        
        ga = galaxy.accel
        
        # if the node has children
        if len(node.children) > 0:
            for child in node.children:
                
                # calculate distance between galaxy and node
                r = np.sqrt((child.rCoM[0] - galaxy.x)**2 + (child.rCoM[1] - galaxy.y)**2 + (child.rCoM[2] - galaxy.z)**2)
                node_size = child.xmax - child.xmin
                
                # if the two references are the same galaxy, skip the node
                if r == 0:
                    continue
                
                # if the node's size is larger than the distance to the node,
                # call the method recursively for this node
                if node_size/r > 1:
                    self.calcAccel(child, galaxy)
                
                # if the node's size is smaller than r, calculate accel using
                # the center of mass of the node and its total mass enclosed
                else:
                    # force softening condition
                    if r < 0.01:
#                        print('softened')
                        M = child.totM
                        e = 0.005 # Mpc
                        ga[0] += G * M * (child.rCoM[0] - galaxy.x)/((r**2 + e**2)*r)
                        ga[1] += G * M * (child.rCoM[1] - galaxy.y)/((r**2 + e**2)*r)
                        ga[2] += G * M * (child.rCoM[2] - galaxy.z)/((r**2 + e**2)*r)
                        
                    else:
#                        print('calc')
                        M = child.totM
                        ga[0] += G * M * (child.rCoM[0] - galaxy.x)/r**3
                        ga[1] += G * M * (child.rCoM[1] - galaxy.y)/r**3
                        ga[2] += G * M * (child.rCoM[2] - galaxy.z)/r**3
        
        # if the node does not have children           
        elif len(node.children) == 0:
            r = np.sqrt((node.rCoM[0] - galaxy.x)**2 + (node.rCoM[1] - galaxy.y)**2 + (node.rCoM[2] - galaxy.z)**2)
            
            # if the two references are the same galaxy, add nothing to accel
            if r == 0:
                ga[0] += 0
                ga[1] += 0
                ga[2] += 0
            
            # if the node is a leaf and therefore has only one galaxy, add to accel
            else: 
                if r < 0.01:
#                    print('softened2')
                    M = node.totM
                    e = 0.005 # Mpc
                    ga[0] += G * M * (node.rCoM[0] - galaxy.x)/((r**2 + e**2)*r)
                    ga[1] += G * M * (node.rCoM[1] - galaxy.y)/((r**2 + e**2)*r)
                    ga[2] += G * M * (node.rCoM[2] - galaxy.z)/((r**2 + e**2)*r)
                    
                else:
#                    print('calc2')
                    M = node.totM
                    ga[0] += G * M * (node.rCoM[0] - galaxy.x)/r**3
                    ga[1] += G * M * (node.rCoM[1] - galaxy.y)/r**3
                    ga[2] += G * M * (node.rCoM[2] - galaxy.z)/r**3
                    
            
            
    def calcPos(self, step):
        
        '''
        method for calculating the new position of the galaxy by calling 
        calcAccel. The step size is the time step to take in years.
        '''
        
        galaxies = self.galaxies
        new_galaxies = []
        step = step * 3.154e7 # convert step in years to step in seconds
        
        for i in galaxies:
            self.calcAccel(self.root, i)
            
            new_x = 2*i.x - i.x_p + (step**2)*(i.accel[0])
            new_y = 2*i.y - i.y_p + (step**2)*(i.accel[1])
            new_z = 2*i.z - i.z_p + (step**2)*(i.accel[2])
            
            new_galaxies.append(galaxy(new_x, new_y, new_z, x_p = i.x, y_p = i.y, z_p = i.z))
        
        return new_galaxies
            
def createGalList(gal_array1, gal_array2):
    
    '''
    Method used to create list of galaxy objects from the two .npy files
    used to initialize the tree. 
    
    gal_array1 is from galaxies0.npy, and gal_array2 is from galaxies1.npy.
    '''
    
    g1 = gal_array1
    g2 = gal_array2
    galaxies = []
    for i in range(len(g1)):
        galaxies.append(galaxy(g2[i,0],g2[i,1],g2[i,2], x_p = g1[i,0], y_p = g1[i,1], z_p = g1[i,2]))
        
    return galaxies                     

# nbody simulation over time
def nbody(num_run, step, X, Y, Z, file1, file2):
    
    '''
    num_run = number of times to run the simulation
    step = step size of the simulation in years
    '''
    import cluster_dynamics as cd
    import numpy as np
    import time
    start = time.time()
    
    # load galaxies in from files
    gal1 = np.load(file1)
    gal2 = np.load(file2)
    
    # create tree
    tree = cd.Tree(X, Y, Z, cd.createGalList(gal1,gal2))
    
    print(f"This will progress the system by {format(num_run*step,'.3e')} years")
    
    # create arrays of galaxy 3D coordinates
    galaxy1 = np.array([[tree.galaxies[0].x, tree.galaxies[0].y, tree.galaxies[0].z]])
    galaxy2 = np.array([[tree.galaxies[-1].x, tree.galaxies[-1].y, tree.galaxies[-1].z]])
    galaxy3 = np.array([[tree.galaxies[225].x, tree.galaxies[225].y, tree.galaxies[225].z]])
    galaxy4 = np.array([[tree.galaxies[675].x, tree.galaxies[675].y, tree.galaxies[675].z]])
    
    
    n = 0
    while n < num_run:
#        print(tree)
        n += 1
        # calculat new galaxy positions
        new_gals = tree.calcPos(step)
        
        # append new positions to arrays
        galaxy1 = np.append(galaxy1, [[new_gals[0].x, new_gals[0].y, new_gals[0].z]], axis = 0)
        galaxy2 = np.append(galaxy2, [[new_gals[-1].x, new_gals[-1].y, new_gals[-1].z]], axis = 0)
        galaxy3 = np.append(galaxy3, [[new_gals[225].x, new_gals[225].y, new_gals[225].z]], axis = 0)
        galaxy4 = np.append(galaxy4, [[new_gals[675].x, new_gals[675].y, new_gals[675].z]], axis = 0)
        
        # create new tree
        tree = cd.Tree((tree.xmin, tree.xmax),(tree.ymin, tree.ymax),(tree.zmin, tree.zmax), new_gals)
        
    end = time.time()
    print((end - start)/3600, 'hours')
    
    return galaxy1, galaxy2, galaxy3, galaxy4