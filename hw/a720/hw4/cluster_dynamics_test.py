# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 20:13:52 2020

@author: Ryan
"""

import cluster_dynamics as cd
import numpy as np
file = 'galaxies0.npy'
file2 = 'galaxies1.npy'

gal1 = np.load(file)
gal2 = np.load(file2)

X = (0,11)
Y = (0,11)
Z = (0,11)



myTree = cd.Tree(X, Y, Z, gal1)

print(len(myTree.galaxies))

gal = myTree.root.children[0].children[0].children[0].children[0].children[0].gals[0]
myTree.calcAccel(myTree.root, gal)



