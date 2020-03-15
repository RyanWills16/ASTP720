# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 20:13:52 2020

@author: Ryan
"""

import cluster_dynamics as cd

file = 'galaxies0.npy'

X = (0,11)
Y = (0,11)
Z = (0,11)

myTree = cd.Tree(X, Y, Z, file)

print(myTree.root.children)

test = myTree.createChild(myTree.root)
print(myTree.root.children)

print(myTree.root.children[6].xmax)

