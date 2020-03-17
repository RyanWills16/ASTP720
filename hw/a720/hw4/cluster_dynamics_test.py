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
print(myTree.root.parent_gal)

print(myTree.root.children[1].numGal)


test2 = myTree.createChild(myTree.root.children[1])

print(myTree.root.children[1].children[2].numGal)
print(myTree.root.children[1].parent)
print(myTree.root)


test3 = myTree.createChild(myTree.root.children[1].children[2])

print(len(myTree.root.children[1].children[2].children))

print(myTree.root.children[1].children[2].children[1].numGal)

myTree2 = cd.Tree(X,Y,Z,file)

myTree2.constructTree(myTree2.root)

print(len(myTree2.root.children))

print(myTree2.root.children[0].children[0].children[0].children[0].children[0].gals[0].M)

print(myTree2.root.children[0].children[0].children[0].children[0].children[0].xmax)

myTree3 = cd.Tree(X, Y, Z, file)
myTree3.constructTree(myTree3.root)
myTree3.root.calcCOM()

print(myTree3.root.children[0].children[0].numGal)
