# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 20:13:52 2020

@author: Ryan
"""

import cluster_dynamics as cd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

file = 'galaxies0.npy'
file2 = 'galaxies1.npy'

gal1 = np.load(file)
gal2 = np.load(file2)

X = (0,11)
Y = (0,11)
Z = (0,11)

myTree = cd.Tree(X, Y, Z, cd.createGalList(gal1,gal2))
print(myTree.root.rCoM)
galaxies = myTree.galArray()

print(myTree.root.children)

gal = myTree.root.children[0].children[0].children[0].children[0].children[0].gals[0]
myTree.calcAccel(myTree.root, gal)
print(gal.accel)


difference = gal2[0,0] - gal1[0,0]
galList= myTree.galaxies


gals2 = myTree.calcPos(1000)
def calca(galaxy, galList):
    accel = [0,0,0]
    G = 4.5172e-48 # Mpc^3/ s^2 / solar mass
    
    for i in galList:
        if i == galaxy:
            continue
        r = float(np.sqrt((i.x - galaxy.x)**2 + (i.y - galaxy.y)**2 + (i.z - galaxy.z)**2))
        if r == 0:
            continue
        M = i.M
        accel[0] += float(G * M * (i.x - galaxy.x)/r**3)
        accel[1] += float(G * M * (i.y - galaxy.y)/r**3)
        accel[2] += float(G * M * (i.z - galaxy.z)/r**3)
    print(accel)
        
    return accel


calca(gal, galList)
a = np.array([[1,2,3]])

a = np.append(a, [[1,2,3]], axis = 0)

fig = plt.figure()
ax = plt.Axes3D(fig)

for i in range(len(gal1)):
    ax.scatter(gal1[i,0], gal1[i,1], gal1[i,2], zdir = 'z')
