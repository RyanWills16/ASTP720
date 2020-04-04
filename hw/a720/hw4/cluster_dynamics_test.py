# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 20:13:52 2020

@author: Ryan
"""

import cluster_dynamics as cd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# loading in files with galaxy position info
file = 'galaxies0.npy'
file2 = 'galaxies1.npy'
gal1 = np.load(file)
gal2 = np.load(file2)

# defining bounds of the box
X = (0,11)
Y = (0,11)
Z = (0,11)

# creating tree from list of galaxy objects, checking that children are 
# generated and that the galaxies are loaded
myTree = cd.Tree(X, Y, Z, cd.createGalList(gal1,gal2))
print(len(myTree.galaxies))
print(myTree.root.children)

# testing calcPos method, checking that the first galaxy in each list is
# the same galaxy, but with a slight changed position
print(myTree.galaxies[0].x, myTree.galaxies[0].y, myTree.galaxies[0].z)
gals2 = myTree.calcPos(1000)
print(gals2[0].x, gals2[0].y, gals2[0].z)

# checking change in x position at different times
print(gal2[0,0] - gal1[0,0])
print(gals2[0].x - myTree.galaxies[0].x)

# testing my Barnes-Hut by using direct summation
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


galaxy = myTree.galaxies[0]
galList= myTree.galaxies

# direct summation acceleration
sum_accel = calca(galaxy, galList)

# my Barnes-Hut method acceleration
myTree.calcAccel(myTree.root, galaxy)

# this is in units of Mpc/s^2
print(sum_accel, myTree.galaxies[0].accel)

# testing 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

ax.set_xlim(0,12)
ax.set_ylim(0,12)
ax.set_zlim(0,12)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.plot(gal1[:,0], gal1[:,1], gal1[:,2], 'o', c = 'gray', ms = 1)
plt.show()
        

# running simulation
num_run = 100
step = 1e7
X = (0,11)
Y = (0,11)
Z = (0,11)
file1 = 'galaxies0.npy'
file2 = 'galaxies1.npy'

g1, g2, g3, g4 = cd.nbody(num_run, step, X, Y, Z, file1, file2)

# write each of the galaxies to a text file for later reference
# galaxy 1

step_years = str(format(step, '.0e'))
with open('galaxy1' + step + '.dat', 'w') as f:
    for i in range(len(g1)):
        f.write(f"{format(g1[i,0], '15.15f')} \t {format(g1[i,1], '15.15f')} \t {format(g1[i,2], '15.15f')}\n")
    f.close()

# galaxy 2
with open('galaxy2' + step + '.dat', 'w') as f:
    for i in range(len(g1)):
        f.write(f"{format(g2[i,0], '15.15f')} \t {format(g2[i,1], '15.15f')} \t {format(g2[i,2], '15.15f')}\n")
    f.close()
    
# galaxy 3 
with open('galaxy3' + step + '.dat', 'w') as f:
    for i in range(len(g1)):
        f.write(f"{format(g3[i,0], '15.15f')} \t {format(g3[i,1], '15.15f')} \t {format(g3[i,2], '15.15f')}\n")
    f.close()

# galaxy 4
with open('galaxy4' + step + '.dat', 'w') as f:
    for i in range(len(g1)):
        f.write(f"{format(g4[i,0], '15.15f')} \t {format(g4[i,1], '15.15f')} \t {format(g4[i,2], '15.15f')}\n")
    f.close()

# load the files of galaxy trajectories and plot in 3D
g1 = np.loadtxt('galaxy1' + step_years + '.dat')
g2 = np.loadtxt('galaxy2' + step_years + '.dat')   
g3 = np.loadtxt('galaxy3' + step_years + '.dat')   
g4 = np.loadtxt('galaxy4' + step_years + '.dat')   
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.set_xlim(0,12)
ax.set_zlim(0,12)
ax.set_ylim(0,12)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.tight_layout()

ax.plot([], [], [], label = "10 billion years")
ax.plot(g1[:,0], g1[:,1], g1[:,2], 'o', ms = 2, color = 'r')
ax.plot(g2[:,0], g2[:,1], g2[:,2], 'o', ms = 2, color = 'b')
ax.plot(g3[:,0], g3[:,1], g3[:,2], 'o', ms = 2, color = 'g')
ax.plot(g4[:,0], g4[:,1], g4[:,2], 'o', ms = 2, color = 'black')
plt.legend()

    
plt.show()


# plot nodes bounds.
print(myTree.numNode)

def nodeList(node, nodes):
    nodes.append(node)
    if len(node.children) > 0:
        for child in node.children:
            nodes.append(child)
            nodeList(child, nodes)

nodes = []
nodeList = nodeList(myTree.root,nodes)

plt.figure(figsize = (5.5,5.5))
plt.plot(gal2[:,0], gal2[:,1], '^', color = 'r', markersize = 1)
plt.axis('equal')
plt.ylabel('Y')
plt.xlabel('X')
for node in nodes:
    plt.plot((node.xmin, node.xmax), (node.ymin, node.ymin), 'black', linewidth = 0.5)
    plt.plot((node.xmin, node.xmax), (node.ymax, node.ymax), 'black', linewidth = 0.5)
    plt.plot((node.xmin, node.xmin), (node.ymin, node.ymax), 'black', linewidth = 0.5)
    plt.plot((node.xmax, node.xmax), (node.ymin, node.ymax), 'black', linewidth = 0.5)
    
    
