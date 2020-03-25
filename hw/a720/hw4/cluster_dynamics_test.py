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
    
    gal1 = np.load(file1)
    gal2 = np.load(file2)
    tree = cd.Tree(X, Y, Z, cd.createGalList(gal1,gal2))
    
    print(f"This will progress the system by {format(num_run*step,'.3e')} years")
    
    galaxy1 = np.array([[tree.galaxies[0].x, tree.galaxies[0].y, tree.galaxies[0].z]])
    galaxy2 = np.array([[tree.galaxies[-1].x, tree.galaxies[-1].y, tree.galaxies[-1].z]])
    galaxy3 = np.array([[tree.galaxies[225].x, tree.galaxies[225].y, tree.galaxies[225].z]])
    galaxy4 = np.array([[tree.galaxies[675].x, tree.galaxies[675].y, tree.galaxies[675].z]])
    
    
    n = 0
    while n < num_run:
#        print(tree)
        n += 1
        print(f"{num_run - n}")
        new_gals = tree.calcPos(step)
        galaxy1 = np.append(galaxy1, [[new_gals[0].x, new_gals[0].y, new_gals[0].z]], axis = 0)
        galaxy2 = np.append(galaxy2, [[new_gals[-1].x, new_gals[-1].y, new_gals[-1].z]], axis = 0)
        galaxy3 = np.append(galaxy3, [[new_gals[225].x, new_gals[225].y, new_gals[225].z]], axis = 0)
        galaxy4 = np.append(galaxy4, [[new_gals[675].x, new_gals[675].y, new_gals[675].z]], axis = 0)
        
        tree = cd.Tree((tree.xmin, tree.xmax),(tree.ymin, tree.ymax),(tree.zmin, tree.zmax), new_gals)
        
    end = time.time()
    print((end - start)/3600, 'hours')
    
    return galaxy1, galaxy2, galaxy3, galaxy4
        
    
num_run = 10000
step = 1000000
X = (0,11)
Y = (0,11)
Z = (0,11)
file1 = 'galaxies0.npy'
file2 = 'galaxies1.npy'

g1, g2, g3, g4 = nbody(num_run, step, X, Y, Z, file1, file2)

# write each of the galaxies to a text file for later reference
# galaxy 1

step_years = '_1mil'
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
    
g1 = np.loadtxt('galaxy1' + step_years + '.dat')
g2 = np.loadtxt('galaxy2' + step_years + '.dat')   
g3 = np.loadtxt('galaxy3' + step_years + '.dat')   
g4 = np.loadtxt('galaxy4' + step_years + '.dat')   
fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
#ax.set_xlim(0,13)
#ax.set_zlim(0,13)
#ax.set_ylim(0,13)
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
