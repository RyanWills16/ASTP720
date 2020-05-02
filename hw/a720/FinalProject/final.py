# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 00:48:28 2020

@author: Ryan
"""

class particle:
    
    def __init__(self, x, y, vx, vy):
        
        # position in kpc
        self.x = x
        self.y = y
        
        # velocity in kpc/s
        self.velx = vx
        self.vely = vy
        
        # position history
        self.x_locations = [self.x]
        self.y_locations = [self.y]
        
        
        
class galaxy:
    
    def __init__(self, x, y, vx, vy, Mass, name = ''):
        
        # position in kpc
        self.x = x
        self.y = y
        
        # velocity in kpc/s
        self.velx = vx*(1/3.086e16)
        self.vely = vy*(1/3.086e16)
        
        # Mass in solar Masses
        self.Mass = Mass
        
        self.name = name
        
        self.particles = []
        
        # position history
        self.x_locations = [self.x]
        self.y_locations = [self.y]
    
    def createParticle(self, R, num, rot = 'counter'):
        import numpy as np
        
        G = 4.513e-39 # (kpc^3) / (s^2 solar mass)
        theta = (2*np.pi)/num
        
        for i in range(0, num+1):
            angle = i*theta
            
            x = R*np.cos(angle) + self.x
            y = R*np.sin(angle) + self.y
            
            if rot == 'counter':
                vx = -1*np.sqrt(G * self.Mass / R)*np.sin(angle) + self.velx
                vy = np.sqrt(G * self.Mass / R)*np.cos(angle) + self.vely
                
            if rot == 'clockwise':
                vx = np.sqrt(G * self.Mass / R)*np.sin(angle) + self.velx
                vy = -1*np.sqrt(G * self.Mass / R)*np.cos(angle) + self.vely
            
            self.addParticle(particle(x,y,vx,vy))
    
            
    def addParticle(self, particle):
        
        self.particles.append(particle)
        return

def calc_acceleration(galaxy_list, body, time_step):
    import numpy as np
    
    # convert time step to seconds
    time_step = time_step * 3.154e7
    
    velx = 0
    vely = 0
    
    posx = 0
    posy = 0
    
    for i in galaxy_list:
        
        # do not calculate forces on yourself
        if i == body:
            continue
        
        else:
            G = 4.513e-39 # (kpc^3) / (s^2 solar mass)
            r = np.sqrt((i.x - body.x)**2 + (i.y - body.y)**2)
            
            C = G * i.Mass / r**3
            
            # k1
            k1x = C * (i.x - body.x)
            k1y = C * (i.y - body.y)
            
            # k2
            vel2x = body.velx + k1x*0.5*time_step
            vel2y = body.vely + k1y*0.5*time_step
            
            p2x = body.x + vel2x*0.5*time_step
            p2y = body.y + vel2y*0.5*time_step
            
            k2x = (i.x - p2x)*C
            k2y = (i.y - p2y)*C
            
            # k3
            vel3x = body.velx + k2x*0.5*time_step
            vel3y = body.vely + k2y*0.5*time_step
            
            p3x = body.x + vel3x*0.5*time_step
            p3y = body.y + vel3y*0.5*time_step
            
            k3x = (i.x - p3x)*C
            k3y = (i.y - p3y)*C

            # k4
            vel4x = body.velx + k3x*time_step
            vel4y = body.vely + k3y*time_step
            
            p4x = body.x + vel4x*time_step
            p4y = body.y + vel4y*time_step
            
            k4x = (i.x - p4x)*C
            k4y =  (i.y - p4y)*C
            
            # calculate new velocity
            velx += time_step * (k1x + k2x + k3x + k4x) / 6
            vely += time_step * (k1y + k2y + k3y + k4y) / 6
            
            # calculate new position
            posx += time_step * (body.velx + vel2x + vel3x + vel4x) / 6
            posy += time_step * (body.vely + vel2y + vel3y + vel4y) / 6
            
            # update particle location history
            body.x_locations.append(body.x + posx)
            body.y_locations.append(body.y + posy)
            
            # update particle position
            body.x += posx
            body.y += posy
            
            # update particle velocity
            body.velx += velx
            body.vely += vely

def progress_system(galaxy_list,  time_step, num_step):
    
    print(f"This will progress the system by {format(time_step*num_step, '.3e')} years")
    
    n = 0
    
    # make a list of all bodies in simulation
    bodies = []
    
    for i in galaxy_list:
        bodies.append(i)
        
        if len(i.particles) == 0:
            continue
        
        else:
            for j in i.particles:
                bodies.append(j)
                
    while n < num_step:
        n += 1
        
        for i in bodies:
            calc_acceleration(galaxy_list, i, time_step)
            
def plotGal(galaxy):
    plt.figure()
    plt.scatter(t.x, t.y, s = 20, c='black', zorder = 2)
    plt.axvline(x = 0, c = 'gray', zorder = 0)
    plt.axhline(y = 0, c = 'gray', zorder = 0)
    for i in galaxy.particles:
        plt.scatter(i.x_locations, i.y_locations, s = 10, c = 'r')
        plt.scatter(i.x_locations[0], i.y_locations[0], s = 10, c = 'black')
        
def plotInteraction(galaxy_list, index_list):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    fig, ax = plt.subplots(6,2, sharex = True, sharey = True, wspace = 0, hspace = 0)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which = 'both')
        
            
import matplotlib.pyplot as plt

t = galaxy(1,1, 0, 0, 1e11, name = 'G1')

t.createParticle(2, 4)

progress_system([t], 1000, 4000)

plotGal(t)




    
