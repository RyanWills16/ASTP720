# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 00:48:28 2020

@author: Ryan
"""

class particle:
    
    '''
    
    particle object
    
    '''
    
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
    
    '''
    
    Galaxy object
    
    '''
    
    def __init__(self, x, y, vx, vy, Mass, name = ''):
        
        # position in kpc
        self.x = x
        self.y = y
        
        # velocity in kpc/s
        self.velx = vx*(1/3.086e16) # convert km/s to kpc/s
        self.vely = vy*(1/3.086e16)
        
        # Mass in solar Masses
        self.Mass = Mass
        
        # name of galaxy
        self.name = name
        
        # particles in galaxy
        self.particles = []
        
        # position history
        self.x_locations = [self.x]
        self.y_locations = [self.y]
    
    def createParticle(self, R, num, rot = 'counter'):
        
        '''
        
        Method for creating particles in a galaxy
        
        R:
            Radius in kpc
            
        num:
            number of particles to generate
            
        rot: which way the particles should orbit
        
        '''
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
    
    '''
    
    updates positions and velocities for all objects using RK4
    
    Takes the list of galaxies, and a single body to make calculations on
    
    '''
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
            
            # but still update my new position
            body.x += body.velx * time_step
            body.y += body.vely * time_step
            
            # update particle history
            body.x_locations.append(body.x)
            body.y_locations.append(body.y)
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
    
    '''
    Builds list of all particles and galaxies and then evolves them forwards 
    in time. Also builds time array.
    '''
    
    print(f"This will progress the system by {format(time_step*num_step, '.3e')} years")
    
    n = 0
    t = 0
    time = [0]
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
        t += time_step
        time.append(t)
        
        for i in bodies:
            calc_acceleration(galaxy_list, i, time_step)
            
    return time
            
def plotGal(galaxy):
    
    '''
    plots all particle histories for a single galaxy, recommended only for 
    use in showing particle orbits, not galaxy movement. 
    
    '''
    
    plt.figure()
    plt.scatter(galaxy.x_locations, galaxy.y_locations, s = 20, c='black', zorder = 2)
    plt.gca().axis('equal')
    plt.axvline(x = 0, c = 'gray', zorder = 0)
    plt.axhline(y = 0, c = 'gray', zorder = 0)
    plt.xlabel('X (kpc)')
    plt.ylabel('Y (kpc)')
    for i in galaxy.particles:
        plt.scatter(i.x_locations, i.y_locations, s = 10, c = 'r')
        plt.scatter(i.x_locations[0], i.y_locations[0], s = 10, c = 'black')
        
def plotInteraction(galaxy_list, index_list, time):
    
    '''
    galaxy_list: 
        List of galaxies to plot along with their particles
    
    index_list: 
        The indicies of the times you want to plot, must be a list of 
    six indicies corresponding to inidicies in time
    
    time: time array as put out by progress_system
    
    Output:
        plot the galaxy and particle positions as a function of time in 
        6 subplots
    
    must be 6 instances in time to plot
    

    '''
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    fig, ax = plt.subplots(2,3,figsize = (12, 4.5), sharex = True, sharey = True, gridspec_kw = {'wspace':0,'hspace':0})
    ax[0,0].set_ylabel('Y (kpc)')
    ax[1,0].set_ylabel('Y (kpc)')
    ax[1,0].set_xlabel('X (kpc)')
    ax[1,1].set_xlabel('X (kpc)')
    ax[1,2].set_xlabel('X (kpc)')
    plt.tight_layout()
    
    figs = [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2)]
    
    for ind, i in enumerate(index_list):
        ax[figs[ind][0],figs[ind][1]].plot([], [], ' ', label = f"{format(time[i],'.3e')}")
        for k in galaxy_list:
            ax[figs[ind][0],figs[ind][1]].scatter(k.x_locations[i], k.y_locations[i], s = 24, c = 'black')
            
            ax[figs[ind][0],figs[ind][1]].xaxis.set_minor_locator(AutoMinorLocator())
            ax[figs[ind][0],figs[ind][1]].xaxis.set_minor_locator(AutoMinorLocator())
            ax[figs[ind][0],figs[ind][1]].tick_params(which = 'both')
            for l in k.particles:
                ax[figs[ind][0],figs[ind][1]].scatter(l.x_locations[i],l.y_locations[i], s = 12, c = 'r')
                ax[figs[ind][0],figs[ind][1]].legend()
        
            
import matplotlib.pyplot as plt

# create first galaxy
gal1 = galaxy(0,0, 50, 300, 1e11, name = 'G1')
gal1.createParticle(2, 10)
gal1.createParticle(4, 15)
gal1.createParticle(6, 20)

# create second galaxy
gal2 = galaxy(8, 8, -20, -250, 1e10)
gal2.createParticle(1.5, 8)
gal2.createParticle(2.25, 12)
gal2.createParticle(3.5, 16)

# show just plots of galaxies and particles
plotGal(gal1)
plotGal(gal2)

# progress the system forward in time
# this takes between 10 - 20 minutes to run depending on number of iterations and the step size.
time = progress_system([gal1, gal2], 200, 200000)

# plot the system at different times
plotInteraction([gal1, gal2], [0, 40000, 80000, 120000, 160000,200000], time)