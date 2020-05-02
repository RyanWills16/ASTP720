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
        self.vx = vx
        self.vy = vy
        
        self.x_locations = [self.x]
        self.y_locations = [self.y]
        
        
        
class galaxy:
    
    def __init__(self, x, y, vx, vy, Mass, name = ''):
        
        # position in kpc
        self.x = x
        self.y = y
        
        # velocity in kpc/s
        self.vx = vx*(1/3.086e16)
        self.vy = vy*(1/3.086e16)
        
        self.Mass = Mass
        
        self.name = name
        
        self.particles = []
        
        self.x_locations = [self.x]
        self.y_locations = [self.y]
    
    def createParticle(self, R, num):
        import numpy as np
        
        G = 4.513e-39 # (kpc^3) / (s^2 solar mass)
        theta = (2*np.pi)/num
        
        for i in range(0, num+1):
            angle = i*theta
            
            x = R*np.cos(angle) + self.x
            y = R*np.sin(angle) + self.y
            
            vx = np.sqrt(G * self.Mass / R)*np.sin(angle) + self.vx
            vy = np.sqrt(G * self.Mass / R)*np.cos(angle) + self.vy
            
            self.addParticle(particle(x,y,vx,vy))
    
            
    def addParticle(self, particle):
        
        self.particles.append(particle)
        return

def calc_acceleration(galaxy_list, body, time_step):
    import numpy as np
    
    accelx = 0
    accely = 0
    
    for i in galaxy_list:
        
        if i == body:
            continue
        
        
        else:
            G = 4.513e-39 # (kpc^3) / (s^2 solar mass)
            r = np.sqrt((i.x - body.x)**2 + (i.y - body.y)**2)
            
            C = G * i.Mass / r**3
            
            # k1
            k1x = C * (i.x - body.x)
            k1y = C * (i.y - body.y)
            
            # print("k1 = ", k1x)
            
            # k2
            vel2x = body.vx + k1x*0.5*time_step
            vel2y = body.vy + k1y*0.5*time_step
            
            p2x = body.x + vel2x*0.5*time_step
            p2y = body.y + vel2y*0.5*time_step
            
            k2x = (i.x - p2x)*C
            k2y = (i.y - p2y)*C
            
            # print("k2 = ", k2x)
            
            # k3
            vel3x = body.vx + k2x*0.5*time_step
            vel3y = body.vy + k2y*0.5*time_step
            
            p3x = body.x + vel3x*0.5*time_step
            p3y = body.y + vel3y*0.5*time_step
            
            k3x = (i.x - p3x)*C
            k3y = (i.y - p3y)*C
            # print("k3 = ", k3x)
            
            # k4
            vel4x = body.vx + k3x*time_step
            vel4y = body.vy + k3y*time_step
            
            p3x = body.x + vel4x*time_step
            p3y = body.y + vel4y*time_step
            
            k4x = (i.x - p3x)*C
            k4y =  (i.y - p3y)*C
            # print("k4 = ", k4x)
            # acceleration
            accelx += (k1x + k2x + k3x + k4x) / 6
            accely += (k1y + k2y + k3y + k4y) / 6
        
    return accelx, accely

def calc_velocity(galaxy_list, body, time_step):
    
    accelerationx, accelerationy = calc_acceleration(galaxy_list, body, time_step)
    # print("vx = ", body.vx)
    # print("vy = ", body.vy)
    
    body.vx += accelerationx * time_step
    body.vy += accelerationy * time_step
    
    # print("name = ", body.name, '\n')
    # print(accelerationx, accelerationy, body.vx, body.vy)
    
def update_position(body, time_step):
    # print("xb = ", body.x)
    # print("yb = ", body.y)
    
    body.x_locations.append(body.x + body.vx*time_step)
    body.x += body.vx*time_step  # convert v to kpc/s
    
    body.y_locations.append(body.y + body.vy*time_step)
    body.y += body.vy*time_step  # convert v to kpc/s
    
    # print("xa = ", body.x)
    # print("ya = ", body.y)

def all_calculations(galaxy_list, body, time_step):
    time_step = time_step * 3.154e7 # convert years to seconds
    
    calc_velocity(galaxy_list, body, time_step)
    
    update_position(body, time_step)
    
    
def progress_system(galaxy_list,  time_step, num_step):
    
    print(f"This will progress the system by {format(time_step*num_step, '.3e')} years")
    
    n = 0
    
    bodies = []
    
    for i in galaxy_list:
        bodies.append(i)
        
        if len(i.particles) == 0:
            continue
        
        else:
            for j in i.particles:
                bodies.append(j)
                
    # print('systems = ', bodies)
    
    while n < num_step:
        n += 1
        
        for i in bodies:
            all_calculations(galaxy_list, i, time_step)
        
            
import matplotlib.pyplot as plt

t = galaxy(1,1, 0, 0, 1e11, name = 'G1')

t.createParticle(2, 1)
# t.createParticle(2.3, 20)
# t.createParticle(3.5, 25)

# t2 = galaxy(10,10, 0, -100, 1e11, name = 'G2')

progress_system([t], 100, 200000)

# plt.figure()
# # plt.plot(t2.x_locations, t2.y_locations, 'blue', zorder = 1)
# plt.scatter(1, 1)
# # plt.scatter(t2.x_locations[-1], t2.y_locations[-1], s =20, c= 'black', zorder = 2)
# plt.scatter(t.particles[0].x, t.particles[0].y)

# # plt.scatter(t2.x_locations[0], t2.y_locations[0], s =20, c= 'black', marker = 'v', zorder = 2, label = 'start')
# plt.scatter(t.x_locations[0], t.y_locations[0], s = 20, c='black', marker = '^', label = 'start')

# plt.axvline(x = 0, c = 'gray', zorder = 0)
# plt.axhline(y = 0, c = 'gray', zorder = 0)
# # plt.legend()

def plotGal(galaxy):
    plt.figure()
    plt.scatter(t.x, t.y, s = 20, c='black', zorder = 2)
    plt.axvline(x = 0, c = 'gray', zorder = 0)
    plt.axhline(y = 0, c = 'gray', zorder = 0)
    for i in galaxy.particles:
        plt.scatter(i.x_locations, i.y_locations, s = 10, c = 'r')


plotGal(t)



# from matplotlib import animation

# fig = plt.figure()
# ax = plt.axes(xlim=(-10, 10), ylim = (-10, 10))
# line, = ax.plot([],[], lw = 2)

# def initialize():
#     line.set_data([],[])
#     return line,

# def animate(i):
#     x = t.x_locations
#     y = t.y_locations
    
#     line.set_data(x,y)
#     return line,

# anim = animation.FuncAnimation(fig, animate, init_func = initialize, interval = 10000)

# plt.show()

# plotmovie(t.x_locations, t.y_locations, -10, 10, -10, 10)


    
