# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 17:08:24 2021

@author: Bethan
"""

import numpy as np
from scipy.constants import k
import matplotlib
import matplotlib.pyplot as plt
import scipy as scipy
from scipy import optimize
from matplotlib import gridspec

class Particle: 
    """A class representing the 3D particle"""
    
    def __init__(self, x, y, z, vx, vy, vz, radius=0.01, mass=10**(-26), styles=None):
        
        self.r = np.array((x, y, z)) # array that stores the current postion of the particle
        self.v = np.array((vx, vy, vz))# array that stores the current velocity of the particle
        self.radius = radius
        self.mass= mass
    
    #functions for setting things for the positions
    @property
    def x(self):
        return self.r[0]
    @x.setter
    def x(self, value):
        self.r[0] = value
    @property
    def y(self):
        return self.r[1]
    @y.setter
    def y(self, value):
        self.r[1] = value
    @property
    def z(self):
        return self.r[2]
    @z.setter
    def z(self, value):
        self.r[2] = value
    
    #fucntions for setting things for the velocities
    @property
    def vx(self):
        return self.v[0]
    @vx.setter
    def vx(self, value):
        self.v[0] = value
    @property
    def vy(self):
        return self.v[1]
    @vy.setter
    def vy(self, value):
        self.v[1] = value
    @property
    def vz(self):
        return self.v[2]
    @vz.setter
    def vz(self, value):
        self.v[2] = value
    
    # the change in position for the particle due to drift movement after time step   
    def drift(self, dt):
        self.x += self.vx*dt
        self.y += self.vy*dt
        self.z += self.vz*dt
        
    def set_rand_rv(self, L, T): #here, L is the cube length, T is the temperature
        self.r = np.random.random(3)*L #sets each r value to a random value between 0 and L
        self.v = np.random.normal(0, np.sqrt(k*T/self.mass), 3)
        #setting velocities to random values selected from a Gaussian (Maxwell-Boltzmann) distribution, mean 0, standard deviation sqrt(kT/m) 
        #also thinking maybe this should be in another class? Like 'Environment' or something, where we specify the cube length
 
    def potential_v_change(self, omega_x, omega_y, omega_z, dt):
        self.vx += -(omega_x**2)*self.x*dt
        self.vy += -(omega_y**2)*self.y*dt
        self.vz += -(omega_z**2)*self.z*dt
        #not so sure about this part, the changes to velocity are super big and kinda gross- could be that my parameters are just way off??
        
#well this part i've commented out was me starting to try to build a box to put the particles in- fully killed the kernel though so I think it's probably not quite right
        
#class Environment: 
#    """A class representing the environment to put the particles in"""
#    
#    def __init__(self, L, Ti): 
#    #just got initial temperature and length of cube side for now, can add in more parameters as needed
#        self.L = L
#        self.Ti = Ti
        
#    @property
#    def Ti(self):
#        return self.Ti
#    @Ti.setter
#    def Ti(self, value):
#        self.Ti = value
#    @property
#    def L(self):
#        return self.L
#    @L.setter
#    def L(self, value):
#        self.L = value

#    def Create_Particle(self):
#        test_p = Particle(1,1,1,1,1,1)
#        Particle.set_rand_rv(test_p,self.L,self.T)
#        print(test_p.r)
        
        
#okay this stuff from here on out is just testing what's going on!
        
test_p = Particle(1,2,3,4,4,4)     
Particle.set_rand_rv(test_p,1,10**(-6))  
test_p1 = Particle(1,2,3,4,4,4)     
Particle.set_rand_rv(test_p1,1,10**(-6))
#creating 2 test particles

tau = 1 #mean free time, would normally calculate using the density etc
N_tau = 5 #number of mean free times we want to run through
Nt = N_tau*1 #total number of time steps, ie whatever we put in time steps per mean free time
dt = N_tau*tau/Nt #length of time step
#tau*N_tau is the same as Nt*dt- ie total amount of time

#just creating a function that'll plot the positions of the two particles over the time we specified
#mostly to test that the trajectory is straight and that the faster particle does go further in the time
def Plot_trajectory(Nt):
    fig = plt.figure(figsize=(5,5))
    gs = gridspec.GridSpec(1,1)
    ax1 = fig.add_subplot(gs[0])
    for i in range(Nt):
        ax1.scatter(test_p.x, test_p.y, c='b')
        Particle.drift(test_p, dt)
        ax1.scatter(test_p1.x, test_p1.y, c='r')
        Particle.drift(test_p1, dt)
    print(test_p.v)
    print(test_p1.v)
    print(np.sqrt(test_p.vx**2 + test_p.vy**2 + test_p.vz**2))
    print(np.sqrt(test_p1.vx**2 + test_p1.vy**2 + test_p1.vz**2))
          
Plot_trajectory(Nt)