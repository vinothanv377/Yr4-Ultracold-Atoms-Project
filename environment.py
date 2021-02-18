# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 15:39:41 2021

@author: Bethan
"""


#import numpy as np
#from scipy.constants import k
#import matplotlib 
import matplotlib.pyplot as plt
#import scipy as scipy
#from scipy import optimize
from matplotlib import gridspec
from matplotlib.animation import FuncAnimation
import imageio
import os

from particle import Particle   
     
class Environment: 
    """A class representing the environment to put the particles in"""
    
    def __init__(self, L, Ti, omega_x, omega_y, omega_z): 
    #just got initial temperature and length of cube side for now, can add in more parameters as needed
        self.L = L
        self.Ti = Ti
        self.omega_x = omega_x
        self.omega_y = omega_y
        self.omega_z = omega_z
        self.particles = []
            
    def get_details(self): #this is just to replace the following code until I work out why it's not working
        print(f'Side length is {self.L}')
        print(f'Initial Temperature is {self.Ti}')
        
    #for some reason, the following code kept killing the kernel when I ran it, not sure why
    # @property
    # def Ti(self):
    #    return self.Ti
    # @Ti.setter
    # def Ti(self, value):
    #    self.Ti = value
    #@property
    #def L(self):
    #    return self.L
    #@L.setter
    #def L(self, value):
    #    self.L = value
    
    #creates a single particle, random r and v 
    def Create_Particle(self):
        test_p = Particle(1,1,1,1,1,1)
        Particle.set_rand_rv(test_p,self.L,self.Ti)
        print(test_p.r) #this is mostly a check that it was doing something
  
    #attempting to create a function to replicate the 'plot trajectory' function, 
    #in a more streamlined way ie don't have to type out the thing multiple times,
    #just specify time step, number of particles and total time
    def Create_many_particles(self, N, Nt, dt):
        #this allows multiple particles to be created, can specify the time step and total time loop as well
        #the plotting part was mostly to visually confirm what was happening!
        #How applicable any of this is to the actual simulation is highly questionable
        fig = plt.figure(figsize=(5,5))
        gs = gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0])
        ax1.set_ylim(-0.5,0.5)
        ax1.set_xlim(-0.5,0.5)
        for n in range(N):
            test_p = Particle(1,1,1,1,1,1)
            Particle.set_rand_rv(test_p,self.L,self.Ti)
            #print(test_p.r) #just to check we were looping through okay
            for i in range(Nt):
                #ax1.text(test_p.x,test_p.y, f'{i}') # this part doesn't plot scatter points, it plots the time index where the particle was, mostly for me to see the direction of movement
                ax1.scatter(test_p.x, test_p.y)
                Particle.drift(test_p, dt)
                Particle.potential_v_change(test_p, 2, 2, 2,dt)
    #definitely not the way to actually do the simulation i don't think, it was mainly 
    #just to check everything was working so far
    #THINGS TO FIX- ideally want to have all the dots related to one particle be the
    #same colour, and also the particles aren't trapped in the box!! maybe incorporate 
    #the bouncing in the drift function? as that's the part where they actually move
    
    def test_sim(self, N, Nt, dt):
        '''N is the number of particles you want to simulate, Nt is the number of timesteps to 
        loop through and dt is the size of the timestep.'''
        #okay so I was having issues with the loops using the test_p object so I've resorted to arrays for now and it's clunky as hell- not a fan but there you go
        #completely defeats the point of all the functions I had written in the particle class
        
        r = [[],[],[]] #creating empty arrays to store the r and v values in
        v = [[],[],[]]
        save_images = False #set to true if you want to save all the images as they are made
        create_graph = False #set to true if you want to plot the xy positions for each time step
        
        for n in range(N):
            #looping through the number of particles, creating a particle for each and setting the r and v to be random (gaussian) values
            test_p = Particle(1,1,1,1,1,1)
            Particle.set_rand_rv(test_p,self.Ti, self.omega_x, self.omega_y, self.omega_z)
            #appending the initial r and v values into these arrays
            r[0].append(test_p.r[0])
            r[1].append(test_p.r[1])
            r[2].append(test_p.r[2])
            v[0].append(test_p.v[0])
            v[1].append(test_p.v[1])
            v[2].append(test_p.v[2])
        plt.hist(r[0], bins=100)
        for t in range(Nt): #looping through all the timesteps
            #create a graph
            if create_graph:
                fig = plt.figure(figsize=(3,3))
                gs = gridspec.GridSpec(1,1)
                ax1 = fig.add_subplot(gs[0])
                ax1.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
                ax1.set_xlim(-self.L/2,self.L/2)
                ax1.text(self.L/3,-self.L/3, f'{t+1}') #print the timestep number on the graph
            for n in range(N): #loop through all the particles for the specific timestep
                #plot the x y positions of the particles
                if create_graph:
                    ax1.scatter(r[0][n], r[1][n], c='b', s=2)
                #these following lines play the role of the drift function- updates r values in the arrays
                r[0][n] += v[0][n]*dt
                r[1][n] += v[1][n]*dt
                r[2][n] += v[2][n]*dt
                #these following lines play the role of the potential_v_drift function- updates v values in the arrays
                v[0][n] += -(self.omega_x**2)*r[0][n]*dt 
                v[1][n] += -(self.omega_y**2)*r[1][n]*dt
                v[2][n] += -(self.omega_z**2)*r[2][n]*dt
            if save_images: #this will save all the images into the specified file, labelled by their timesteps
                plt.savefig(r'C:\\Users\Bethan\Documents\evaporative cooling\test sim\images\timestep{t}.png'.format(t=t))
        
        
        #attempting to plot a histogram of the last value
        plt.hist(r[0], bins=50)
                        
env = Environment(0.1,10**-6,20,20,20) #here L= 10cm, T = 1uK, note these freq are the omega, included the 2pi factors 
env.test_sim(100,10,0.05)
