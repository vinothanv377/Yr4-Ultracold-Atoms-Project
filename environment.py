# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 15:39:41 2021

@author: Bethan
"""


import numpy as np
from scipy.constants import k
import matplotlib
import matplotlib.pyplot as plt
import scipy as scipy
from scipy import optimize
from matplotlib import gridspec

from particle import Particle   
     
class Environment: 
    """A class representing the environment to put the particles in"""
    
    def __init__(self, L, Ti): 
    #just got initial temperature and length of cube side for now, can add in more parameters as needed
        self.L = L
        self.Ti = Ti
        
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
            print(test_p.r) #just to check we were looping through okay
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
                        
env = Environment(1,10**-6)
env.Create_many_particles(1,100,0.05)