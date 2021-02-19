# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 15:39:41 2021

@author: Bethan
"""


import numpy as np
#from scipy.constants import k
#import matplotlib 
import matplotlib.pyplot as plt
#import scipy as scipy
#from scipy import optimize
from matplotlib import gridspec
from scipy.constants import k
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
    
    #creates a single particle, random r and v 
    def Create_Particle(self, N):
        for i in range(N):
            test_p = Particle(1,1,1,1,1,1)
            Particle.set_rand_rv(test_p, self.Ti, self.omega_x, self.omega_y, self.omega_z)
            self.particles.append(test_p)
            
    def Make_graph(self):
        fig = plt.figure(figsize=(3,3))
        gs = gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0])
        ax1.set_ylim(-0.5,0.5) #this part is specific to 
        ax1.set_xlim(-0.5,0.5)
        for i in self.particles:
            #print(i.x)
            ax1.scatter(i.x, i.y, c='b')
            
    def time_evolve(self,dt,Nt):
        for t in range(Nt):
            env.Make_graph()
            for i in self.particles:
                Particle.drift(i, dt)
                #Particle.potential_v_change(i, self.omega_x, self.omega_y, self.omega_z, dt)
        
        
    
    
                        
env = Environment(0.1,10**-6,20,20,20)
#env.test_sim(15,10,0.05)
env.Create_Particle(3)
env.Make_graph()
env.time_evolve(1,5)