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
#import imageio
import os
from mpl_toolkits import mplot3d
import scipy.stats
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
            
    def histogram(self):
        '''Function to plot the histogram'''
        to_plot = [] #empty array 
        for i in self.particles:
            to_plot.append(i.vy)# appending the values to an array
        fig, (ax2, ax3) = plt.subplots(2)
        result1 = ax2.hist(to_plot, bins=100)
        ax2.set_xlim(-0.02,0.02)
        mean1 = np.mean(to_plot)
        sigma1 = np.sqrt(np.var(to_plot))
        ax2.text(0.01,10, 'mean {:.3g}'.format(mean1))
        ax2.text(0.01,20, 'sd {:.3g}'.format(sigma1))
        #ax2.set_title(f'N = {N}, Nt ={Nt}, dt={dt} T={self.Ti}, omega={self.omega_x}')
        x1 = np.linspace(min(to_plot), max(to_plot), 100)
        dx1 = result1[1][1] - result1[1][0]
        scale1 = len(to_plot)*dx1
        #plotting a Gaussian of the results, display mean and sd of the data on the graph
        ax2.plot(x1, scipy.stats.norm.pdf(x1, mean1, sigma1)*scale1, color='g')
    
    '''def Make_graph(self):
        fig = plt.figure(figsize=(3,3))
        gs = gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0], projection='3d')
        ax1.set_ylim(-0.5,0.5) #this part is specific to 
        ax1.set_xlim(-0.5,0.5)
        
        for i in self.particles:
            #print(i.x)
            ax1.scatter3D(i.x, i.y, c='b', s=2)'''
            
    def time_evolve(self,dt,Nt):
        self.histogram()
        for t in range(Nt):
            #fig = plt.figure(figsize=(3,3))
            #gs = gridspec.GridSpec(1,1)
            #ax1 = fig.add_subplot(gs[0], projection='3d') # making the plot 3D capable 
            #ax1.set_ylim(-0.1,0.1) #this part is specific to 
            #ax1.set_xlim(-0.1,0.1)
            #ax1.set_zlim(-0.1,0.1)
            #ax1.text(self.L/3,-self.L/3, self.L/3, f'{t+1}') # stating on the graph the time-step of the simulation graph
            for i in self.particles:
                #print(i.x)
                #ax1.scatter3D(i.x, i.y, i.z, c='b', s=2)
                #ax1.view_init(0, 60)
                Particle.drift(i, dt)
                Particle.potential_v_change(i, self.omega_x, self.omega_y, self.omega_z, dt)
        self.histogram()
        
        
    
    
                        
env = Environment(0.1,10**-6,20,20,20)
#env.test_sim(15,10,0.05)
env.Create_Particle(1000)
env.time_evolve(0.01,20)