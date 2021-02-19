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
            
        #attempting to plot a histogram of the initial values (for x at the minute)
        to_plot = v[0]
        fig, (ax2, ax3) = plt.subplots(2)
        result1 = ax2.hist(to_plot, bins=100)
        ax2.set_xlim(-0.02,0.02)
        mean1 = np.mean(to_plot)
        sigma1 = np.sqrt(np.var(to_plot))
        ax2.text(0.01,10, 'mean {:.3g}'.format(mean1))
        ax2.text(0.01,20, 'sd {:.3g}'.format(sigma1))
        ax2.set_title(f'N = {N}, Nt ={Nt}, dt={dt} T={self.Ti}, omega={self.omega_x}')
        x1 = np.linspace(min(to_plot), max(to_plot), 100)
        dx1 = result1[1][1] - result1[1][0]
        scale1 = len(to_plot)*dx1
        #plotting a Gaussian of the results, display mean and sd of the data on the graph
        ax2.plot(x1, scipy.stats.norm.pdf(x1, mean1, sigma1)*scale1)
        
        
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
        
        
        #attempting to plot a histogram of the last values (for x at the minute)
        result = ax3.hist(to_plot, bins=100)
        ax3.set_xlim(-0.02,0.02)
        #plotting a Gaussian of the results, display mean and sd of the data on the graph
        mean = np.mean(to_plot)
        sigma = np.sqrt(np.var(to_plot))
        ax3.text(0.01,10, 'mean {:.3g}'.format(mean))
        ax3.text(0.01,20, 'sd {:.3g}'.format(sigma))
        x = np.linspace(min(to_plot), max(to_plot), 100)
        dx = result[1][1] - result[1][0]
        scale = len(to_plot)*dx
        ax3.plot(x, scipy.stats.norm.pdf(x, mean, sigma)*scale)
        


                
env = Environment(0.1,10**-6,40,40,40) #here L= 10cm, T = 1uK, note these freq are the omega, included the 2pi factors 
env.test_sim(1000,50,0.01)
