# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 15:39:41 2021

@author: Bethan
"""


import numpy as np
from scipy.constants import k
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits import mplot3d
import scipy.stats

#comment this in to set the random number generators, keep it repeatable
#np.random.seed(17)

####################################################################################

class Particle: 
    """A class representing the 3D particle"""
    
    def __init__(self, x, y, z, vx, vy, vz, mass=85*10**(-27), styles=None):
        #mass of Rb-85 atom here 
        
        self.r = np.array((x, y, z))
        self.v = np.array((vx, vy, vz))
        self.mass= mass
    
    #setting things for the positions
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
    
    #setting things for the velocities
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
        
    def drift(self, dt):
        '''Moves the particles through distances appropriate to their velocities.'''
        self.x += self.vx*dt
        self.y += self.vy*dt
        self.z += self.vz*dt
        
    def set_rand_rv(self, T, omega_x, omega_y, omega_z): #here, L is the cube length, T is the temperature
        '''Sets the position components of the particles to a random value between -L/2 and L/2.
        Sets the velocity components of the particles to a random value from a Maxwell-Boltzmann 
        distribution, with mean 0 and standard deviation sqrt(kT/m).'''
        self.v = np.random.normal(0, np.sqrt(k*T/self.mass), 3)
        self.r = np.array((1.0,1.0,1.0)) 
        #so this is a bit weird, not sure why but for some reason the following code sets r components
        #as int unless I put the last line in where the components are given as floats
        self.r[0] = np.random.normal(0, np.sqrt((k*T)/(self.mass*(omega_x**2))))
        self.r[1] = np.random.normal(0, np.sqrt((k*T)/(self.mass*omega_y**2)))
        self.r[2] = np.random.normal(0, np.sqrt((k*T)/(self.mass*omega_z**2)))


    def potential_v_change(self, omega_x, omega_y, omega_z, dt):
        '''Changes the velocity components of the particles by an amount appropriate
        to the 3D harmonic trapping potential.'''
        self.vx += -(omega_x**2)*self.x*dt
        self.vy += -(omega_y**2)*self.y*dt
        self.vz += -(omega_z**2)*self.z*dt
        
        
########################################################################################################        
     
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
    
    
    def Create_Particle(self, N):
        '''Function to create N particles and initialise their r and v coordinates. 
        The particles are then appended into a Nx1 dimensional array.'''
        
        for i in range(N):
            test_p = Particle(1,1,1,1,1,1)
            Particle.set_rand_rv(test_p, self.Ti, self.omega_x, self.omega_y, self.omega_z)
            self.particles.append(test_p)
            
            
    def histogram(self, N, Nt, dt, label):
        '''Function to plot a histogram with associated Gaussian distribution.
        Displays the mean and standard deviation of this Gaussian distribution on the graph.'''
        
        to_plot = [] #empty array 
        for i in self.particles:
            to_plot.append(i.x) #select which values to plot, and append to the array
        fig = plt.figure(figsize=(8,3))
        gs = gridspec.GridSpec(1,1)
        ax2 = fig.add_subplot(gs[0])
        ax2.set_xlim(-0.02,0.02)
        ax2.set_xlabel('x position of particles')
        ax2.text(-0.015, N/40, '{}'.format(label))
        #set title with the simulation parameters to keep track of the data
        ax2.set_title(f'T={self.Ti}, omega={(self.omega_x, self.omega_y, self.omega_z)}, N={N}, Nt={Nt}, dt= {dt}')
        
        result1 = ax2.hist(to_plot, bins=100) 
        mean1 = np.mean(to_plot)
        sigma1 = np.sqrt(np.var(to_plot))
        ax2.text(0.01,N/40, 'mean = {:.3g}'.format(mean1))
        ax2.text(0.01,N/45, 'sd = {:.3g}'.format(sigma1))
        x1 = np.linspace(min(to_plot), max(to_plot), 100)
        dx1 = result1[1][1] - result1[1][0]
        scale1 = len(to_plot)*dx1
        #plotting a Gaussian of the results, display mean and sd of the data on the graph
        ax2.plot(x1, scipy.stats.norm.pdf(x1, mean1, sigma1)*scale1, color='g')
    
            
    def time_evolve(self,dt,Nt,N):
        '''Runs the simulation through Nt timesteps, of individual size dt. For each timestep
        alters the positions of the particles according to their velocities. Then 
        changes the velocities of the particles to account for the influence of the trapping 
        potential. Plots histograms labelled with the simulation parameters both before and
        after the simulation.
        Can plot 2D or 3D graphs of the particle positions, and save these locally if desired.'''
        
        create_graph2d = False #set to true to plot a xy slice of the particle positions for each time step
        create_graph3d = False #set to true to plot the particle positions in 3D for each time step
        save_images = False #set to true to save the associated images to the specified file path
        create_histograms = False #set to true to plot histograms before and after the time loop
        plot_sigma = True #set to true to plot the standard deviation of r or v over time
        
        sigma = [] #creating an empty array to contain the standard deviation values
        
        if create_histograms:
            self.histogram(N, Nt, dt, 't = 0') #create a histogram of the initial values
        
        for t in range(Nt): #loop through the timesteps
            if create_graph3d:
                fig = plt.figure(figsize=(10,10))
                gs = gridspec.GridSpec(1,1)
                ax1 = fig.add_subplot(gs[0], projection='3d') # making the plot 3D capable 
                ax1.set_ylim(-self.L/2,self.L/2) 
                ax1.set_xlim(-self.L/2,self.L/2)
                ax1.set_zlim(-self.L/2,self.L/2)
                ax1.text(self.L/3,-self.L/3, self.L/3, f'{t+1}') # stating on the graph the time-step of the simulation graph
            if create_graph2d:
                fig = plt.figure(figsize=(3,3))
                gs = gridspec.GridSpec(1,1)
                ax1 = fig.add_subplot(gs[0])
                ax1.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
                ax1.set_xlim(-self.L/2,self.L/2)
                ax1.text(self.L/3,-self.L/3, f'{t+1}') #print the timestep number on the graph
                
            sigma_variable = [] #creating an empty array to store the r/v values from which to calc an sd
            
            for i in self.particles: #for each timestep, loop through the particles in the array
                #plot the positions if needed
                if create_graph3d:
                    ax1.scatter3D(i.x, i.y, i.z, c='b', s=2)
                    ax1.view_init(30, 35)
                if create_graph2d:
                    ax1.scatter(i.x, i.y, c='b', s=2)
                    
                sigma_variable.append(i.x) #append the r/v values into the empty array
        
                #alter the positions of the particles according to their velocities
                Particle.drift(i, dt) 
                #alter the velocities of the particles, account for the trapping potential
                Particle.potential_v_change(i, self.omega_x, self.omega_y, self.omega_z, dt)
                
            sigmadt = np.sqrt(np.var(sigma_variable)) #calculate the sd of the r/v values for this timestep
            sigma.append(sigmadt) #append it to the empty array created at the start
                
            if save_images: #save the graphs as they are created to the specified file
                 plt.savefig(r'C:\\Users\Bethan\Documents\evaporative cooling\test sim\trapthennot\timestep{t}.png'.format(t=t))
                 
        if create_histograms:
            self.histogram(N, Nt, dt, 't = Ntdt') #create a histogram of the final values
        
        #create graph, and plot the values in the 'sigma' array against time
        if plot_sigma:
            fig = plt.figure(figsize=(8,3))
            gs = gridspec.GridSpec(1,1)
            ax = fig.add_subplot(gs[0])
            ax.set_xlabel('Time /s')
            ax.set_ylabel('Sigma_x /m')
            time = np.arange(len(sigma))*dt
            ax.plot(time, sigma , color='g')
            #midpoint = (np.max(sigma)-np.min(sigma))/2 +np.min(sigma)
            #ax.axvline(x=time[np.argmax(sigma[0:1000])], ymin=0, ymax= 1, linestyle='--', color='k')
            #ax.axvline(x=time[np.argmax(sigma)], ymin=0, ymax= 1, linestyle='--', color='b')
            #ax.axvline(x=time[np.argmin(sigma)], ymin=0, ymax= 1, linestyle='--', color='r')
            #ax.axhline(y= midpoint, xmin=0, xmax= 1, linestyle='--', color='m')
            #print(f'Max sigma value: {np.max(sigma)} (blue)')
            #print(f'Min sigma value: {np.min(sigma)} (red)')
            #print(f'Oscillating around {midpoint} (pink)')
            #print(f'Time value for first maxima (black): {time[np.argmax(sigma[0:1000])]}')
            #print(f'Time value for maxima (blue): {time[np.argmax(sigma)]}')
            

    
                        
env = Environment(0.01,10**-6,60,60,60)
env.Create_Particle(100)
env.time_evolve(0.0001, 20000, 100)