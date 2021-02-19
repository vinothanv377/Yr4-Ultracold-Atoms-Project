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
from mpl_toolkits import mplot3d


class Particle: 
    """A class representing the 3D particle"""
    
    def __init__(self, x, y, z, vx, vy, vz, radius=0.01, mass=10**(-26), styles=None):
        
        self.r = np.array((x, y, z))
        self.v = np.array((vx, vy, vz))
        self.radius = radius
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
        self.x += self.vx*dt
        self.y += self.vy*dt
        self.z += self.vz*dt
        
    def set_rand_rv(self, L, T):
        #here, L is the cube length, T is the temperature
        #self.r = [0.01,0.01,0.2]
        self.r = (np.random.random(3)-0.5)*L #sets each r value to a random value between -L/2 and L/2
        self.v = np.random.normal(0, np.sqrt(k*T/self.mass), 3) 
        #setting velocities to random values selected from a Gaussian (Maxwell-Boltzmann) distribution
        #mean 0, standard deviation sqrt(kT/m) 
        #also thinking maybe this should be in another class? Like 'Environment' or something,
        #where we specify the cube length
 
    def potential_v_change(self, omega_x, omega_y, omega_z, dt):
        self.vx += -(omega_x**2)*self.x*dt
        self.vy += -(omega_y**2)*self.y*dt
        self.vz += -(omega_z**2)*self.z*dt
        #not so sure about this part, the changes to velocity are super big and kinda gross-
        #could be that my parameters are just way off??
        


     
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
    
    def Create_many_particles_3D(self, N, Nt, dt):
        #this allows multiple particles to be created, can specify the time step and total time loop as well
        #the plotting part was mostly to visually confirm what was happening!
        #How applicable any of this is to the actual simulation is highly questionable
        fig = plt.figure(figsize=(10,10))
        gs = gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0], projection='3d')
        #ax1.set_ylim(-0.5,0.5)
        #ax1.set_xlim(-0.5,0.5)
        #ax1.set_zlim(-0.5,0.5)
        
        
        '''--------Plotting the harmonic potential-------------'''
        def quad_func(x, y, z, w_x, w_y, w_z):
            ''' a function to describe the quadratic hamronic potential'''
            m=10**-27
            return (1/2)*m*(w_x**2*x**2 +w_y**2*y**2 +w_z**2*z**2)
        
        x = np.linspace(-0.5, 0.5, 30)# defining the range of x values
        y = np.linspace(-0.5, 0.5, 30)
        z = np.linspace(-0.5, 0.5, 30)
        #print(z.ndim)# defiening the range of y values
        X ,Y, Z = np.meshgrid(x, y, z)
        #print(X)
        #print(Z.ndim)
        #print(Y)
        #print(z)#produces a matrix of the x and y coordiantes for the 3D plot
        U = (quad_func(X, Y, Z, 20, 20, 20)) #
        har_pot = ax1.scatter3D(X, Y, Z, c=U)
        fig.colorbar(har_pot, shrink=0.5, aspect=5)
        '''----------------------------------------------------'''
        
        
        for n in range(N):
            test_p = Particle(1,1,1,1,1,1)
            Particle.set_rand_rv(test_p,self.L,self.Ti)
            print(test_p.r)
            print(test_p.v)#just to check we were looping through okay
            for i in range(Nt):
                #ax1.text(test_p.x,test_p.y, f'{i}') # this part doesn't plot scatter points, it plots the time index where the particle was, mostly for me to see the direction of movement
                #ax1.scatter3D(test_p.x, test_p.y, test_p.z)
                #ax1.view_init(0, 90)
                ax1.view_init(20,35)
                #ax1.view_init(60, 65)
                Particle.drift(test_p, dt)
                Particle.potential_v_change(test_p, 2, 2, 2,dt)
    #definitely not the way to actually do the simulation i don't think, it was mainly 
    #just to check everything was working so far
    #THINGS TO FIX- ideally want to have all the dots related to one particle be the
    #same colour, and also the particles aren't trapped in the box!! maybe incorporate 
    #the bouncing in the drift function? as that's the part where they actually move
                        
env = Environment(1,10**-6)
#env.Create_many_particles(1,100,0.05)
env.Create_many_particles_3D(1,100,0.05)