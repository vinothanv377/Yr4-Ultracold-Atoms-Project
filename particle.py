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
    
    def __init__(self, x, y, z, vx, vy, vz, radius=0.01, mass=85*10**(-26), styles=None):
        #mass of Rb-85 atom here 
        
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
        '''Moves the particles through distances appropriate to their velocities.'''
        self.x += self.vx*dt
        self.y += self.vy*dt
        self.z += self.vz*dt
        #put in the box constraints here? REMEMBER BETHAN- initial pos are from -L/2 to L/2
        
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
        self.vx += -(omega_x**2)*self.x*dt
        self.vy += -(omega_y**2)*self.y*dt
        self.vz += -(omega_z**2)*self.z*dt
        
        
#test_p = Particle(1,1,1,1,1,1)
#Particle.set_rand_rv(test_p, 10**-6, 20, 20, 20)

        
