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
        
    def set_rand_rv(self, L, T): #here, L is the cube length, T is the temperature
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
        

        
