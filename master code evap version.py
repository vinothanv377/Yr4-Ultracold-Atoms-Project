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
np.random.seed(19)

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
        
    def mag_r(self):
        '''Determines the magnitude of the radius for the particle's position'''
        return np.sqrt((self.r[0])**2+(self.r[1])**2+(self.r[2])**2)
    
    def mag_vr(self):
        '''Determines the magnitude of the radial velocity for the particle's position'''
        return np.sqrt((self.v[0])**2+(self.v[1])**2+(self.v[2])**2) #just the speed!
        
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
        
#########################################################################################################
        
        
class Cell: 
    """A class representing a cell to discretise the space"""
    
    def __init__(self, x, y, z): 
        self.stored_atoms = []
        self.centre_x = x
        self.centre_y = y
        self.centre_z = z
        self.centre_coords = np.array([x,y,z])
        
        
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
        self.mass = 85*10**(-27)
      #  self.temp_array = []
        self.space = []
    
    def Create_Particle(self, N):
        '''Function to create N particles and initialise their r and v coordinates. 
        The particles are then appended into a Nx1 dimensional array.'''
        
        for i in range(N):
            test_p = Particle(1,1,1,1,1,1)
            Particle.set_rand_rv(test_p, self.Ti, self.omega_x, self.omega_y, self.omega_z)
            self.particles.append(test_p)
            
            
            
    def Create_Cell(self, Ncell_i):
        #along_L = int(Ncell**(1/3))
        #totalN = Ncell_i**3
        
        for x in range(Ncell_i):
            for y in range(Ncell_i):
                for z in range(Ncell_i):
                    x_pos = -self.L/2 + self.L/(2*Ncell_i) + (x*self.L)/(Ncell_i)
                    y_pos = -self.L/2 + self.L/(2*Ncell_i) + (y*self.L)/Ncell_i
                    z_pos = -self.L/2 + self.L/(2*Ncell_i) + (z*self.L)/Ncell_i
                    test_cell = Cell(x_pos,y_pos,z_pos)
                    self.space.append(test_cell)
                    #print(test_cell.centre_x)
        print(np.size(self.space))

        
    def Check_cells(self):
        fig = plt.figure(figsize=(3,3))
        gs = gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0])
        ax1.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
        ax1.set_xlim(-self.L/2,self.L/2)
        for n in self.space:
            ax1.axhline(y=n.centre_y, xmin=0, xmax=1, linestyle='--', color='g')
            ax1.axvline(x=n.centre_x, ymin=0, ymax= 1, linestyle='--', color='b')
            #print(n.centre_coords)
        
        
    def sort_atoms(self, Ncell_i):
        for c in self.space:
            for n in self.particles:
                if c.centre_x - self.L/(2*Ncell_i) < n.x < c.centre_x + self.L/(2*Ncell_i) and c.centre_y - self.L/(2*Ncell_i) < n.y < c.centre_y + self.L/(2*Ncell_i) and c.centre_z - self.L/(2*Ncell_i) < n.z < c.centre_z + self.L/(2*Ncell_i):
                    c.stored_atoms.append(n)
            print('cell done!')
           
            
    def histogram(self, N, Nt, dt, label):
        '''Function to plot a histogram with associated Gaussian distribution.
        Displays the mean and standard deviation of this Gaussian distribution on the graph.'''
        
        to_plot = [] #empty array
        for i in self.particles:
            to_plot.append(i.x)#Particle.mag_r(i)) #select which values to plot, and append to the array
        fig = plt.figure(figsize=(8,3))
        gs = gridspec.GridSpec(1,1)
        ax2 = fig.add_subplot(gs[0])
        #set title with the simulation parameters to keep track of the data
        ax2.set_title(f'T={self.Ti}, omega={(self.omega_x, self.omega_y, self.omega_z)}, N={N}, Nt={Nt}, dt= {dt}')
        result1 = ax2.hist(to_plot, bins=100) #plotting the histogram
        mean1 = np.mean(to_plot)
        sigma1 = np.sqrt(np.var(to_plot))
        x1 = np.linspace(min(to_plot), max(to_plot), 100)
        dx1 = result1[1][1] - result1[1][0]
        scale1 = len(to_plot)*dx1 #scaling the curves to match the histogram
        
        #plotting a Gaussian of the results, for single direction 
        ax2.plot(x1, scipy.stats.norm.pdf(x1, mean1, sigma1)*scale1, color='g')
        ax2.set_xlabel('x position of particles')
        #ax2.set_ylim(0,300)
        ax2.set_xlim(-0.001,0.001)
        ax2.text(0.0005, 250, 'mean = {:.3g}'.format(mean1))
        ax2.text(0.0005, 200, 'sd = {:.3g}'.format(sigma1))
        ax2.text(-0.0009, max(result1[0])*0.95, '{}'.format(label))
 
        #fitting a maxwell distribution to the data, for radial coords
        #params = scipy.stats.maxwell.fit(to_plot, floc=0) 
        #ax2.plot(x1, scipy.stats.maxwell.pdf(x1, *params)*scale1)
        #ax2.set_xlabel('Radial position of particles')
        #ax2.text(max(result1[1])*0.85, max(result1[0])*0.95, 'mean = {:.3g}'.format(mean1))
        #ax2.text(max(result1[1])*0.85, max(result1[0])*0.75, 'sd = {:.3g}'.format(sigma1))
        #ax2.set_xlim(0,0.06)
        #ax2.set_ylim(0,300)
        
        
    def histogram_TOF(self, N, Nt, dt, label, TOF):
        '''Function to plot a histogram with associated Gaussian distribution.
        Displays the mean and standard deviation of this Gaussian distribution on the graph.'''
        
        to_plot = [] #empty array
        for i in self.particles:
            to_plot.append(Particle.mag_vr(i)) #select which values to plot, and append to the array
        fig = plt.figure(figsize=(8,3))
        gs = gridspec.GridSpec(1,1)
        ax2 = fig.add_subplot(gs[0])
        
        ax2.set_xlabel('x position of particles')
        #set title with the simulation parameters to keep track of the data
        ax2.set_title(f'T={self.Ti}, omega={(self.omega_x, self.omega_y, self.omega_z)}, N={N}, Nt={Nt}, dt= {dt}')
       # ax2.set_xlim(-0.001,0.001)
        result1 = ax2.hist(to_plot, bins=100)
        mean1 = np.mean(to_plot)
        sigma1 = np.sqrt(np.var(to_plot))
        ax2.text(max(result1[1])*0.25, max(result1[0])*0.95, '{}'.format(label))
        ax2.text(max(result1[1])*0.55, max(result1[0])*0.95, 'mean = {:.3g}'.format(mean1))
        ax2.text(max(result1[1])*0.55, max(result1[0])*0.75, 'sd = {:.3g}'.format(sigma1))
        x1 = np.linspace(min(to_plot), max(to_plot), 100)
        dx1 = result1[1][1] - result1[1][0]
        scale1 = len(to_plot)*dx1
        #plotting a Gaussian of the results, display mean and sd of the data on the graph
        ax2.plot(x1, scipy.stats.norm.pdf(x1, mean1, sigma1)*scale1, color='r')
        temp = (self.mass/k)*((self.omega_x*sigma1)**2/(1+(self.omega_x*TOF)**2))
        ax2.text(max(result1[1])*0.55, max(result1[0])*0.25, 'Temp = {:.3g}K'.format(temp))
        
        
        
        

                
                 
    def time_evolve_standard(self,dt,Nt,N):
        '''Runs the simulation through Nt timesteps, of individual size dt. For each timestep
        alters the positions of the particles according to their velocities. Then 
        changes the velocities of the particles to account for the influence of the trapping 
        potential. Plots histograms labelled with the simulation parameters both before and
        after the simulation.
        Can plot 2D or 3D graphs of the particle positions, and save these locally if desired.'''
        
        create_graph2d = False #set to true to plot a xy slice of the particle positions for each time step
        create_graph3d = False #set to true to plot the particle positions in 3D for each time step
        save_images = False #set to true to save the associated images to the specified file path
        create_histograms = True #set to true to plot histograms before and after the time loop
        plot_sigma = False #set to true to plot the standard deviation of r or v over time
        
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
            self.histogram(N, Nt, dt, 't = Ntdt') #create a histogram of the final value
        
        #create graph, and plot the values in the 'sigma' array against time
        if plot_sigma:
            fig = plt.figure(figsize=(8,3))
            gs = gridspec.GridSpec(1,1)
            ax = fig.add_subplot(gs[0])
            ax.set_xlabel('Time /s')
            ax.set_ylabel('Sigma_x /m')
            time = np.arange(len(sigma))*dt
            ax.plot(time, sigma , color='g')
            #print(f'Mean sigma x: {np.mean(sigmax)}')
            #print(f'Peak to peak amplitude (sigma x): {np.max(sigmax)-np.min(sigmax)}')
            #print(f'Mean average sigma: {np.mean(sigma)}')
            #print(f'Peak to peak amplitude (avg sigma): {np.max(sigma)-np.min(sigma)}')
            
##################################################################################################################################################################
            
            
    def time_evolve_TOF(self,dt,Nt,N):
        '''Runs the simulation through Nt timesteps, of individual size dt. For each timestep
        alters the positions of the particles according to their velocities. Then 
        changes the velocities of the particles to account for the influence of the trapping 
        potential. Plots histograms labelled with the simulation parameters both before and
        after the simulation.
        Can plot 2D or 3D graphs of the particle positions, and save these locally if desired.'''
        
        create_graph2d = False #set to true to plot a xy slice of the particle positions for each time step
        create_graph3d = False #set to true to plot the particle positions in 3D for each time step
        save_images = False #set to true to save the associated images to the specified file path
        create_histograms = True #set to true to plot histograms before and after the time loop
        plot_sigma = False #set to true to plot the standard deviation of r or v over time
        
        sigma = [] #creating an empty array to contain the standard deviation values
        
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
                #no potential involved as the trap is off- particles in free expansion
                
            sigmadt = np.sqrt(np.var(sigma_variable)) #calculate the sd of the r/v values for this timestep
            sigma.append(sigmadt) #append it to the empty array created at the start
                
            if save_images: #save the graphs as they are created to the specified file
                 plt.savefig(r'C:\\Users\Bethan\Documents\evaporative cooling\test sim\trapthennot\timestep{t}.png'.format(t=t))
                 
        if create_histograms:
            self.histogram_TOF(N, Nt, dt, 't = Ntdt', Nt*dt) #create a histogram of the final values, and prints calc temp 
        
        #create graph, and plot the values in the 'sigma' array against time
        if plot_sigma:
            fig = plt.figure(figsize=(8,3))
            gs = gridspec.GridSpec(1,1)
            ax = fig.add_subplot(gs[0])
            ax.set_xlabel('Time /s')
            ax.set_ylabel('Sigma_x /m')
            time = np.arange(len(sigma))*dt
            ax.plot(time, sigma , color='g')
            

###################################################################################################################################################################            


    def time_evolve_evap(self,dt,Nt,N,therm_time,rate_of_evap):
        '''Runs the simulation through Nt timesteps, of individual size dt. For each timestep
        alters the positions of the particles according to their velocities. Then 
        changes the velocities of the particles to account for the influence of the trapping 
        potential. Plots histograms labelled with the simulation parameters both before and
        after the simulation.
        This simulation also takes into account evaporative cooling effects and with a defined choice
        of the minimum spactial distribution all particles outside of the minimum are evaporatively
        cooled. therm_time is measure in number of time steps dt how long we leave for rether
        Can plot 2D or 3D graphs of the particle positions, and save these locally if desired.'''
        
        create_graph2d = False #set to true to plot a xy slice of the particle positions for each time step
        create_graph3d = False #set to true to plot the particle positions in 3D for each time step
        save_images = False #set to true to save the associated images to the specified file path
        create_histograms = False #set to true to plot histograms before and after the time loop
        create_histograms_start_finish = True
        plot_sigma = False #set to true to plot the standard deviation of r or v over time
        plot_N = False
        plot_T = False
        
        sigma = [] #creating an empty array to contain the standard deviation values
        sigma_avg = []
        num_atoms = []
        r_vals = [] #defining an empty array for the magnitude radius values
        
        t = float(0) #ensures the value for t is treated as a float to avoid compatibility conditons
        
        for i in self.particles: #looping through all the particles in the system
            r_vals.append(Particle.mag_r(i))  #to append the magnitude radius to the r_vals array

        cut_off = max(r_vals) #setting the cut-off point for minimum radius from the paricles as an arbitrary max value of the r_vals
        evap_times = np.linspace(0, Nt, int(Nt/(therm_time)), endpoint=False)#defining the specific timesteps we want to do evaporative cooling at
        #print(int(Nt/(therm_time-1)))
        #print(f'The Evap Times are {evap_times}')
        #print(f'Starting time t ={t}')
        print(f'Starting Cut-Off is {cut_off}')
        evap_atoms = []
        
        if create_histograms_start_finish:
            self.histogram(N, Nt, dt, 't = 0') #create a histogram of the initial values
        
        while t != Nt: #loop through the timesteps
            sigma_variable = [] #creating an empty array to store the r/v values from which to calc an sd
            sigma_variablex = [] 
            sigma_variabley = []
            sigma_variablez = []
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
                
            for i in self.particles: #for each timestep, loop through the particles in the array
                    #plot the positions if needed
                    if create_graph3d:
                        ax1.scatter3D(i.x, i.y, i.z, c='b', s=2)
                        ax1.view_init(30, 35)
                    if create_graph2d:
                        ax1.scatter(i.x, i.y, c='b', s=2)
                    
                    sigma_variable.append(Particle.mag_r(i))
                    sigma_variablex.append(i.x) #append the r/v values into the empty array
                    sigma_variabley.append(i.y)
                    sigma_variablez.append(i.z)
        
                    #alter the positions of the particles according to their velocities
                    Particle.drift(i, dt) 
                    #alter the velocities of the particles, account for the trapping potential
                    Particle.potential_v_change(i, self.omega_x, self.omega_y, self.omega_z, dt)
                    
            if np.any(t == evap_times): #if statement, a conditon to see t is equal to any of the evap_times
                print(t)
                for i in self.particles: #looping through all partilces
                    if Particle.mag_r(i) >= cut_off: #if statement, a condition to see if magnitude radius of the partilce is greater than or equal to the cut_off radius
                        i.x = 0.004 #if the statment is satisified we set the position of the particle to an arbirtary position, to mimic the effects of evaporative cooling
                        i.y = 0.004 #""
                        i.z = 0.004 #""
                        atom = Particle(i.x, i.y, i.z, i.vx, i.vy, i.vz)
                        evap_atoms.append(atom)
                        self.particles.remove(i)
                   
                    
                cut_off = cut_off*rate_of_evap #outside of this if-else statment, once all particles have been looped/checked against the statment, we redefine the cut-off radius at a rate defined in the funtion variable
                #print(f'Setting new cut-off as {cut_off}')
                if create_histograms:
                    self.histogram(N, Nt, dt, f'{t+1}') #create a histogram of the final value
                if save_images: #save the graphs as they are created to the specified file
                    plt.savefig(r'E:\Yr 4 Cold Atoms Project\03.03.2021 Yr4-Ultracold-Atoms-Project-vinothanv377-patch-1\Yr4-Ultracold-Atoms-Project-vinothanv377-patch-1\images\velocities\timestep{t}.png'.format(t=t))
              
            t+=1 #ensures that time-steps progress by one every full cycle of checks 
                
            sigmaxdt = np.sqrt(np.var(sigma_variablex)) 
            sigmaydt = np.sqrt(np.var(sigma_variabley))
            sigmazdt = np.sqrt(np.var(sigma_variablez))
            sigmadt_avg = (1/3)*(sigmaxdt + sigmaydt + sigmazdt)
            sigmadt = np.sqrt(np.var(sigma_variable)) #calculate the sd of the r/v values for this timestep
            sigma.append(sigmadt) #append it to the empty array created at the start
            sigma_avg.append(sigmadt_avg)
            num_atoms.append(np.size(self.particles))
        
                  
        if create_histograms_start_finish:
            self.histogram(N, Nt, dt, 't = Ntdt') #create a histogram of the final values
            #print(evap_times)
        if save_images: #save the graphs as they are created to the specified file
                 plt.savefig(r'E:\Yr 4 Cold Atoms Project\03.03.2021 Yr4-Ultracold-Atoms-Project-vinothanv377-patch-1\Yr4-Ultracold-Atoms-Project-vinothanv377-patch-1\images\velocities\timestep{t}.png'.format(t=t))
            
        
        #create graph, and plot the values in the 'sigma' array against time
        if plot_sigma:
            fig = plt.figure(figsize=(8,3))
            gs = gridspec.GridSpec(1,1)
            ax = fig.add_subplot(gs[0])
            ax.set_xlabel('Time /s')
            ax.set_ylabel('Sigma_r /m')
            ax.set_title(f'T={self.Ti}, omega={(self.omega_x, self.omega_y, self.omega_z)}, N={N}, Nt={Nt}, dt= {dt}, therm_time= {therm_time}, evap_rate= {rate_of_evap}')
            time = np.arange(len(sigma))*dt
            ax.plot(time, sigma , color='g')
            ax.plot(time, sigma_avg , color='m')
            
           
        if plot_N:
            fig = plt.figure(figsize=(8,3))
            gs = gridspec.GridSpec(1,1)
            ax = fig.add_subplot(gs[0])
            ax.set_xlabel('Time /s')
            ax.set_ylabel('N')
            ax.set_title(f'T={self.Ti}, omega={(self.omega_x, self.omega_y, self.omega_z)}, N={N}, Nt={Nt}, dt= {dt}, therm_time= {therm_time}, evap_rate= {rate_of_evap}')
            time = np.arange(len(sigma))*dt
            ax.plot(time, num_atoms , color='r')
            
        if plot_T:
            fig = plt.figure(figsize=(8,3))
            gs = gridspec.GridSpec(1,1)
            ax = fig.add_subplot(gs[0])
            ax.set_xlabel('Time /ms')
            ax.set_ylabel('T')
            #ax.set_title(f'T={self.Ti}, omega={(self.omega_x, self.omega_y, self.omega_z)}, N={N}, Nt={Nt}, dt= {dt}, therm_time= {therm_time}, evap_rate= {rate_of_evap}')
            time = np.arange(len(sigma))*dt
            ax.plot(evap_times, self.temp_array, color='b')
            
    
    
#########################################################################################
            

        

    
                        
env = Environment(0.01,10**-6,60,60,60)
env.Create_Particle(100)
env.Create_Cell(10)
env.Check_cells()
env.sort_atoms(10)
#env.time_evolve_standard(0.0001, 200, 10000)
#env.time_evolve_evap(0.0001, 200, 10000, 5, 0.97)
#env.time_evolve_TOF(0.0001, 1000, np.size(env.particles))

