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
from matplotlib.patches import Rectangle
import random


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
        self.coll_atoms = [] 
        self.ens_avg_N_array = [] #array to store the number of atoms in the cell over each timestep, used to determine the time averaged value for N in this cell.
        self.subcellno = 0
        self.particleno = []
        self.subcellspace = []
        
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
        self.space = [] #array to store the subcells
        self.cross_sect = 10**-16 #just placeholder values- this is from 4*pi*a^2 using a for Rb 
        self.vr_max = 0.05 #from those very brief tests
        self.N_collisions = [] #an empty array to store the total number of collisions per timestep
        
    
    def Create_Particle(self, N):
        '''Function to create N particles and initialise their r and v coordinates. 
        The particles are then appended into a Nx1 dimensional array.'''
        
        for i in range(N):
            test_p = Particle(1,1,1,1,1,1)
            Particle.set_rand_rv(test_p, self.Ti, self.omega_x, self.omega_y, self.omega_z)
            self.particles.append(test_p)
            
            
            
    def Create_Cell(self, Ncell_i):
        for x in range(Ncell_i):
            for y in range(Ncell_i):
                for z in range(Ncell_i):
                    x_pos = -self.L/2 + self.L/(2*Ncell_i) + (x*self.L)/(Ncell_i)
                    y_pos = -self.L/2 + self.L/(2*Ncell_i) + (y*self.L)/Ncell_i
                    z_pos = -self.L/2 + self.L/(2*Ncell_i) + (z*self.L)/Ncell_i
                    test_cell = Cell(x_pos,y_pos,z_pos)
                    self.space.append(test_cell)
        print(np.size(self.space))

        
    def Check_cells(self):
        fig = plt.figure(figsize=(3,3))
        gs = gridspec.GridSpec(1,1)
        ax1 = fig.add_subplot(gs[0])
        ax1.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
        ax1.set_xlim(-self.L/2,self.L/2)
        for n in self.space: #plots lines that intersect at the centre of the cells
            ax1.axhline(y=n.centre_y, xmin=0, xmax=1, linestyle='--', color='g')
            ax1.axvline(x=n.centre_x, ymin=0, ymax= 1, linestyle='--', color='b')
            
        
    def sort_atoms(self, Ncell_i): #initially sorting the atoms into the appropriate cells
        for c in self.space:
            for n in self.particles:
                if c.centre_x - self.L/(2*Ncell_i) < n.x < c.centre_x + self.L/(2*Ncell_i) and c.centre_y - self.L/(2*Ncell_i) < n.y < c.centre_y + self.L/(2*Ncell_i) and c.centre_z - self.L/(2*Ncell_i) < n.z < c.centre_z + self.L/(2*Ncell_i):
                    c.stored_atoms.append(n)

            
    def sort_atoms_again(self, Ncell_i): #sorting the atoms again, doesn't append them so total atom no. remains constant
        for c in self.space:
            for n in self.particles:
                c.stored_atoms = [n for n in self.particles if c.centre_x - self.L/(2*Ncell_i) < n.x < c.centre_x + self.L/(2*Ncell_i) and c.centre_y - self.L/(2*Ncell_i) < n.y < c.centre_y + self.L/(2*Ncell_i) and c.centre_z - self.L/(2*Ncell_i) < n.z < c.centre_z + self.L/(2*Ncell_i)]
            c.particleno.append(np.size(c.stored_atoms)) #append the number of particles present in this cell in this timestep to this array
                
    

    def create_subcells(self, Ncell_i, t):
        countatomsappend= []
        countatomsloop= [] #again these arrays are just to check
        
        self.sort_atoms_again(Ncell_i) #sort atoms into the big cells- should also update the cell population no 
        
        for c in self.space: #loop through all the large cells   
            #Only make new cells if we needed to- would need to access infomation about previous timestep
            if t==0 or c.particleno[int(t-1)] != c.particleno[int(t)]: 
                c.subcellspace.clear() #only clear subcells and make new ones if the number of atoms in cell is different than it was before, otherwise just use the old cells
                
                if np.size(c.stored_atoms) > 20: #determining no. of subcells needed based on how many atoms in cell
                    c.subcellno = 5
                elif 20 >= np.size(c.stored_atoms) > 15:
                    c.subcellno = 4
                elif 15 >= np.size(c.stored_atoms) > 10:
                    c.subcellno = 3
                elif 10 >= np.size(c.stored_atoms) > 5:
                    c.subcellno = 2
                elif 5 >= np.size(c.stored_atoms)>= 0:
                    c.subcellno = 1 
                
                #then need to actually create the subcells using these allocations
                for x in range(c.subcellno):
                    for y in range(c.subcellno):
                        for z in range(c.subcellno):
                            x_pos = c.centre_x - self.L/(2*Ncell_i) +self.L/(2*Ncell_i*c.subcellno) + (x*self.L)/(Ncell_i*c.subcellno)
                            y_pos = c.centre_y - self.L/(2*Ncell_i) +self.L/(2*Ncell_i*c.subcellno)+ (y*self.L)/(Ncell_i*c.subcellno)
                            z_pos = c.centre_z - self.L/(2*Ncell_i) +self.L/(2*Ncell_i*c.subcellno)+ (z*self.L)/(Ncell_i*c.subcellno)
                            test_cell = Cell(x_pos,y_pos,z_pos)
                            c.subcellspace.append(test_cell) #appending these subcells to the subcell array not the big cell one
            
                 
                for s in c.subcellspace: #loop through all subcells in the cell, don't need to clear as these are fresh subcells
                    for n in c.stored_atoms: #loop through atoms in the specific cell we're considering
                        if s.centre_x - self.L/(2*Ncell_i*c.subcellno) < n.x <= s.centre_x + self.L/(2*Ncell_i*c.subcellno)  and s.centre_y - self.L/(2*Ncell_i*c.subcellno) < n.y <= s.centre_y + self.L/(2*Ncell_i*c.subcellno) and s.centre_z - self.L/(2*Ncell_i*c.subcellno) < n.z <= s.centre_z + self.L/(2*Ncell_i*c.subcellno):
                            s.stored_atoms.append(n) #append the atom into the subcell's stored atoms array if it's in the subcell

                    countatomsappend.append(np.size(s.stored_atoms)) #just here to checl

                    
            else: #ie if we're not making fresh subcells
                for s in c.subcellspace: #loop through all subcells in the cell
                    s.stored_atoms.clear() #clear the subspace array to solve the over counting iss
                    for n in c.stored_atoms: #loop through atoms in the specific cell we're considering
                        if s.centre_x - self.L/(2*Ncell_i*c.subcellno) < n.x <= s.centre_x + self.L/(2*Ncell_i*c.subcellno)  and s.centre_y - self.L/(2*Ncell_i*c.subcellno) < n.y <= s.centre_y + self.L/(2*Ncell_i*c.subcellno) and s.centre_z - self.L/(2*Ncell_i*c.subcellno) < n.z <= s.centre_z + self.L/(2*Ncell_i*c.subcellno):
                                s.stored_atoms.append(n) #append atoms into the right subcells
                    
                    countatomsloop.append(np.size(s.stored_atoms)) #just here to check
              
        print(f'APPEND  atoms as counted in the subcell array: {np.sum(countatomsappend)}')
        print(f'entries in subcell count atoms array: {np.size(countatomsappend)}')
        print(f'LOOP    atoms as counted in the subcell array: {np.sum(countatomsloop)}')
        print(f'entries in subcell count atoms array: {np.size(countatomsloop)}')
        

    def check_subcells(self, Ncell_i):
        counting_subcells = []
        #fig = plt.figure(figsize=(3,3))
        #gs = gridspec.GridSpec(1,1)
        #ax1 = fig.add_subplot(gs[0])
        #ax1.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
        #ax1.set_xlim(-self.L/2,self.L/2)
        #print(self.space)
        #b = 3
        #for c in self.space[b].subcellspace:
        #    ax1.scatter(c.centre_x, c.centre_y)
        for x in self.space:
            nosubcells = np.size(x.subcellspace)
            counting_subcells.append(nosubcells)
        print(f'no of entries in array: {np.size(counting_subcells)}')
        print(f'no of subcells total: {np.sum(counting_subcells)}')
        
        #ax1.add_patch(Rectangle((self.space[b].centre_x - self.L/(2*Ncell_i), self.space[b].centre_y - self.L/(2*Ncell_i)), self.L/(Ncell_i), self.L/(Ncell_i), fill=False, ls='--'))
        #ax1.text(0,-1.5*(self.L/2), f'atoms = {np.size(self.space[b].stored_atoms)}')
        #ax1.text(0,-1.5*(self.L/2)-0.0003, f'subcells = {np.size(self.space[b].subcellspace)}')
        #cenx = self.space[b].subcellspace[0].centre_x - self.L/(2*Ncell_i*self.space[b].subcellno) #+self.L/(2*Ncell_i*self.space[b].subcellno) 
        #ceny = self.space[b].subcellspace[0].centre_y - self.L/(2*Ncell_i*self.space[b].subcellno) #+self.L/(2*Ncell_i*self.space[b].subcellno) 
        #ax1.add_patch(Rectangle((cenx, ceny), self.L/(Ncell_i*self.space[b].subcellno), self.L/(Ncell_i*self.space[b].subcellno), fill=False, ls='--'))
        
        
        
##########################################################################################################################################################################################################################################################################################################
    
    def rel_vel_final(self, atom_1, atom_2, c_m, c_r_mag):
        '''determining the final velocities of the 2 particles after their collision'''
        m1 = self.mass
        m2 = self.mass
        n_vector = np.random.rand(3) #producing a randomly generated post-collision vector in the CoM frame
        unit_n_vector = n_vector/np.linalg.norm(n_vector) #normalising the vector so that it would be a vector produced on the unit sphere
        #atom_1.v = c_m + (m2/(m1 + m2))*c_r*unit_n_vector #determining the post velocities back into general frame of refernce from the CoM frame
        #atom_2.v = c_m - (m1/(m1 + m2))*c_r*unit_n_vector 
        atom_1.v = c_m + (m2/(m1 + m2))*c_r_mag*unit_n_vector #determining the post velocities back into general frame of refernce from the CoM frame
        atom_2.v = c_m - (m1/(m1 + m2))*c_r_mag*unit_n_vector
        
        return atom_1.v, atom_2.v #outputting the velocities of 1st and 2nd atom in general frame 
        
        
    
    def collisons(self, dt, Ncell_i):
        ''' function that will evalute the possibility and outcomes of two-body particle collisons in a cell over 
        all cells in the environment'''
        
        m1 = self.mass #mass of atom 1
        m2 = self.mass #mass of atom 2
        V_c = (self.L/Ncell_i)**3 #volume of cell
        N_coll_dt = [] #create empty array to store the collisions for each cell for this timestep
        
        for c in self.space: #looping over all cells in the environment i.e. space
            n_coll_cell = 0 #resetting the collision counter for the new cell
            c.coll_atoms.clear() #clears the array for any particles stored in it from the last time-step's collided atoms
            
            print(np.size(c.stored_atoms))
            
            if np.size(c.stored_atoms) > 1: #minimum condition needed to even consider the cell, i.e. you need at least two atoms in the cell in order to do a two-body collision
                pair_cand_counter = 0 #setting the current N of occured collision to zero at the start of the sequence
                N = np.size(c.stored_atoms) #number of atoms in the cell
                c.ens_avg_N_array.append(N) #appending the number of atoms in the cell from this time step into the time/ensemble averaged array, used for calaclauting the time/ensemble average
                N_ens_avg = np.average(c.ens_avg_N_array) #time/ensemble averaged number of atoms in the cell over the timesteps so far
                
                pair_cand = np.ceil(((N-1)*(N_ens_avg)*(self.cross_sect*self.vr_max)*dt)/(2*V_c)) #calculate the number of pair candidates in the cell
                print(f'no. of candidates ={pair_cand}')
                
                
                while pair_cand_counter != pair_cand: #statement ensuring that the required number of possible collisions are considered before moving to the next cell
                    atom_1 = random.choice(c.stored_atoms) #choosing the two atoms
                    c.stored_atoms.remove(atom_1) #removing from the array so that the same atom cannot be chosen, i.e. cannot collide with itself
                    atom_2 = random.choice(c.stored_atoms) 
                    c.stored_atoms.append(atom_1) #once atom_2 is chosen without it being the same atoms as atom_1, it must be reappended to the array to ensure that it is equally considered for the next round of particle collisons, if it hasn't already collided 
                    
                    c_r = atom_1.v-atom_2.v # determining the relative velocity for the pair of particles
                    c_r_mag = np.sqrt(c_r[0]**2 + c_r[1]**2 + c_r[2]**2)
                    R_f = random.random() #choosing a random value between 0 and 1 used for the acceptance-rejection method for collisions 
                    
                    if (c_r_mag)/(self.vr_max) < R_f: # acceptance-rejection method condition for collision to be accepted
                        
                        c.coll_atoms.append(atom_1) #appending the accepted atom_1 for collision to a seperate array for atoms that have/will collide
                        c.stored_atoms.remove(atom_1) #removing the collided atom_1 from the array of atoms in the cell to ensure they are not selected again for collisions
                        c.coll_atoms.append(atom_2) #appending the accepted atom_2 for collision to a seperate array for atoms that have/will collide
                        c.stored_atoms.remove(atom_2) #removing the collided atom_2 from the array of atoms in the cell to ensure they are not selected again for collisions
                        
                        c_m = (m1*atom_1.v + m2*atom_2.v)/(m1+m2) #determining the CoM velocity 
                        
                        #atom_1.v = env.rel_vel_final(atom_1, atom_2, c_m, c_r_mag)[0] #determing the post-collison velocity vectors and transforming back into the general frame of reference 
                        #atom_2.v = env.rel_vel_final(atom_1, atom_2, c_m, c_r_mag)[0] 
                        
                        atom_1.v, atom_2.v = env.rel_vel_final(atom_1, atom_2, c_m, c_r_mag)
                        
                        
                        
                        print('collision goes ahead!')
                        n_coll_cell +=1 #increasing the collision counter for this cell
                        
                    else:
                        print('failed to collide')
                        
                    pair_cand_counter += 1 #ensuring to increase the number for the current number of possible collisions in this time step considered by 1 each time a pair of randomly selected particles are cosidered for a collision
            #self.N_collisions.append(n_collisions)
            
            N_coll_dt.append(n_coll_cell) #append the number of collisions for this cell into the array for this timestep
            
        self.N_collisions.append(np.sum(N_coll_dt)) #calculate the sum of the collisions for this timestep (over all the cells) and append to the N_colllsions array
    
    
    
    def collisons_subcell(self, dt, Ncell_i):
        ''' function that will evalute the possibility and outcomes of two-body particle collisons in a cell over 
        all cells in the environment'''
        
        m1 = self.mass #mass of atom 1
        m2 = self.mass #mass of atom 2
        
        N_coll_dt = [] #create empty array to store the collisions for each cell for this timestep
        
        for w in self.space: #looping over all cells in the environment i.e. space
            V_c = (self.L/(Ncell_i*w.subcellno))**3 #volume of cell- have to recalc for each cell as the subcells are diff sizes
            for c in w.subcellspace: # for each cell, loop through the subcells in it
                n_coll_cell = 0 #resetting the collision counter for the new cell
                c.coll_atoms.clear() #clears the array for any particles stored in it from the last time-step's collided atoms
                
                print(np.size(c.stored_atoms))
                
                if np.size(c.stored_atoms) > 1: #minimum condition needed to even consider the cell, i.e. you need at least two atoms in the cell in order to do a two-body collision
                    pair_cand_counter = 0 #setting the current N of occured collision to zero at the start of the sequence
                    N = np.size(c.stored_atoms) #number of atoms in the cell
                    c.ens_avg_N_array.append(N) #appending the number of atoms in the cell from this time step into the time/ensemble averaged array, used for calaclauting the time/ensemble average
                    N_ens_avg = np.average(c.ens_avg_N_array) #time/ensemble averaged number of atoms in the cell over the timesteps so far
                    
                    pair_cand = np.ceil(((N-1)*(N_ens_avg)*(self.cross_sect*self.vr_max)*dt)/(2*V_c)) #calculate the number of pair candidates in the cell
                    print(f'no. of candidates ={pair_cand}')
                    
                    
                    while pair_cand_counter != pair_cand: #statement ensuring that the required number of possible collisions are considered before moving to the next cell
                        atom_1 = random.choice(c.stored_atoms) #choosing the two atoms
                        c.stored_atoms.remove(atom_1) #removing from the array so that the same atom cannot be chosen, i.e. cannot collide with itself
                        atom_2 = random.choice(c.stored_atoms) 
                        c.stored_atoms.append(atom_1) #once atom_2 is chosen without it being the same atoms as atom_1, it must be reappended to the array to ensure that it is equally considered for the next round of particle collisons, if it hasn't already collided 
                        
                        c_r = atom_1.v-atom_2.v # determining the relative velocity for the pair of particles
                        c_r_mag = np.sqrt(c_r[0]**2 + c_r[1]**2 + c_r[2]**2)
                        R_f = random.random() #choosing a random value between 0 and 1 used for the acceptance-rejection method for collisions 
                        
                        if (c_r_mag)/(self.vr_max) < R_f: # acceptance-rejection method condition for collision to be accepted
                            
                            c.coll_atoms.append(atom_1) #appending the accepted atom_1 for collision to a seperate array for atoms that have/will collide
                            c.stored_atoms.remove(atom_1) #removing the collided atom_1 from the array of atoms in the cell to ensure they are not selected again for collisions
                            c.coll_atoms.append(atom_2) #appending the accepted atom_2 for collision to a seperate array for atoms that have/will collide
                            c.stored_atoms.remove(atom_2) #removing the collided atom_2 from the array of atoms in the cell to ensure they are not selected again for collisions
                            
                            c_m = (m1*atom_1.v + m2*atom_2.v)/(m1+m2) #determining the CoM velocity 
                             
                            atom_1.v, atom_2.v = env.rel_vel_final(atom_1, atom_2, c_m, c_r_mag)
                            
                            print('collision goes ahead!')
                            n_coll_cell +=1 #increasing the collision counter for this cell
                            
                        else:
                            print('failed to collide')
                            
                        pair_cand_counter += 1 #ensuring to increase the number for the current number of possible collisions in this time step considered by 1 each time a pair of randomly selected particles are cosidered for a collision
                #self.N_collisions.append(n_collisions)
                
                N_coll_dt.append(n_coll_cell) #append the number of collisions for this cell into the array for this timestep
                
        self.N_collisions.append(np.sum(N_coll_dt)) #calculate the sum of the collisions for this timestep (over all the cells) and append to the N_colllsions array
                        
                        
            
                
                
################################################################################################################################################################################################################################################################################# 



           
            
    def histogram(self, N, Nt, dt, label):
        '''Function to plot a histogram with associated Gaussian distribution.
        Displays the mean and standard deviation of this Gaussian distribution on the graph.'''
        
        to_plot = [] #empty array
        for i in self.particles:
            to_plot.append(Particle.mag_vr(i)) #select which values to plot, and append to the array
        fig = plt.figure(figsize=(8,3))
        gs = gridspec.GridSpec(1,1)
        ax2 = fig.add_subplot(gs[0])
        #set title with the simulation parameters to keep track of the data
        ax2.set_title(f'T={self.Ti}, omega={(self.omega_x, self.omega_y, self.omega_z)}, N={np.size(self.particles)}, Nt={Nt}, dt= {dt}')
        result1 = ax2.hist(to_plot, bins=100) #plotting the histogram
        mean1 = np.mean(to_plot)
        sigma1 = np.sqrt(np.var(to_plot))
        x1 = np.linspace(min(to_plot), max(to_plot), 100)
        dx1 = result1[1][1] - result1[1][0]
        scale1 = len(to_plot)*dx1 #scaling the curves to match the histogram
        
        #plotting a Gaussian of the results, for single direction 
        #ax2.plot(x1, scipy.stats.norm.pdf(x1, mean1, sigma1)*scale1, color='g')
        #ax2.set_xlabel('x position of particles')
        #ax2.set_ylim(0,300)
        #ax2.set_xlim(-0.001,0.001)
        #ax2.text(0.0005, 250, 'mean = {:.3g}'.format(mean1))
        #ax2.text(0.0005, 200, 'sd = {:.3g}'.format(sigma1))
        #ax2.text(-0.0009, max(result1[0])*0.95, '{}'.format(label))
 
        #fitting a maxwell distribution to the data, for radial coords
        params = scipy.stats.maxwell.fit(to_plot, floc=0) 
        ax2.plot(x1, scipy.stats.maxwell.pdf(x1, *params)*scale1)
        ax2.set_xlabel('Speed of particles')
        ax2.text(max(result1[1])*0.85, max(result1[0])*0.95, 'mean = {:.3g}'.format(mean1))
        ax2.text(max(result1[1])*0.85, max(result1[0])*0.75, 'sd = {:.3g}'.format(sigma1))
        ax2.text(max(result1[1])*0.85, max(result1[0])*0.55, '{}'.format(label))
        ax2.set_xlim(0,0.06)
        ax2.set_ylim(0,55)
        
        
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
        
        
#############################################################################################################
        
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
    
    
    def time_evolve_sorting_and_evap(self,dt,Nt,N, Ncell_i, therm_time, rate_of_evap):
            
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
            checking_cube = False #set to true to plot xy, xz,yz slices with the cell chosen marked- particles in that cell are marked in pink
            c=0 # the index of the cell we're checking
            
            
            sigma = [] #creating an empty array to contain the standard deviation values
            sigma_avg = []
            num_atoms = []
            r_vals = [] #defining an empty array for the magnitude radius values
            
            t = float(0) #ensures the value for t is treated as a float to avoid compatibility conditons
            
            for i in self.particles: #looping through all the particles in the system
                r_vals.append(Particle.mag_r(i))  #to append the magnitude radius to the r_vals array
    
            cut_off = max(r_vals) #setting the cut-off point for minimum radius from the paricles as an arbitrary max value of the r_vals
            evap_times = np.linspace(0, Nt, int(Nt/(therm_time)), endpoint=False)#defining the specific timesteps we want to do evaporative cooling at

            print(f'Starting Cut-Off is {cut_off}')
            evap_atoms = []
            
            if create_histograms:
                self.histogram(N, Nt, dt, 't = 0') #create a histogram of the initial values
                
            self.sort_atoms(Ncell_i) #initially sorting the atoms into the appropriate arrays associated with cells
            
            while t != Nt: #loop through the timesteps
                print(f'~~~~~~~~~~~~~~~~~~~Start of New Timestep {t} of {Nt}~~~~~~~~~~~~~~~~~')
                sigma_variable = [] #creating an empty array to store the r/v values from which to calc an sd
                sigma_variablex = [] 
                sigma_variabley = []
                sigma_variablez = []
                
                
                if create_graph3d:# and t%10 == 0: #can choose to only print graphs every few timesteps
                    fig = plt.figure(figsize=(10,10))
                    gs = gridspec.GridSpec(1,1)
                    ax1 = fig.add_subplot(gs[0], projection='3d') # making the plot 3D capable 
                    ax1.set_ylim(-self.L/2,self.L/2) 
                    ax1.set_xlim(-self.L/2,self.L/2)
                    ax1.set_zlim(-self.L/2,self.L/2)
                    ax1.text(self.L/3,-self.L/3, self.L/3, f'{t+1}') # stating on the graph the time-step of the simulation graph
                if create_graph2d and t%10 == 0:
                    fig = plt.figure(figsize=(3,3))
                    gs = gridspec.GridSpec(1,1)
                    ax1 = fig.add_subplot(gs[0])
                    ax1.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
                    ax1.set_xlim(-self.L/2,self.L/2)
                    ax1.set_xlabel('x')
                    ax1.set_ylabel('y')
                    ax1.text(self.L/3,-self.L/3, f'{t+1}') #print the timestep number on the graph
                    curr_num_atoms = np.size(self.particles)
                    ax1.text(0,-1.5*(self.L/2), f'current total atoms = {curr_num_atoms}')
                
                    if checking_cube:
                        #plots a square showing the boundaries of the chosen cell. The atoms in that cell's array are indicated in pink, cell centre marked with black cross.
                        ax1.scatter(self.space[c].centre_x, self.space[c].centre_y,marker='x',c='k', s=24)
                        ax1.add_patch(Rectangle((self.space[c].centre_x - self.L/(2*Ncell_i), self.space[c].centre_y - self.L/(2*Ncell_i)), self.L/(Ncell_i), self.L/(Ncell_i), fill=False, ls='--'))
                    
                    
                        fig1 = plt.figure(figsize=(3,3))
                        gs1 = gridspec.GridSpec(1,1)
                        ax2 = fig1.add_subplot(gs1[0])
                        ax2.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
                        ax2.set_xlim(-self.L/2,self.L/2)
                        ax2.set_xlabel('x')
                        ax2.set_ylabel('z')
                        ax2.text(self.L/3,-self.L/3, f'{t+1}') #print the timestep number on the graph
                        ax2.scatter(self.space[c].centre_x, self.space[c].centre_z, marker='x', c='k', s=24)
                        ax2.add_patch(Rectangle((self.space[c].centre_x - self.L/(2*Ncell_i), self.space[c].centre_z - self.L/(2*Ncell_i)), self.L/(Ncell_i), self.L/(Ncell_i), fill=False, ls='--'))
                    
                        fig3 = plt.figure(figsize=(3,3))
                        gs3 = gridspec.GridSpec(1,1)
                        ax3 = fig3.add_subplot(gs3[0])
                        ax3.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
                        ax3.set_xlim(-self.L/2,self.L/2)
                        ax3.set_xlabel('y')
                        ax3.set_ylabel('z')
                        ax3.text(self.L/3,-self.L/3, f'{t+1}') #print the timestep number on the graph
                        ax3.scatter(self.space[c].centre_y, self.space[c].centre_z, marker='x', c='k', s=24)
                        ax3.add_patch(Rectangle((self.space[c].centre_y - self.L/(2*Ncell_i), self.space[c].centre_z - self.L/(2*Ncell_i)), self.L/(Ncell_i), self.L/(Ncell_i), fill=False, ls='--'))
                        
                
                for i in self.particles: #for each timestep, loop through the particles in the array
                    #plot the positions if needed
                    if create_graph3d and t%10 == 0:
                        ax1.scatter3D(i.x, i.y, i.z, c='b', s=2)
                        ax1.view_init(40, 0)
                    if create_graph2d and t%10 == 0:
                        ax1.scatter(i.x, i.y, c='b', s=2)
                        if checking_cube: #this is the section of the code that marks the atoms in the cell's array pink
                            for n in self.space[c].stored_atoms:
                                ax1.scatter(n.x, n.y, c='m', s=6)
                                
                            ax2.scatter(i.x, i.z, c='g', s=2)
                            for n in self.space[c].stored_atoms:
                                ax2.scatter(n.x, n.z, c='m', s=6)
                                
                            ax3.scatter(i.y, i.z, c='y', s=2)
                            for n in self.space[c].stored_atoms:
                                ax3.scatter(n.y, n.z, c='m', s=6)
                    
                    sigma_variable.append(Particle.mag_r(i))
                    sigma_variablex.append(i.x) #append the r/v values into the empty array
                    sigma_variabley.append(i.y)
                    sigma_variablez.append(i.z)
                    
            
                    #alter the positions of the particles according to their velocities
                    Particle.drift(i, dt)
                    #alter the velocities of the particles, account for the trapping potential
                    Particle.potential_v_change(i, self.omega_x, self.omega_y, self.omega_z, dt)
                    
                #self.sort_atoms_again(Ncell_i) #sorting the atoms- only keeping atoms that are still in the cell at the end of the timestep
                #print(f'atoms sorted t={t}!')               
                #env.collisons(dt, Ncell_i) #dt, Ncell
                self.create_subcells(Ncell_i, t)
                print(f'sorted atoms t={t}!')
                self.collisons_subcell(dt, Ncell_i)
                
                if np.any(t == evap_times): #if statement, a conditon to see t is equal to any of the evap_times
                    print(f'time to evap t={t}!')
                    for i in self.particles: #looping through all partilces
                        if Particle.mag_r(i) >= cut_off: #if statement, a condition to see if magnitude radius of the partilce is greater than or equal to the cut_off radius
                            i.x = 0.004 #if the statment is satisified we set the position of the particle to an arbirtary position, to mimic the effects of evaporative cooling
                            i.y = 0.004 #""
                            i.z = 0.004 #""
                            atom = Particle(i.x, i.y, i.z, i.vx, i.vy, i.vz)
                            evap_atoms.append(atom)
                            self.particles.remove(i)
                    N_evap_atoms = np.size(evap_atoms)
                    print(f'{N_evap_atoms} got evaporated!')
                       
                        
                    cut_off = cut_off*rate_of_evap #outside of this if-else statment, once all particles have been looped/checked against the statment, we redefine the cut-off radius at a rate defined in the funtion variable
                    #print(f'Setting new cut-off as {cut_off}')
                if create_histograms and t%2 == 0:
                        self.histogram(N, Nt, dt, f'{t+1}') #create a histogram of the final value
                if save_images and t%2 == 0: #save the graphs as they are created to the specified file
                        plt.savefig(r'C:\\Users\Bethan\Documents\evaporative cooling\test sim\bigrunapril54\timestep{t}.png'.format(t=t))
                
                
                print(f'~~~~~~~~~~~~~~~~~~~End of Timestep {t} of {Nt}~~~~~~~~~~~~~~~~~~~~~~~')
                t+=1 #ensures that time-steps progress by one every full cycle of checks 
                    
                sigmaxdt = np.sqrt(np.var(sigma_variablex)) 
                sigmaydt = np.sqrt(np.var(sigma_variabley))
                sigmazdt = np.sqrt(np.var(sigma_variablez))
                sigmadt_avg = (1/3)*(sigmaxdt + sigmaydt + sigmazdt)
                sigmadt = np.sqrt(np.var(sigma_variable)) #calculate the sd of the r/v values for this timestep
                sigma.append(sigmadt) #append it to the empty array created at the start
                sigma_avg.append(sigmadt_avg)
                num_atoms.append(np.size(self.particles))
                    
                if save_images: #save the graphs as they are created to the specified file
                     plt.savefig(r'C:\\Users\Bethan\Documents\evaporative cooling\test sim\bigrunapril4\timestep{t}.png'.format(t=t))
                
            print(f'total number of collsions: {np.sum(self.N_collisions)}')
            print(f'total number of particles remaining: {np.size(env.particles)}')
                
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
                        
#####################################################################################################################################        

    def time_evolve_sorting(self,dt,Nt,N, Ncell_i):
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
        checking_cube = False #set to true to plot xy, xz,yz slices with the cell chosen marked- particles in that cell are marked in pink
        c=0 # the index of the cell we're checking
        
        sigma = [] #creating an empty array to contain the standard deviation values
        
        if create_histograms:
            self.histogram(N, Nt, dt, 't = 0') #create a histogram of the initial values
            
        self.sort_atoms(Ncell_i) #initially sorting the atoms into the appropriate arrays associated with cells
        
        for t in range(Nt): #loop through the timesteps
            if create_graph3d:# and t%10 == 0: #can choose to only print graphs every few timesteps
                fig = plt.figure(figsize=(10,10))
                gs = gridspec.GridSpec(1,1)
                ax1 = fig.add_subplot(gs[0], projection='3d') # making the plot 3D capable 
                ax1.set_ylim(-self.L/2,self.L/2) 
                ax1.set_xlim(-self.L/2,self.L/2)
                ax1.set_zlim(-self.L/2,self.L/2)
                ax1.text(self.L/3,-self.L/3, self.L/3, f'{t+1}') # stating on the graph the time-step of the simulation graph
            if create_graph2d and t%10 == 0:
                fig = plt.figure(figsize=(3,3))
                gs = gridspec.GridSpec(1,1)
                ax1 = fig.add_subplot(gs[0])
                ax1.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
                ax1.set_xlim(-self.L/2,self.L/2)
                ax1.set_xlabel('x')
                ax1.set_ylabel('y')
                ax1.text(self.L/3,-self.L/3, f'{t+1}') #print the timestep number on the graph
            
                if checking_cube:
                    #plots a square showing the boundaries of the chosen cell. The atoms in that cell's array are indicated in pink, cell centre marked with black cross.
                    ax1.scatter(self.space[c].centre_x, self.space[c].centre_y,marker='x',c='k', s=24)
                    ax1.add_patch(Rectangle((self.space[c].centre_x - self.L/(2*Ncell_i), self.space[c].centre_y - self.L/(2*Ncell_i)), self.L/(Ncell_i), self.L/(Ncell_i), fill=False, ls='--'))
                
                
                    fig1 = plt.figure(figsize=(3,3))
                    gs1 = gridspec.GridSpec(1,1)
                    ax2 = fig1.add_subplot(gs1[0])
                    ax2.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
                    ax2.set_xlim(-self.L/2,self.L/2)
                    ax2.set_xlabel('x')
                    ax2.set_ylabel('z')
                    ax2.text(self.L/3,-self.L/3, f'{t+1}') #print the timestep number on the graph
                    ax2.scatter(self.space[c].centre_x, self.space[c].centre_z, marker='x', c='k', s=24)
                    ax2.add_patch(Rectangle((self.space[c].centre_x - self.L/(2*Ncell_i), self.space[c].centre_z - self.L/(2*Ncell_i)), self.L/(Ncell_i), self.L/(Ncell_i), fill=False, ls='--'))
                
                    fig3 = plt.figure(figsize=(3,3))
                    gs3 = gridspec.GridSpec(1,1)
                    ax3 = fig3.add_subplot(gs3[0])
                    ax3.set_ylim(-self.L/2,self.L/2) #sets the dimensions of the axes to be the same as the box
                    ax3.set_xlim(-self.L/2,self.L/2)
                    ax3.set_xlabel('y')
                    ax3.set_ylabel('z')
                    ax3.text(self.L/3,-self.L/3, f'{t+1}') #print the timestep number on the graph
                    ax3.scatter(self.space[c].centre_y, self.space[c].centre_z, marker='x', c='k', s=24)
                    ax3.add_patch(Rectangle((self.space[c].centre_y - self.L/(2*Ncell_i), self.space[c].centre_z - self.L/(2*Ncell_i)), self.L/(Ncell_i), self.L/(Ncell_i), fill=False, ls='--'))
                    
                    
                
                
            sigma_variable = [] #creating an empty array to store the r/v values from which to calc an sd
            
            for i in self.particles: #for each timestep, loop through the particles in the array
                #plot the positions if needed
                if create_graph3d and t%10 == 0:
                    ax1.scatter3D(i.x, i.y, i.z, c='b', s=2)
                    ax1.view_init(40, 0)
                if create_graph2d and t%10 == 0:
                    ax1.scatter(i.x, i.y, c='b', s=2)
                    if checking_cube: #this is the section of the code that marks the atoms in the cell's array pink
                        for n in self.space[c].stored_atoms:
                            ax1.scatter(n.x, n.y, c='m', s=6)
                            
                        ax2.scatter(i.x, i.z, c='g', s=2)
                        for n in self.space[c].stored_atoms:
                            ax2.scatter(n.x, n.z, c='m', s=6)
                            
                        ax3.scatter(i.y, i.z, c='y', s=2)
                        for n in self.space[c].stored_atoms:
                            ax3.scatter(n.y, n.z, c='m', s=6)
                    
                sigma_variable.append(i.x) #append the r/v values into the empty array
        
                #alter the positions of the particles according to their velocities
                Particle.drift(i, dt)
                #alter the velocities of the particles, account for the trapping potential
                Particle.potential_v_change(i, self.omega_x, self.omega_y, self.omega_z, dt)
                
        
            #self.sort_atoms_again(Ncell_i) #sorting the atoms- only keeping atoms that are still in the cell at the end of the timestep
            #self.collisons(dt, Ncell_i)
            print(f'about to sort atoms t={t}!')
            self.create_subcells(Ncell_i, t)
            print(f'sorted atoms t={t}!')
            self.collisons_subcell(dt, Ncell_i)
            
            sigmadt = np.sqrt(np.var(sigma_variable)) #calculate the sd of the r/v values for this timestep
            sigma.append(sigmadt) #append it to the empty array created at the start
              
            if create_histograms:
                self.histogram(N, Nt, dt, 't = {t}')
            if save_images: #save the graphs as they are created to the specified file
                 plt.savefig(r'C:\\Users\Bethan\Documents\evaporative cooling\test sim\trapthennot\timestep{t}.png'.format(t=t))
            
            
            
        print(f'total number of collsions: {np.sum(self.N_collisions)}') #prints the total number of collisions that occurred over the whole timeloop
            
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
            
N = 5000
Nt = 1
Ncell = 3
                        
env = Environment(0.002,10**-6,60,60,60) #L, T, omega
env.Create_Particle(N) #N
env.Create_Cell(Ncell) #Ncell
#env.Check_cells()
#env.time_evolve_sorting(0.0001, Nt, N, Ncell) #dt, Nt, N, Ncell
env.time_evolve_sorting_and_evap(0.0005, Nt, N, Ncell, 1, 0.95) #dt, Nt, N, Ncell, therm_time, rate of evap
#env.time_evolve_evap(0.0001, 200, 10000, 5, 0.97) #dt, Nt, N, evap_timestep, rate_of_evap
#env.time_evolve_TOF(0.0001, 1000, np.size(env.particles)) #dt, Nt, N
print(env.N_collisions)
