#!/usr/bin/env python
# coding: utf-8

# In[1]:


# inputs
import numpy as np
import astropy.units as u
import astropy.table as tbl
import matplotlib.pyplot as plt
from Readfile import Read
from CenterOfMass import CenterOfMass as COM
from astropy.constants import G
# convert G to correct units
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

# create class MassProfile
class MassProfile:
    # Class to define the mass profile of a galaxy at a snapshot time
    
    def __init__(self, galaxy, snap):
        
        # inputs:
        #       galaxy: a string with galaxy name, such as "MW" or "M31"
        #       snap:  the snapshot number, such as 0, 1, etc
        
        # Initialize instance of this class with the following properties:
        
        # store the name of the galaxy
        self.gname = galaxy
        
        # add a string of the filenumber to the value "000"
        ilbl = "000"+str(snap)
        
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "/home/astr400b/VLowRes/%s_"%(galaxy)+ilbl+".txt"
       
        # read in the data of our desired file
        self.time, self.total, self.data = Read(self.filename) 
    
        # read in the data for the x, y, z positions and the mass
        self.m = self.data["m"]
        self.x = self.data["x"]*u.kpc
        self.y = self.data["y"]*u.kpc
        self.z = self.data["z"]*u.kpc
        
    def MassEnclosed(self, component, r):
        # function that computes the mass enclosed
        #      within a given radius of the COM position
        #      for a specified galaxy & a specified component of galaxy
        
        
        # inputs:
        #      component: the component you want the mass of
        #           either 1 (Halo), 2 (Disk) or 3 (Bulge)
        #      r: the radii (kpc) in which you want to calculate the mass
        #           input r as an array
        
        # returns:
        #      an array of the mass enclosed (Msun) within a given radius
        
        
        # Determine COM Position
        CoM = COM(self.filename,component).COM_P(0.00001,VolDec = 4)
        # distance to the center of mass
        COM_D = np.sqrt(CoM[0]**2 + CoM[1]**2 + CoM[2]**2)
        
        # create an index to extract only the component we are using
        index = np.where(self.data["type"]==component)
        
        # get mass, x, y z  of only the desired component
        part_mass = self.m[index]
        part_x = self.x[index].value - CoM[0]
        part_y = self.y[index].value - CoM[1]
        part_z = self.z[index].value - CoM[2]
        
        # calculate radii for our component
        part_r = np.sqrt( (part_x**2)+(part_y**2)+(part_z**2) )
        
        # find Mass Enclosed within that radius of the COM
        Mass = np.zeros_like(r)
        
        r_index = np.zeros_like(r)
        
        
        i = 0
        for i in range(len(r)):
            # create index to get all points within given radius
            #r_index[i] = np.where(part_r <= r[i])
            Mass[i] = np.sum( part_mass[np.where( part_r < r[i])] )
            
        return Mass*1e10*u.Msun
    
    
    
    def MassEnclosedTotal(self, r):
        # function calculates the total mass enclosed within the given radii
        
        # inputs:
        #      r: a 1D array of radii (kpc)
        
        # returns:
        #      a 1D array of the total masses enclosed within the radii
        
        
        # make a caveat for M33, galaxy with no bulge
        if self.gname == "M33":
            iterate = 2 # M33 only has 2 components
        else:
            iterate = 3 # The other galaxies have 3 components
        
        # make an array to store the masses in
        part_mass = np.zeros((iterate,len(r)))
        
        # get the masses of each component
        i = 0
        for i in range(iterate):
            part_mass[i,:] = self.MassEnclosed(i+1,r)
      
        return np.sum(part_mass,axis=0)
            
    def HernquistMass(self, r, a, Mhalo):
        # function computes the mass enclosed within a given radius
        #     given the theoretical Hernquist profile
        
        
        # inputs:
        #     r: 1D array of radii (kpc) from Center of mass position
        #     a: scale factor
        #     Mhalo: the halo mass of the galaxy (Msun)
        
        # returns:
        #      Returns halo mass (Msun)
        
        return (Mhalo*(r**2))/((a+r)**2)
            
        
    def CircularVelocity(self, component, r):
        # function that computes the circular velocity of the enclosed mass
        # assumes spherical symmetry
        # Vc = (GM/R)^1/2
        
        # inputs:
        #     component: the type of particle you're working with
        #           1 (halo), 2 (disk), 3 (bulge)
        #     r: 1D array of radii (kpc)
        
        # returns:
        #     1D array of circular speeds in units of km/s
        
        M = self.MassEnclosed(component, r)
        
        return np.sqrt(G*M/r)

    def CircularVelocityTotal(self, r):
        # function that calculates the circular velocity of the total enclosed mass
        
        # inputs: r: 1D array of radii (kpc)
        
        # returns: circular velocity of total enclosed mass (km/s)
        
        M = self.MassEnclosedTotal(r)
        
        return np.sqrt(G*M/r)
    
    def HernquistVCirc(self,r, a, Mhalo):
        # computes the circular speed of the Hernquist mass
        
        # inputs:
        #     r: 1D array of radii (kpc) from Center of mass position
        #     a: scale factor
        #     Mhalo: the halo mass of the galaxy (Msun)
        
        # returns: the circular speed of the Hernquist Mass (km/s)
    
        M = self.HernquistMass(r,a,Mhalo)
        return np.round((np.sqrt(G*M/r)),2)


# In[ ]:




