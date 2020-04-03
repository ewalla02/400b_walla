#!/usr/bin/env python
# coding: utf-8

# In[1]:


# make necessary imports
import numpy as np
import astropy.units as u
from Readfile import Read
import pandas as pd


# In[ ]:


def ComponentMass(filename, Ptype):
    # inputs:
    #      filename: the name of the txt file you wish to read
    #      Ptype: the type of particle you wish to find the mass of.
    #             options for Ptype = 1, 2, or 3.
    #             1 = Halo matter, 2 = Disk material, 3 = Bulge stuff
    
    # returns:
    #      The total mass of your desired galactic component, in 1e12Msun
    
    # read in the file
    time, number, data = Read(filename)
    
    # select only that data of your selected material type
    select_data = data[np.where(data["type"]==Ptype)]
    
    # extract mass data from your selected data
    mass_array = select_data["m"]
    
    # find total mass of the selected material
    total_mass = np.around(np.sum(mass_array)*u.Msun, 3)*1e10
    
    return np.around(total_mass/1e12, 3)
    
    

