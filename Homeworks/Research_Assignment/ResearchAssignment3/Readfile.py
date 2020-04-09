#!/usr/bin/env python
# coding: utf-8

# # Readfile
# 
# Program has function that reads in a file and returns 

# In[1]:


# make necessary imports
import numpy as np
import astropy.units as u


# In[2]:


# define the Read function that takes the file name as input:

def Read(filename):
    # ope  the file
    file = open(filename,'r')
    
    # read the first line and store the time in units of Myr
    line1 = file.readline()
    label,value = line1.split()
    time = float(value)*u.Myr
    
    # read the second line and store the to total number of particles
    line2 = file.readline()
    label2, number = line2.split()
    
    
    #close the file
    file.close()
    
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)
    
    #return the desired outputs
    return value, number, data

