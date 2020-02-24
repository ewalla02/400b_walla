#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Homework 4
# Center of Mass Position and Velocity
# Emily Walla


# In[2]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl
from Readfile import Read


# In[3]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.m = self.data["m"][self.index]
        self.x = self.data["x"][self.index]
        self.y = self.data["y"][self.index]
        self.z = self.data["z"][self.index]
        self.vx = self.data["vx"][self.index]
        self.vy = self.data["vy"][self.index]
        self.vz = self.data["vz"][self.index]


    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)
    
    #Eq 1: COM: X_com = sum(position*m_i)/sum(m_i)

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Xcom = np.sum(a*m)/np.sum(m)
        # ycomponent Center of mass
        Ycom = np.sum(b*m)/np.sum(m)
        # zcomponent Center of mass
        Zcom = np.sum(c*m)/np.sum(m)
        
        return Xcom, Ycom, Zcom
    
    
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        RCOM = np.sqrt((XCOM**2)+(YCOM**2)+(ZCOM**2))


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt((xNew**2)+(yNew**2)+(zNew**2))

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW<=RCOM)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = np.sqrt((XCOM2**2)+(YCOM2**2)+(ZCOM2**2))
            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            #print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this.                                                                                              
            #print ("maxR", RMAX)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            xNew = self.x-XCOM2
            yNew = self.y-YCOM2
            zNew = self.z-ZCOM2
            RNEW = np.sqrt((xNew**2)+(yNew**2)+(zNew**2))

            # set the center of mass positions to the refined values                                                   
            XCOM = np.around((XCOM2),2)
            YCOM = np.around((YCOM2),2)
            ZCOM = np.around((ZCOM2),2)
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = np.array([XCOM, YCOM, ZCOM])

        # set the correct units usint astropy and round all values
        # and then return the COM positon vector
        # write your own code below
            return COMP
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x - COMX
        yV = self.y - COMY
        zV = self.z - COMZ
        RV = np.sqrt((xV**2)+(yV**2)+(zV**2))
       
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV<=RVMAX.value)

        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew =  self.m[indexV]
      
        # compute the center of mass velocity using those particles
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew,mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        # write your own code below
        COMV = [VXCOM, VYCOM, VZCOM]

        # return the COM vector                                                                                        
        return COMV

