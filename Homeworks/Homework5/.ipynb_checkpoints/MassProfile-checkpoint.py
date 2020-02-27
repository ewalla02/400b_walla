{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "# create class MassProfile\n",
    "class MassProfile:\n",
    "    # Class to define the mass profile of a galaxy at a snapshot time\n",
    "    \n",
    "    def __init__(self, galaxy, snap):\n",
    "        \n",
    "        # inputs:\n",
    "        #       galaxy: a string with galaxy name, such as \"MW\" or \"M31\"\n",
    "        #       snap:  the snapshot number, such as 0, 1, etc\n",
    "        \n",
    "        # Initialize instance of this class with the following properties:\n",
    "        \n",
    "        # store the name of the galaxy\n",
    "        self.gname = galaxy\n",
    "        \n",
    "        # add a string of the filenumber to the value \"000\"\n",
    "        ilbl = \"000\"+str(snap)\n",
    "        \n",
    "        # remove all but the last 3 digits\n",
    "        ilbl = ilbl[-3:]\n",
    "        self.filename = \"%s_\"%(galaxy)+ilbl+\".txt\"\n",
    "      \n",
    "        # read in the data of our desired file\n",
    "        self.time, self.total, self.data = Read(self.filename) \n",
    "    \n",
    "        # read in the data for the x, y, z positions and the mass\n",
    "        self.m = self.data[\"m\"]\n",
    "        self.x = self.data[\"x\"]*u.kpc\n",
    "        self.y = self.data[\"y\"]*u.kpc\n",
    "        self.z = self.data[\"z\"]*u.kpc\n",
    "        \n",
    "    def MassEnclosed(self, component, r):\n",
    "        # function that computes the mass enclosed\n",
    "        #      within a given radius of the COM position\n",
    "        #      for a specified galaxy & a specified component of galaxy\n",
    "        \n",
    "        \n",
    "        # inputs:\n",
    "        #      component: the component you want the mass of\n",
    "        #           either 1 (Halo), 2 (Disk) or 3 (Bulge)\n",
    "        #      r: the radii (kpc) in which you want to calculate the mass\n",
    "        #           input r as an array\n",
    "        \n",
    "        # returns:\n",
    "        #      an array of the mass enclosed (Msun) within a given radius\n",
    "        \n",
    "        \n",
    "        # Determine COM Position\n",
    "        CoM = COM(self.filename,component).COM_P(0.00001)\n",
    "        # distance to the center of mass\n",
    "        COM_D = np.sqrt(CoM[0]**2 + CoM[1]**2 + CoM[2]**2)\n",
    "        \n",
    "        # create an index to extract only the component we are using\n",
    "        index = np.where(self.data[\"type\"]==component)\n",
    "        \n",
    "        # get mass, x, y z  of only the desired component\n",
    "        part_mass = self.m[index]\n",
    "        part_x = self.x[index].value - CoM[0]\n",
    "        part_y = self.y[index].value - CoM[1]\n",
    "        part_z = self.z[index].value - CoM[2]\n",
    "        \n",
    "        # calculate radii for our component\n",
    "        part_r = np.sqrt( (part_x**2)+(part_y**2)+(part_z**2) )\n",
    "        \n",
    "        # find Mass Enclosed within that radius of the COM\n",
    "        Mass = np.zeros_like(r)\n",
    "        \n",
    "        r_index = np.zeros_like(r)\n",
    "        \n",
    "        \n",
    "        i = 0\n",
    "        for i in range(len(r)):\n",
    "            # create index to get all points within given radius\n",
    "            #r_index[i] = np.where(part_r <= r[i])\n",
    "            Mass[i] = np.sum( part_mass[np.where( part_r < r[i])] )\n",
    "            \n",
    "        return Mass*1e10*u.Msun\n",
    "    \n",
    "    \n",
    "    \n",
    "    def MassEnclosedTotal(self, r):\n",
    "        # function calculates the total mass enclosed within the given radii\n",
    "        \n",
    "        # inputs:\n",
    "        #      r: a 1D array of radii (kpc)\n",
    "        \n",
    "        # returns:\n",
    "        #      a 1D array of the total masses enclosed within the radii\n",
    "        \n",
    "        \n",
    "        # make a caveat for M33, galaxy with no bulge\n",
    "        if self.gname == \"M33\":\n",
    "            iterate = 2 # M33 only has 2 components\n",
    "        else:\n",
    "            iterate = 3 # The other galaxies have 3 components\n",
    "        \n",
    "        # make an array to store the masses in\n",
    "        part_mass = np.zeros((iterate,len(r)))\n",
    "        \n",
    "        # get the masses of each component\n",
    "        i = 0\n",
    "        for i in range(iterate):\n",
    "            part_mass[i,:] = self.MassEnclosed(i+1,r)\n",
    "      \n",
    "        return np.sum(part_mass,axis=0)\n",
    "            \n",
    "    def HernquistMass(self, r, a, Mhalo):\n",
    "        # function computes the mass enclosed within a given radius\n",
    "        #     given the theoretical Hernquist profile\n",
    "        \n",
    "        \n",
    "        # inputs:\n",
    "        #     r: 1D array of radii (kpc) from Center of mass position\n",
    "        #     a: scale factor\n",
    "        #     Mhalo: the halo mass of the galaxy (Msun)\n",
    "        \n",
    "        # returns:\n",
    "        #      Returns halo mass (Msun)\n",
    "        \n",
    "        return (Mhalo*(r**2))/((a+r)**2)\n",
    "            \n",
    "        \n",
    "    def CircularVelocity(self, component, r):\n",
    "        # function that computes the circular velocity of the enclosed mass\n",
    "        # assumes spherical symmetry\n",
    "        # Vc = (GM/R)^1/2\n",
    "        \n",
    "        # inputs:\n",
    "        #     component: the type of particle you're working with\n",
    "        #           1 (halo), 2 (disk), 3 (bulge)\n",
    "        #     r: 1D array of radii (kpc)\n",
    "        \n",
    "        # returns:\n",
    "        #     1D array of circular speeds in units of km/s\n",
    "        \n",
    "        M = self.MassEnclosed(component, r)\n",
    "        \n",
    "        return np.sqrt(G*M/r)\n",
    "\n",
    "    def CircularVelocityTotal(self, r):\n",
    "        # function that calculates the circular velocity of the total enclosed mass\n",
    "        \n",
    "        # inputs: r: 1D array of radii (kpc)\n",
    "        \n",
    "        # returns: circular velocity of total enclosed mass (km/s)\n",
    "        \n",
    "        M = self.MassEnclosedTotal(r)\n",
    "        \n",
    "        return np.sqrt(G*M/r)\n",
    "    \n",
    "    def HernquistVCirc(self,r, a, Mhalo):\n",
    "        # computes the circular speed of the Hernquist mass\n",
    "        \n",
    "        # inputs:\n",
    "        #     r: 1D array of radii (kpc) from Center of mass position\n",
    "        #     a: scale factor\n",
    "        #     Mhalo: the halo mass of the galaxy (Msun)\n",
    "        \n",
    "        # returns: the circular speed of the Hernquist Mass (km/s)\n",
    "    \n",
    "        M = self.HernquistMass(r,a,Mhalo)\n",
    "        return np.round((np.sqrt(G*M/r)),2)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
