#!/usr/bin/env python
# coding: utf-8

# # Prepare the environment

# In[9]:


from matplotlib import rcParams

rcParams["savefig.dpi"] = 120
rcParams["figure.dpi"] = 120
#rcParams["font.size"] = 20


# In the following tutorial, we will use **pyodine** to compute RVs from SONG spectra for the star !!!!XYZ!!!! This includes three basic steps:
# 
# - using I2-free observations of the star, plus B-star spectra *with* the I2 cell, to :doc:`create a deconvolved stellar template <template>`;
# 
# - with that template, :doc:`modelling the I2 observations <observation>` of the star to arrive at best-fit parameters (including, most importantly, the chunk velocities);
# 
# - :doc:`weighting and combining the chunk velocities <velocities>` to compute the RV timeseries of the observations.
# 
# In each of these steps, we will present how to use built-in methods of **pyodine** to visualize and analyze the results.
# 
# 
# But, first of all, we need to prepare everything.
# 
# One possibility to use **pyodine** from any location in your filesystem is by appending its root path to the Python system path, and then importing like this:

# In[2]:


import sys

sys.path.append('/home/paul/pyodine/')

import pyodine


# In the following, we will use some capabilities of **pyodine** to examine the I2 atlas and the utilities module.

# ## The I2 atlas

# -> if you have your own, bring it into the correct hdf5-format (find out more about `HDF5/h5py <http://www.h5py.org/>`_)
# 
# -> here, we just use the one by SONG
# 
# We use functions contained in `pyodine.lib.h5quick` to print information about the I2 atlas: `h5print` to show the structure of the HDF5 file/h5py object, and `h5data` to load the wavelength and flux vectors needed in the later modelling.

# In[12]:


from pyodine.lib.h5quick import h5print, h5data

I2_path = '/home/paul/data_lick/iodine_atlas/song_iodine_cell_01_65C.h5'

h = h5py.File(I2_path, 'r')
h5print(h)

flux_normalized = h5data(h['flux_normalized'])
wavelength_air  = h5data(h['wavelength_air'])

print('\nLength of the wavelength/flux vectors:', len(wavelength_air))
print('\nFull wavelength range (air): {} - {} Angstrom'.format(wavelength_air[0], wavelength_air[-1]))


# As you can see, the I2 atlas contains the wavelength range between roughly 4989 and 6498 Å.

# In[13]:


import matplotlib.pyplot as plt
import numpy as np

ind = np.where(np.logical_and(wavelength_air>5500., wavelength_air<5505.))

fig, axs = plt.subplots(2)
axs[0].set_title('I2 atlas')
axs[0].plot(wavelength_air, flux_normalized, alpha=0.7)
axs[1].plot(wavelength_air[ind], flux_normalized[ind], alpha=0.7)
axs[1].set_xlabel(r'Wavelength [$\AA$]')
axs[0].set_ylabel('Norm. flux')
axs[1].set_ylabel('Norm. flux')


# ## The utilities module

# -> prepare all instrument-specific code, such as I/O of spectra, modelling parameters etc.
# 
# -> again, we just use the existing ``utilities_song`` module

# In[ ]:




