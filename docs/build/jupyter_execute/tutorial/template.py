#!/usr/bin/env python
# coding: utf-8

# # Create a deconvolved stellar template

# Before modelling observations of the star to extract RVs, we need to create a high-S/N stellar template
# 
# 
# What you need:
# 
# - observation spectra of hot stars, which do not show any spectral features in the I2 wavelength range between ~ 5000 -- 6000 Å, obtained *with* the I2 cell in the light path; the recorded spectra then basically only consist of the I2 features, and can be used to determine the LSF of the instrument. 
# 
# - observation spectra of the star of interest *without* any I2 absorption features;

# First, we need to set up the path and import the required pyodine modules:

# In[ ]:


import sys
import os
import matplotlib.pyplot as plt

sys.path.append('/home/paul/pyodine/')

import pyodine
import pyodine_create_templates


# Also, we import the utilities module for the instrument we are working with - in this case for SONG. Upon import we call it just ``utilities``.

# In[ ]:


import utilities_song as utilities


# Additionally, a parameter input object of type ``Template_Parameters`` is needed. By default, there should be a module called ``pyodine_parameters`` within each ``utilities``, but you could also create a different one and import it if you wish you experiment with changing parameters. However, we will stick with the well-tested default parameters here:

# In[ ]:


Pars = utilities.pyodine_parameters.Template_Parameters()


# Load the pathnames to the observations of the O-star with I2, and to the observations of the star without I2:

# In[ ]:


# O-star observations to use for the modelling
ostar_files = []

# Stellar observations to use for the deconvolution
temp_files = []

# Output pathname for the template
temp_outname = '/home/paul/data_song/templates/blub.h5'

# Output directory for plots and pathnames for modelling results
plot_dir = ''
res_files = ['',
             '']


# And now, we can kick off the template creation:

# In[ ]:


pyodine_create_templates.create_template(utilities, Pars, ostar_files, temp_files, temp_outname, 
                    plot_dir=plot_dir, res_files=res_files)


# In[ ]:





# In[ ]:





# Open the saved results file and inspect it

# In[ ]:




