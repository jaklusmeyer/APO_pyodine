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

# In[1]:


# Automatic reloading of imports
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import sys
import os
import matplotlib.pyplot as plt

sys.path.append('/home/paul/pyodine/')

import pyodine
import pyodine_create_templates


# Also, we import the utilities module for the instrument we are working with - in this case for SONG. It contains all the instrument-specific code, particularly routines to import the spectra from fits along with useful information from the fits headers. Upon import we call the SONG-specific ``utilities_song`` simply ``utilities``:

# In[2]:


import utilities_song as utilities


# Additionally, a parameter input object of type ``Template_Parameters`` is needed. That one contains all the instruments-specific parameters such as oversampling used, how to chunk up the spectrum, and also the exact definition of the workflow for the analysis routine. By default, there should be a module called ``pyodine_parameters`` within each ``utilities`` module, but you could also create a different one and import it if you wish to experiment with changing parameters. However, we will stick with the well-tested default parameters here:

# In[3]:


Pars = utilities.pyodine_parameters.Template_Parameters()


# Finally, we need to define the pathnames to the observations of the O-star with I2, and to the observations of the star without I2. We simply collect them from the respective directories of the tutorial data. Also, we define the output pathname for the deconvolved template, the pathname for a directory where to save analysis plots, and the output pathnames of the collected results. In our analysis, we model the O-star data in two runs (first with a Single-Gaussian LSF to establish good parameter guesses, and then with the final Multi-Gaussian LSF), and we save the results from both runs - the first as '.h5' (HDF5 format), the second as '.pkl' (through the **dill** Python package).

# In[4]:


# O-star observations to use for the modelling
ostar_dir = '/home/paul/data_song/data_ext/sigdra_template/2018-05-16/obs_ostar'
ostar_files = [os.path.join(ostar_dir, f) for f in os.listdir(ostar_dir)]
ostar_files.sort()

# Stellar observations to use for the deconvolution
temp_dir = '/home/paul/data_song/data_ext/sigdra_template/2018-05-16/obs_temp'
temp_files = [os.path.join(temp_dir, f) for f in os.listdir(temp_dir)]
temp_files.sort()

# Output pathname for the template
temp_outname = '/home/paul/data_song2/templates/temp_sigdra_2018-05-16.h5'

# Output directory for plots and pathnames for modelling results
plot_dir = '/home/paul/data_song2/data_res/sigdra/'
res_files = ['/home/paul/data_song2/data_res/sigdra/sigdra_res0.h5',
             '/home/paul/data_song2/data_res/sigdra/sigdra_res1.pkl']


# And now, we can kick off the template creation:

# In[5]:


pyodine_create_templates.create_template(utilities, Pars, ostar_files, temp_files, temp_outname, 
                    plot_dir=plot_dir, res_files=res_files)


# Great, everything went fine!

# In[ ]:





# Open the saved results file and inspect it

# In[6]:


fit_results_1 = pyodine.fitters.results_io.load_results(res_files[1], filetype='dill')

chunks = pyodine.components.ChunkArray()
for r in fit_results_1:
    chunks.append(r.chunk)


# In[7]:


residuals = pyodine.plot_lib.plot_residual_hist(fit_results_1, title='Residuals histogram', show_plot=True)


# In[8]:


pyodine.plot_lib.plot_chunkmodel(fit_results_1, chunks, 270, template=False, show_plot=True)


# In[9]:


lsf_model = fit_results_1[0].model.lsf_model
lsf_x = lsf_model.generate_x(6, conv_width=6.)

lsfs = []
for i in range(len(fit_results_1)):
    lsf_pars = fit_results_1[i].params.filter('lsf')
    lsfs.append(lsf_model.eval(lsf_x, lsf_pars))

pyodine.plot_lib.plot_lsfs_grid(lsfs, chunks, x_lsf=lsf_x, x_nr=3, y_nr=3, alpha=0.7, xlim=(-4,4), show_plot=True)


# In[10]:


wave_slopes_model = [r.params['wave_slope'] for r in fit_results_1]
wave_slopes_data = [(ch.wave[-1]-ch.wave[0])/len(ch) for ch in chunks]

pyodine.plot_lib.plot_chunk_scatter(scatter=[wave_slopes_model,wave_slopes_data], 
                                    scatter_fmt='.', scatter_label=['model', 'data'], 
                                    ylabel=r'wave_slope [$\AA$/pix]', show_plot=True)


# In[ ]:





# In[ ]:





# In[19]:


chunks, fit_results_0 = pyodine.fitters.results_io.restore_results_object(utilities, res_files[0])


# In[20]:


lsf_model = fit_results_0[0].model.lsf_model
lsf_x = lsf_model.generate_x(6, conv_width=6.)

lsfs = []
for i in range(len(fit_results_0)):
    lsf_pars = fit_results_0[i].params.filter('lsf')
    lsfs.append(lsf_model.eval(lsf_x, lsf_pars))

pyodine.plot_lib.plot_lsfs_grid(lsfs, chunks, x_lsf=lsf_x, x_nr=3, y_nr=3, alpha=0.7, xlim=(-4,4), show_plot=True)


# In[ ]:




