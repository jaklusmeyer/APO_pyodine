#!/usr/bin/env python
# coding: utf-8

# # Weight and combine the chunk velocities

# We have modelled all observations of our star and saved the results to file. To retrieve the RV timeseries of the star, we still need to combine the individual best-fit velocities of all 522 chunks for each observation. One might think of a simple mean or median over the chunk velocities to achieve this - but some chunks contain more spectral information than others, deliver therefore more realistic results and should be weighted higher than poorly constraint outlier chunks.
# 
# In this script, we show how to run the combination routine. For a sound mathematical description, please check out our !!!paper!!!.
# 
# First of, as before, we set up the path structure and import the velocity combination routine:

# In[6]:


# Automatic reloading of imports
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import sys
import os
import matplotlib.pyplot as plt

sys.path.append('/home/paul/pyodine/')

import pyodine
import pyodine_combine_vels         # <- the velocity combination routine


# Important parameters for the velocity combination routine are defined in a timeseries parameter object, which we also need to load:

# In[7]:


import utilities_song as utilities

Pars = utilities.timeseries_parameters.Timeseries_Parameters()


# Next, we need to specify the pathnames for the saved model results to use. We can also define pathnames for a directory where to save analysis plots, for a log-file with diagnosis information, and a text-file where to save the RV timeseries results in a human-readable format. Also, the combined model results (including best-fit results for all observations and chunks, and the RV timeseries along with additional information), represented as a :class:`CombinedResults` object, can be saved to a HDF5-file whose pathname we also specify below.

# In[8]:


# Individual result files to use
parent_dir = '/home/paul/data_song2/data_res/sigdra_obs/'
res_files = [os.path.join(parent_dir, f, 'sigdra_res1.h5') for f in os.listdir(parent_dir)]
res_files.sort()

# The directory for analysis plots
plot_dir = '/home/paul/data_song2/data_res/sigdra_obs/vel_combined'

# The diagnosis file to write analysis info into
diag_file = os.path.join(plot_dir, 'diag_file.log')

# The output name of the CombinedResults object
comb_res_out = os.path.join(plot_dir, 'sigdra_combined.h5')

# The output name of the velocity results text file
vels_out = os.path.join(plot_dir, 'sigdra.vels')


# Finally, we run the velocity combination routine, which loads all the individual model results from file, computes the RV timeseries, and saves the :class:`CombinedResults` object as well as some analysis plots and additional data as specified above. In the end the :class:`CombinedResults` object is returned:

# In[9]:


Results = pyodine_combine_vels.combine_velocity_results(
    Pars, res_files=res_files, diag_file=diag_file, plot_dir=plot_dir, 
    comb_res_out=comb_res_out, vels_out=vels_out)


# Great!

# In[25]:


fig = plt.figure(figsize=(10,6))
plt.errorbar(Results.bary_date, Results.rv_bc, yerr=Results.rv_err, fmt='o', alpha=0.7)
plt.xlabel('Time [JD]')
plt.ylabel('Velocity [m/s]')
plt.title('{}: RV timeseries'.format(Results.info['star_name']))


# In[ ]:





# In[23]:


fig = plt.figure(figsize=(10,6))
plt.errorbar(Results.bary_date, Results.crx, yerr=Results.crx_err, fmt='o', alpha=0.7)
plt.xlabel('Time [JD]')
plt.ylabel('Chromatic index [(m/s)/Np]')
plt.title('{}: Chromatix index variation'.format(Results.info['star_name']))


# In[ ]:





# In[27]:


fig = plt.figure(figsize=(10,6))
plt.plot(Results.bary_date, Results.c2c_scatter, 'o', alpha=0.7)
plt.xlabel('Time [JD]')
plt.ylabel('Velocity [m/s]')
plt.title('{}: Chunk velocity scatter in observations'.format(Results.info['star_name']))


# In[ ]:





# In[5]:


obs_ind = 10
velocities_fit = pyodine.timeseries.combine_vels.velocity_from_chromatic_index(
    Results.params['wave_intercept'][obs_ind], Results.rv[obs_ind], Results.RV_wave[obs_ind], Results.crx[obs_ind])
velocities_fit_lower = pyodine.timeseries.combine_vels.velocity_from_chromatic_index(
    Results.params['wave_intercept'][obs_ind], Results.rv[obs_ind]-Results.rv_err[obs_ind], 
    Results.RV_wave[obs_ind]-Results.RV_wave_err[obs_ind], Results.crx[obs_ind]-Results.crx_err[obs_ind])
velocities_fit_upper = pyodine.timeseries.combine_vels.velocity_from_chromatic_index(
    Results.params['wave_intercept'][obs_ind], Results.rv[obs_ind]+Results.rv_err[obs_ind], 
    Results.RV_wave[obs_ind]+Results.RV_wave_err[obs_ind], Results.crx[obs_ind]+Results.crx_err[obs_ind])

fig = plt.figure(figsize=(10,6))
plt.plot(Results.params['wave_intercept'][obs_ind], Results.params['velocity'][obs_ind], '.', alpha=0.7)
plt.plot(Results.params['wave_intercept'][obs_ind], velocities_fit, alpha=0.7)
plt.plot(Results.params['wave_intercept'][obs_ind], velocities_fit_lower, 'k', alpha=0.3)
plt.plot(Results.params['wave_intercept'][obs_ind], velocities_fit_upper, 'k', alpha=0.3)
plt.ylim(1500., 2400.)
plt.show()
plt.close()


# In[ ]:




