#!/usr/bin/env python
# coding: utf-8

# # Model the stellar I2 observations

# Bla blub

# In[1]:


# Automatic reloading of imports
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import sys
import os
import matplotlib.pyplot as plt

sys.path.append('/home/paul/pyodine/')

import pyodine
import pyodine_model_observations


# Again, import the utilities, and the parameter input object:

# In[2]:


import utilities_song as utilities

Pars = utilities.pyodine_parameters.Parameters()


# Observation filenames, pathnames of deconvolved stellar templates, plot directories and results filenames:

# In[4]:


# Observations to model
obs_dir = '/home/paul/data_song/data_ext/sigdra_obs_tutorial/'
obs_files = [os.path.join(obs_dir, f) for f in os.listdir(obs_dir)]
obs_files.sort()

# Deconvolved stellar template to use
temp_file = '/home/paul/data_song2/templates/temp_sigdra_2018-05-16.h5'
temp_files = [temp_file] * len(obs_files)

# Output directories for plots and output pathnames for modelling results
plot_dir_parent = '/home/paul/data_song2/data_res/sigdra_obs/'
plot_dirs = []
res_files = []
for obs_file in obs_files:
    plot_dir_base = os.path.splitext(os.path.basename(obs_file))[0]
    plot_dirs += [os.path.join(plot_dir_parent, plot_dir_base)]
    
    res_files.append([os.path.join(plot_dirs[-1], 'sigdra_res0.h5'), 
                      os.path.join(plot_dirs[-1], 'sigdra_res1.h5')])


# Now model the observations:

# In[5]:


pyodine_model_observations.model_multi_observations(utilities, Pars, obs_files, temp_files, 
                                                    plot_dirs=plot_dirs, res_files=res_files)


# In[ ]:





# In[ ]:





# In[6]:


chunks, fit_results_1 = pyodine.fitters.results_io.restore_results_object(utilities, res_files[0][1])


# In[7]:


pyodine.plot_lib.plot_chunkmodel(fit_results_1, chunks, 270, template=True, show_plot=True)


# In[ ]:




