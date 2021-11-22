#!/usr/bin/env python
# coding: utf-8

# # Weight and combine the chunk velocities

# ## Background

# We have modelled all observations of our star and saved the results to file. To retrieve the RV timeseries of the star, we still need to combine the individual best-fit velocities of all 522 chunks for each observation. One might think of a simple mean or median over the chunk velocities within each observation to achieve this - but as we have seen [in the previous chapter](./observation.ipynb) we have quite a number of 'bad' chunks, with outlier velocities and huge errors. Also, some chunks contain more spectral information than others, deliver therefore more realistic results and should be weighted higher than poorly constraint outlier chunks.
# 
# In this script, we show how to run the combination routine. For a sound mathematical description, please check out our !!!paper!!!. But basically, weights are created for each chunk not only based on their 'velocity performance' within one observation, but also based on their performance throughout the whole timeseries. I.e. if the velocity timeseries of one chunk with index $i$ shows a highly significant scatter around the mean timeseries of all chunk velocities, that chunk receives lower weight IN ALL OBSERVATIONS!

# ## Run the code

# First of, as before, we set up the path structure and import the velocity combination routine:

# In[1]:


# Automatic reloading of imports
get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import sys
import os
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('/home/paul/pyodine/')  # Put in the pyodine path on your machine here!

import pyodine
import pyodine_combine_vels             # <- the velocity combination routine


# Again we load an object which contains all the important parameters for the velocity combination - this is the `Timeseries_Parameters` object:

# In[2]:


import utilities_song as utilities

Pars = utilities.timeseries_parameters.Timeseries_Parameters()


# Next, we need to specify the pathnames for the saved model results to use. We can also define pathnames for a directory where to save analysis plots, for a log-file with diagnosis information, and a text-file where to save the RV timeseries results in a human-readable format. Also, the combined model results (including best-fit results for all observations and chunks, and the RV timeseries along with additional information), represented as a `CombinedResults` object, can be saved to a HDF5-file whose pathname we also specify below.

# In[3]:


# Individual result files to use
parent_dir = '/home/paul/data_song2/data_res/sigdra_obs/'
res_files = [os.path.join(parent_dir, f, f+'_res1.h5') for f in os.listdir(parent_dir)]
res_files.sort()

# The directory for analysis plots
plot_dir = '/home/paul/data_song2/data_res/sigdra_combined'

# The output name of the CombinedResults object
comb_res_out = os.path.join(plot_dir, 'sigdra_tutorial_comb.h5')

# The output name of the velocity results text file
vels_out = os.path.join(plot_dir, 'sigdra_tutorial.vels')

# Log files
error_file = os.path.join(plot_dir, 'error.log')
info_file  = os.path.join(plot_dir, 'info.log')


# Finally, we run the velocity combination routine, which loads all the individual model results from file, computes the RV timeseries, and saves the `CombinedResults` object as well as some analysis plots and additional data as specified above. In the end the `CombinedResults` object is returned:

# In[4]:


Results = pyodine_combine_vels.combine_velocity_results(
    Pars, res_files=res_files, plot_dir=plot_dir, comb_res_out=comb_res_out, 
    vels_out=vels_out, error_log=error_file, info_log=info_file)


# Great, and again some print output to explain:
# 
# - First up, after the individual model results are loaded, barycentric velocities are computed for each observation. This is done with the Python package [barycorrpy](https://github.com/shbhuk/barycorrpy), and using the stellar coordinates supplied by the first observation in the timeseries (they have been taken from the respective Fits-header and saved in the individual results file).
# 
# - Then, the weighting parameters used in the velocity combination (and supplied through the `Timeseries_Parameters` object) are printed. Also the code reports about the number of observations and chunks that it found in the results.
# 
# - Next, the (robust) mean and standard deviation (std) of the offsets of chunk timeseries from observation means are printed. Below, the mean and std of the velocity scatter within each chunk timeseries is reported ('Chunk sigmas'). The sigmas are used to weight the chunk velocities when computing the RV timeseries.
# 
# - Finally, two 'quality factors' for the velocity weighting are reported. They should always be very similar, and as they basically represent the RV precision achieved in this timeseries, you would want them to be as small as possible.

# ## Inspect the results

# Now finally, after all this modelling, we can plot our RV timeseries. The barycentric dates, corrected RVs and uncertainties are available as attributes of the `CombinedResults` object called 'Results':

# In[7]:


# The barycentric date, RV with barycentric correction, and RV uncertainty
bary_date = Results.bary_date
rv_bc     = Results.rv_bc
rv_err    = Results.rv_err

star_name = Results.info['star_name']     # the star name

# And plot
fig = plt.figure(figsize=(10,6))
plt.errorbar(bary_date, rv_bc, yerr=rv_err, fmt='o', alpha=0.7,
            label='{:.2f} +- {:.2f} m/s'.format(np.mean(rv_bc), np.std(rv_bc)))
plt.legend()
plt.xlabel('Time [JD]')
plt.ylabel('Velocity [m/s]')
plt.title('{}: RV timeseries'.format(star_name))


# In addition to the RV timeseries, the 'Results' object contains meta-information from the velocity weighting algorithm, that can help us to evaluate the trustworthyness of the results. For example, we can look at the scatter of individual chunk velocities (i.e. the robust standard deviation) for each observation:

# In[8]:


c2c_scatter = Results.c2c_scatter

fig = plt.figure(figsize=(10,6))
plt.plot(Results.bary_date, c2c_scatter, 'o', alpha=0.7,
        label='{:.2f} +- {:.2f} m/s'.format(np.mean(c2c_scatter), np.std(c2c_scatter)))
plt.legend()
plt.xlabel('Time [JD]')
plt.ylabel('Velocity [m/s]')
plt.title('{}: Chunk velocity scatter in observations'.format(star_name))


# Even more: We can plot the chunk sigmas (velocity scatter within each chunk timeseries) over the mean counts of each chunk timeseries. We expect a downward correlation, where chunks with higher counts have smaller sigmas (as they should deliver a better estimate of the true RV in each observation):

# In[15]:


# Extract the median counts (medcnts) of each chunk in each observation, and average them over the observations
medcnts_mean = np.mean(Results.medcnts, axis=0)
# Extract the chunk sigmas
chunk_sigma  = Results.auxiliary['chunk_sigma']

fig = plt.figure(figsize=(10,6))
plt.plot(medcnts_mean, chunk_sigma, '.', alpha=0.7)
plt.xlabel('Mean of chunk fluxes [ADU]')
plt.ylabel('Chunk sigma [m/s]')
plt.ylim(0,200)
plt.title('{}: Chunk sigmas over mean chunk fluxes'.format(star_name))


# The correlation is just there, but so are a lot of outliers - most of which are chunks at redder wavelengths, where tellurics and the fringing effect negatively impact the model.
# 
# In addition to the RV timeseries, the velocity weighting algorithm computes the so-called Chromatic Index for each observation - a measure of the slope of Doppler velocity over wavelength. It can be used as an activity indicator, as a varying Chromatic Index can be a sign of stellar spots leading to a RV variation (see e.g. [Reiners +2010](https://ui.adsabs.harvard.edu/abs/2010ApJ...710..432R/abstract) for a physical description, and [Zechmeister +2017](https://ui.adsabs.harvard.edu/abs/2018A%26A...609A..12Z/abstract) for a code implementation and exemplary results).
# 
# *Note*: This is still very experimental - we do not know how significant the Chromatic Index results from I2 cell codes are, and certainly need to test this more!

# In[16]:


crx     = Results.crx
crx_err = Results.crx_err

fig = plt.figure(figsize=(10,6))
plt.errorbar(bary_date, crx, yerr=crx_err, fmt='o', alpha=0.7,
            label='{:.2f} +- {:.2f} (m/s)/Np'.format(np.mean(crx), np.std(crx)))
plt.legend()
plt.xlabel('Time [JD]')
plt.ylabel('Chromatic index [(m/s)/Np]')
plt.title('{}: Chromatix index variation'.format(star_name))


# Looks quite stable! And to give you an idea of what is measured by the Chromatic Index, let's plot the individual chunk velocities for one observation and overplot the fitted Chromatic Index model, which looks like:
# 
# $$
#     v_\mathrm{fit} (\lambda) = RV + \beta \ln \frac{\lambda_\mathrm{RV}}{\lambda} \quad ,
# $$
# 
# where $RV$ is the weighted RV of the observation, $\beta$ is the Chromatic Index, and $\lambda_\mathrm{RV}$ is the effective wavelength of the weighted RV.

# In[30]:


obs_ind = 10     # observation index

# The chunk wavelength intercepts of the observation
wave_int     = Results.params['wave_intercept'][obs_ind]
# The chunk velocities of the observation
velocity     = Results.params['velocity'][obs_ind]
# The weighted RV of the observation
rv_obs       = Results.rv[obs_ind]
# The CRX model effective wavelength of the RV (and error)
RV_wave      = Results.RV_wave[obs_ind]
RV_wave_err  = Results.RV_wave_err[obs_ind]
# The CRX model chromatic index (and error)
crx_obs      = Results.crx[obs_ind]
crx_obs_err  = Results.crx_err[obs_ind]

# Compute the velocities from the CRX model, along with upper and lower uncertainty
velocities_fit = pyodine.timeseries.combine_vels.velocity_from_chromatic_index(
    wave_int, rv_obs, RV_wave, crx_obs)
velocities_fit_lower = pyodine.timeseries.combine_vels.velocity_from_chromatic_index(
    wave_int, rv_obs-rv_err[obs_ind], RV_wave-RV_wave_err, crx_obs-crx_obs_err)
velocities_fit_upper = pyodine.timeseries.combine_vels.velocity_from_chromatic_index(
    wave_int, rv_obs+rv_err[obs_ind], RV_wave+RV_wave_err, crx_obs+crx_obs_err)

fig = plt.figure(figsize=(10,6))
plt.plot(wave_int, velocity, '.', alpha=0.7, label='Chunk velocities')
plt.plot(wave_int, velocities_fit, alpha=0.7, label='CRX model\n' + 
         r'$\beta={:.1f}\pm{:.1f}$ (m/s)/Np'.format(crx_obs, crx_obs_err))
plt.plot(wave_int, velocities_fit_lower, 'k', alpha=0.3)
plt.plot(wave_int, velocities_fit_upper, 'k', alpha=0.3)
plt.legend()
plt.xlabel(r'Wavelength [$\AA$]')
plt.ylabel('Chunk velocities [m/s]')
plt.ylim(400., 1600.)
plt.title('{}, observation {}: Chunk velocities and CRX model'.format(star_name, obs_ind))


# In[ ]:




