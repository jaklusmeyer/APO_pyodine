#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 15:45:44 2021

@author: pheeren
"""

import numpy as np
from .misc import robust_mean, robust_std, reweight
from ..lib.misc import printLog


class Combiner():
    
    
    def __init__(self, results):
        self.results = results
        
    
    def create_timeseries(self):
        pass
        
        
        

def combine_chunk_velocities(velocities, bvc, diag_file):
    
    """
    Here the actual velocity weighting starts.
    """
    
    alpha = 1.8
    beta = 8.0
    sigma = 2.0
    
    # How many observations and chunks?
    (nr_obs, nr_chunks) = velocities.shape
    
    # Set up the barycentric corrected chunk velocities
    vel_bc = np.transpose(np.transpose(velocities) + bvc)
    
    # For each observation, center the chunk velocities around 0 by
    # subtracting the robust mean of that observation
    vel_offset_corrected = np.zeros((nr_obs, nr_chunks))
    for i in range(nr_obs):
        vel_offset_corrected[i,:] = velocities[i,:] - robust_mean(velocities[i,150:350])
    
    # Calculate the mean offset of each individual chunk timeseries from
    # the observation means by taking the robust mean along the 
    # (offset-corrected) chunks
    offsets_chunk = robust_mean(vel_offset_corrected, axis=0)
    
    # Subtract the overall robust mean of all chunk velocities in all 
    # observations (after the chunk velocities were already offset-corrected,
    # so in fact this should barely make a difference)
    offsets_chunk -= robust_mean(vel_offset_corrected)
    
    # Print to file
    printLog(diag_file, '')
    printLog(diag_file, 'Chunk-to-chunk offsets from observation mean:')
    printLog(diag_file, 'Median: {:.2f} +- {:.2f}'.format(
            np.nanmedian(offsets_chunk), np.nanstd(offsets_chunk)))
    
    # Set up the sigma and deviation arrays
    sig = np.zeros(nr_chunks)
    dev = np.zeros((nr_obs, nr_chunks))
    
    # For each (offset-corrected) chunk timeseries, compute its robust scatter
    # -> this is the sigma array
    # Also center each chunk timeseries around 0 by subtracting its robust mean,
    # and then express the deviation of each chunk velocity in this timeseries
    # from 0 in terms of the sigma of this timeseries
    # -> this gives the dev array
    for j in range(nr_chunks):
        sig[j] = robust_std(vel_offset_corrected[:,j])
        dev[:,j] = (vel_offset_corrected[:,j] - robust_mean(vel_offset_corrected[:,j])) / sig[j]
    
    # In each observation, center the chunk deviations around 0 by subtracting
    # their robust mean
    for i in range(nr_obs):
        dev[i,:] = dev[i,:] - robust_mean(dev[i,:])
    
    # Finally, set very low and very high sigmas to a pre-defined value
    ind = np.where(np.logical_or(sig<=4., sig>=1000.))
    sig[ind] = 1000.
    
    # Print to file
    printLog(diag_file, '')
    printLog(diag_file, 'Chunk sigmas:')
    printLog(diag_file, 'Median: {:.2f} +- {:.2f}'.format(
            np.nanmedian(sig), np.nanstd(sig)))
    
    # Prepare the final results arrays
    rv = np.zeros(nr_obs)               # Weighted RV timeseries
    rv_bc = np.zeros(nr_obs)            # Weighted RV timeseries, BV-corrected
    c2c_scatter = np.zeros(nr_obs)      # The chunk-to-chunk velocity scatter in each observation
    ww = np.zeros((nr_obs,nr_chunks))   # 
    mdvel = np.zeros(nr_obs)            # The simple observation median after correcting for chunk timeseries offsets
    errvel = np.zeros(nr_obs)           # The theoretical measurement uncertainty
    
    # The weights for the chunks based on the scatter in their time-series.
    wt0 = 1./sig**2.   # The weight from the 'sigmas' for each chunk.
    wt1 = np.zeros(vel_bc.shape)
    
    # Now loop over the observations to compute the final results
    for i in range(nr_obs):
        
        # [From iSONG:]
        # We use the reweight function on the 'dev' for the chunks. This
        # will give the chunks with small deviations higher weight.        
        wd = reweight(dev[i,:], alpha, beta, sigma)
        
        # Check for NaNs and zeros and correct to 0.01 (why 0.01?!)
        ind = np.where(np.logical_or(wd == 0., np.isnan(wd)))
        wd[ind] = 0.01
        
        # Now correct the sigma accordingly. The chunk weights are then
        # simply the inverse of the square of the sigmas.
        sig_corr = sig / wd
        weight_corr = 1.0 / sig_corr**2
        wt1[i] = weight_corr
        
        # Correct the chunk velocities of this observation by subtracting
        # their respective chunk timeseries offsets
        chunk_vels_corr = velocities[i,:] - offsets_chunk
        
        # A simple, unweighted estimate of the RV timeseries is just a
        # median of these corrected chunk velocities
        mdvel[i] = np.nanmedian(chunk_vels_corr)
        
        # The weighted RV timeseries takes the chunk weights into account
        rv[i] = np.nansum(chunk_vels_corr * weight_corr) / np.nansum(weight_corr)
        
        # The theoretical measurement uncertainty should be the inverse 
        # square-root of the sum of all weights
        errvel[i] = 1. / np.sqrt(np.nansum(weight_corr))
        
        # The chunk-to-chunk velocity scatter is the robust std of the 
        # corrected chunk velocities
        c2c_scatter[i] = robust_std(chunk_vels_corr)
        
        
        rv_bc[i] = rv[i] + bvc[i]
        
        ww[i,:] = weight_corr / np.nansum(weight_corr)
        
        
        