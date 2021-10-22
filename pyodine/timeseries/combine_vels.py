#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 15:45:44 2021

@author: pheeren
"""

import numpy as np
from .misc import robust_mean, robust_std, reweight
from ..lib.misc import printLog


def combine_chunk_velocities(velocities, bvc, diag_file=None, pars=None):
    """Here the actual velocity weighting starts.
    """
    
    if not isinstance(pars, dict):
        pars = {
                'reweight_alpha': 1.8,
                'reweight_beta': 8.0,
                'reweight_sigma': 2.0,
                'weight_correct': 0.01,
                'sig_limits': (4., 1000.),
                'sig_correct': 1000.,
                'good_chunks': (150, 350)
                }
    
    printLog(diag_file, '--------------------------------------------------')
    printLog(diag_file, '- Pyodine chunk combination (based on SONG code) -')
    printLog(diag_file, '--------------------------------------------------')
    printLog(diag_file, '')
    
    printLog(diag_file, 'Weighting parameters used:')
    for key, value in pars.items():
        if isinstance(value, (list,tuple)):
            print_str = '\t{}\t'.format(key)
            for i in range(len(value)):
                print_str += str(value[i]) + '\t'
        else:
            print_str = '\t{}\t{}'.format(key, value)
        printLog(diag_file, print_str)
    printLog(diag_file, '')
    
    # How many observations and chunks?
    (nr_obs, nr_chunks) = velocities.shape
    printLog(diag_file, 'Nr. of obs, chunks per obs: {}, {}'.format(nr_obs, nr_chunks))
    
    # Where are chunk velocities or barycentric velocities nan?
    ind_nan = np.where(np.isnan(velocities))
    ind_nan_bvc = np.where(np.isnan(bvc))
    if len(ind_nan[0]) > 0:
        printLog(diag_file, '')
        printLog(diag_file, 'Nan velocities (obs,chunk):')
        outstring = ''
        for i in range(len(ind_nan[0])):
            outstring += '({},{})  '.format(ind_nan[0][i], ind_nan[1][i])
        printLog(diag_file, outstring)
    if len(ind_nan_bvc[0]) > 0:
        printLog(diag_file, '')
        printLog(diag_file, 'Nan barycentric velocities (obs):')
        outstring = ''
        for i in range(len(ind_nan[0])):
            outstring += '({})\t'.format(ind_nan_bvc[0][i])
        printLog(diag_file, outstring)
    printLog(diag_file, '')
    
    # Set up the barycentric corrected chunk velocities
    vel_bc = np.transpose(np.transpose(velocities) + bvc)
    
    # For each observation, center the chunk velocities around 0 by
    # subtracting the robust mean of that observation
    vel_offset_corrected = np.zeros((nr_obs, nr_chunks))
    for i in range(nr_obs):
        vel_offset_corrected[i,:] = velocities[i,:] - robust_mean(
                velocities[i,pars['good_chunks'][0]:pars['good_chunks'][1]])
    
    # Calculate the mean offset of each individual chunk timeseries from
    # the observation means by taking the robust mean along the 
    # (offset-corrected) chunks
    offsets_chunk = robust_mean(vel_offset_corrected, axis=0)
    
    # Subtract the overall robust mean of all chunk velocities in all 
    # observations (after the chunk velocities were already offset-corrected,
    # so in fact this should barely make a difference)
    offsets_chunk -= robust_mean(vel_offset_corrected)
    
    # Print to file
    printLog(diag_file, 'Chunk-to-chunk offsets from observation mean:')
    printLog(diag_file, 'Median: {:.2f} +- {:.2f}\n'.format(
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
    ind = np.where(np.logical_or(sig<=pars['sig_limits'][0], sig>=pars['sig_limits'][1]))
    sig[ind] = pars['sig_correct']
    
    # Print to file
    printLog(diag_file, 'Chunk sigmas:')
    printLog(diag_file, 'Median: {:.2f} +- {:.2f}\n'.format(
            np.nanmedian(sig), np.nanstd(sig)))
    
    # Prepare the output dict
    rv_dict = {
            'rvs': np.zeros(nr_obs),            # Weighted RV timeseries
            'rvs_bc': np.zeros(nr_obs),         # Weighted RV timeseries, BV-corrected
            'mdvels': np.zeros(nr_obs),         # The simple observation median (after correcting 
                                                # for chunk timeseries offsets)
            'rv_errs': np.zeros(nr_obs),        # The theoretical measurement uncertainty
            'c2c_scatters': np.zeros(nr_obs),   # The chunk-to-chunk velocity scatter in each observation
            'crx': np.zeros(nr_obs),            # The chromatic index in each observation
            }
    
    
    # The weights for the chunks based on the scatter in their time-series.
    wt0 = 1./sig**2.   # The weight from the 'sigmas' for each chunk.
    wt1 = np.zeros(vel_bc.shape)
    
    # Now loop over the observations to compute the final results
    for i in range(nr_obs):
        
        # [From iSONG:]
        # We use the reweight function on the 'dev' for the chunks. This
        # will give the chunks with small deviations higher weight.        
        wd = reweight(dev[i,:], pars['reweight_alpha'], 
                      pars['reweight_beta'], pars['reweight_sigma'])
        
        # Check for NaNs and zeros and correct to 0.01 (why 0.01?!)
        ind = np.where(np.logical_or(wd == 0., np.isnan(wd)))
        wd[ind] = pars['weight_correct']
        
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
        rv_dict['mdvels'][i] = np.nanmedian(chunk_vels_corr)
        
        # The weighted RV timeseries takes the chunk weights into account
        rv_dict['rvs'][i] = np.nansum(chunk_vels_corr * weight_corr) / np.nansum(weight_corr)
        
        # The theoretical measurement uncertainty should be the inverse 
        # square-root of the sum of all weights
        rv_dict['rv_errs'][i] = 1. / np.sqrt(np.nansum(weight_corr))
        
        # The chunk-to-chunk velocity scatter is the robust std of the 
        # corrected chunk velocities
        rv_dict['c2c_scatters'][i] = robust_std(chunk_vels_corr)
        
        # BV-correction instead through actual z and barycorrpy?!
        rv_dict['rvs_bc'][i] = rv_dict['rvs'][i] + bvc[i]
        
    # Some metrics of the quality of the RV timeseries
    rv_quality1 = np.sqrt(1./np.nansum(wt0))
    rv_quality2 = np.sqrt(1./np.nansum(np.nanmedian(wt1, axis=0)))
    
    printLog(diag_file, 'RV quality factor 1 ( sqrt(1/sum(sig**-2)) ): {} m/s'.format(rv_quality1))
    printLog(diag_file, 'RV quality factor 2 ( sqrt(1/sum(med(wt1))) ): {} m/s'.format(rv_quality2))
    
    return rv_dict
        