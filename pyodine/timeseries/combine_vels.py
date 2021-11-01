#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 15:45:44 2021

@author: pheeren
"""

import numpy as np
from .misc import robust_mean, robust_std, reweight
from ..lib.misc import printLog
import lmfit

"""The default _weighting_pars array contains a combination of values for the 
weighting parameters that have proven to work well for the computation of RVs
from SONG spectra. Particularly the 'good_chunks' and 'good_orders' will have 
to be changed when using a different instrument!
"""
_weighting_pars = {
        'reweight_alpha': 1.8,
        'reweight_beta': 8.0,
        'reweight_sigma': 2.0,
        'weight_correct': 0.01,
        'sig_limit_low': 4., 
        'sig_limit_up': 1000.,
        'sig_correct': 1000.,
        'good_chunks': None, #(3, 15), #(150, 350)
        'good_orders': None #(6,14)
        }


def combine_chunk_velocities(velocities, nr_chunks_order, bvc=None, 
                             wavelengths=None, diag_file=None, 
                             weighting_pars=None):
    """Weight and combine the chunk velocities of a modelled timeseries
    
    This routine follows the algorithm used in the original iSONG pipeline 
    code, developed by Frank Grundahl.
    ToDo: Check the plausibility of the final RV uncertainties!
    Added here: The chromatic index (slope of the modelled velocities with
    wavelength) is computed if an array of wavelength zero points for each
    chunk in each observation is supplied.
    
    :param velocities: The modelled velocities for each chunk in each 
        observation of the timeseries.
    :type velocities: ndarray[nr_obs,nr_chunks]
    :param nr_chunks_order: Number of chunks per order.
    :type nr_chunks_order: int
    :param bvc: If barycentric velocity corrections are supplied for all 
        observations, the output also contains bvc-corrected RVs.
    :type bvc: list, ndarray[nr_obs], or None
    :param wavelengths: If an array of wavelength zeropoints for each chunk
        in each observation is supplied, the chromatic indices (crx) of the 
        observations are modelled.
    :type wavelength: ndarray[nr_obs,nr_chunks], or None
    :param diag_file: If the pathname to a text-file is supplied, diagnosis
        output is printed there in addition to terminal output.
    :type diag_file: str, or None
    :param weighting_pars: A dictionary of parameters used in the weighting 
        algorithm. If None is supplied, a dictionary of default values is used.
    :type weighting_pars: dict, or None
    
    :return: A dictionary of results: RVs ('rvs'), BVC-corrected RVs 
        ('rvs_bc', optionally), simple median velocity timeseries ('mdvels'),
        RV uncertainties ('rv_errs'), chunk-to-chunk scatter within each
        observation ('c2c_scatters'), the chromatic index of each observation 
        ('crxs', optional), the uncertainty of the chromatic indices 
        ('crx_errs', optional).
    :rtype: dict
    :return: A dictionary of auxiliary results: The chunk timeseries sigmas 
        ('chunk_sigma'), individual chunk deviations ('chunk_dev'), chunk 
        timeseries offsets from observation medians ('chunk_offsets'), 
        corrected chunk weights ('chunk_weights'), and measures of the achieved 
        RV precision of the timeseries ('RV_precision1', 'RV_precision2').
    :rtype: dict
    """
    
    if not isinstance(weighting_pars, dict):
        pars = _weighting_pars
    else:
        pars = weighting_pars.copy()
    
    printLog(diag_file, '---------------------------------------------------')
    printLog(diag_file, '- Pyodine chunk combination (based on iSONG code) -')
    printLog(diag_file, '---------------------------------------------------\n')
    
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
    
    # Where are chunk velocities nan?
    ind_nan = np.where(np.isnan(velocities))
    if len(ind_nan[0]) > 0:
        printLog(diag_file, '')
        printLog(diag_file, 'Nan velocities (obs,chunk):')
        outstring = ''
        for i in range(len(ind_nan[0])):
            outstring += '({},{})  '.format(ind_nan[0][i], ind_nan[1][i])
        printLog(diag_file, outstring)
    printLog(diag_file, '')
    
    # For each observation, center the chunk velocities around 0 by
    # subtracting the robust mean of that observation
    # For the robust mean: Only use the best chunk velocities
    # (parameters 'good_chunks' & 'good_orders', if available):
    if isinstance(pars['good_orders'], (list,tuple,np.ndarray)) and \
    isinstance(pars['good_chunks'], (list,tuple,np.ndarray)) and \
    len(pars['good_orders']) == 2 and len(pars['good_chunks']) == 2:
        good_ind = []
        for o in range(pars['good_orders'][0], pars['good_orders'][1]+1):
            good_ind += [(i + nr_chunks_order*o) for i in range(pars['good_chunks'][0], pars['good_chunks'][1]+1)]
    else:
        good_ind = [i for i in range(len(velocities))]
        del pars['good_orders']
        del pars['good_chunks']
    
    vel_offset_corrected = np.zeros((nr_obs, nr_chunks))
    for i in range(nr_obs):
        vel_offset_corrected[i,:] = velocities[i,:] - robust_mean(velocities[i,good_ind]) #pars['good_chunks'][0]:pars['good_chunks'][1]])
    
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
    ind = np.where(np.logical_or(sig<=pars['sig_limit_low'], sig>=pars['sig_limit_up']))
    sig[ind] = pars['sig_correct']
    
    # Print to file
    printLog(diag_file, 'Chunk sigmas:')
    printLog(diag_file, 'Median: {:.2f} +- {:.2f}\n'.format(
            np.nanmedian(sig), np.nanstd(sig)))
    
    # Prepare the output dicts
    rv_dict = {
            'rv': np.zeros(nr_obs),             # Weighted RV timeseries
            'rv_bc': np.zeros(nr_obs),          # Weighted RV timeseries, BV-corrected
            'mdvel': np.zeros(nr_obs),          # The simple observation median (after correcting 
                                                # for chunk timeseries offsets)
            'rv_err': np.zeros(nr_obs),         # The theoretical measurement uncertainty
            'c2c_scatter': np.zeros(nr_obs),    # The chunk-to-chunk velocity scatter in each observation
            }
    if isinstance(wavelengths, (list,tuple,np.ndarray)):
        rv_dict['crx'] = np.zeros(nr_obs)           # The chromatic index in each observation
        rv_dict['crx_err'] = np.zeros(nr_obs)       # The fit errors of the chromatic indices
        rv_dict['RV_wave'] = np.zeros(nr_obs)       # The effect. wavelengths at the weighted RVs
        rv_dict['RV_wave_err'] = np.zeros(nr_obs)   # The fit errors of the effect. wavelengths
        rv_dict['crx_redchi'] = np.zeros(nr_obs)    # The red. Chi**2 of the crx fits
        
    
    chunk_weights = np.zeros((nr_obs,nr_chunks))    # The corrected weights for each individual chunk 
                                                    # (after the reweight function)
    
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
        chunk_weights[i] = 1.0 / sig_corr**2
        
        # Correct the chunk velocities of this observation by subtracting
        # their respective chunk timeseries offsets
        chunk_vels_corr = velocities[i,:] - offsets_chunk
        
        # A simple, unweighted estimate of the RV timeseries is just a
        # median of these corrected chunk velocities
        rv_dict['mdvel'][i] = np.nanmedian(chunk_vels_corr)
        
        # The weighted RV timeseries takes the chunk weights into account
        rv_dict['rv'][i] = np.nansum(chunk_vels_corr * chunk_weights[i]) / np.nansum(chunk_weights[i])
        
        # The theoretical measurement uncertainty should be the inverse 
        # square-root of the sum of all weights
        rv_dict['rv_err'][i] = 1. / np.sqrt(np.nansum(chunk_weights[i]))
        
        # The chunk-to-chunk velocity scatter is the robust std of the 
        # corrected chunk velocities
        rv_dict['c2c_scatter'][i] = robust_std(chunk_vels_corr)
        
        # BV-correction instead through actual z and barycorrpy?!
        if isinstance(bvc, (list,np.ndarray)):
            rv_dict['rv_bc'][i] = rv_dict['rv'][i] + bvc[i]
        
        # Compute the chromatic index in the observation, which is the slope
        # of the chunk velocities over chunk wavelengths
        # ToDo: Use the corrected chunk velocities or the raw ones?!
        # ToDo: Use weights?!
        if isinstance(wavelengths, (list,tuple,np.ndarray)):
            crx_dict = chromatic_index_observation(
                    #vel_offset_corrected[i], wavelengths[i], rv_dict['rv'][i], weights=chunk_weights[i])
                    #chunk_vels_corr, wavelengths[i], rv_dict['rv'][i], weights=chunk_weights[i])
                    velocities[i,:], wavelengths[i], rv_dict['rv'][i], weights=chunk_weights[i])
            
            for key, value in crx_dict.items():
                rv_dict[key][i] = value
        
    # Some metrics of the quality of the RV timeseries
    rv_precision1 = np.sqrt(1./np.nansum(1./sig**2.))
    rv_precision2 = np.sqrt(1./np.nansum(np.nanmedian(chunk_weights[i], axis=0)))
    
    printLog(diag_file, 'RV quality factor 1 ( sqrt(1/sum(1/sig**2)) ): {} m/s'.format(
            rv_precision1))
    printLog(diag_file, 'RV quality factor 2 ( sqrt(1/sum(med(wt1))) ): {} m/s'.format(
            rv_precision2))
    
    auxiliary_dict = {
            'chunk_sigma': sig,
            'chunk_dev': dev,
            'chunk_offsets': offsets_chunk,
            'chunk_weights': chunk_weights,
            'rv_precision1': rv_precision1,     # The inverse root of all summed simple weights (1/sigma**2)
            'rv_precision2': rv_precision2,     # The inverse root of all summed observation-median
                                                # corrected weights
            }
    
    return rv_dict, auxiliary_dict, weighting_pars


def velocity_from_chromatic_index(wavelengths, RV, RV_wave, crx):
    """The function to evaluate chunk velocities, given certain crx parameters
    
    :param wavelengths: Input wavelengths at which to evaluate the velocities.
    :type wavelengths: ndarray[nr_chunks]
    :param RV: The weighted RV of the observation.
    :type RV_zero: float
    :param RV_wave: The effective wavelength at the weighted observation RV.
    :type RV_wave: float
    :param crx: The chromatic index of the observation.
    :type crx: float
    
    :return: The evaluated chunk velocities from the crx model.
    :rtype: ndarray[nr_chunks]
    """
    return RV + crx * np.log(wavelengths / RV_wave)
    


def chromatic_index_observation(velocities, wavelengths, RV, weights=None):
    """Model the chromatic index (crx) of an observation
    
    This follows the idea as implemented in SERVAL (Zechmeister et al., 2018; 
    see Section 3.1 and Equation 21).
    
    :param velocities: The modelled velocities of all chunk.
    :type velocities: ndarray[nr_chunks]
    :param wavelengths: The modelled wavelength intercepts (zeropoints) of all
        chunks.
    :type wavelengths: ndarray[nr_chunks]
    :param RV: The weighted RV of the observation.
    :type RV: float
    :param weights: An optional array of chunk weights to use in the modelling
        of the crx (e.g. to downweight chunks which are inherently bad). If 
        None, no weights are used in the fitting.
    :type weights: ndarray[nr_chunks], or None
    
    :return: A dictionary of results: chromatic index ('crx') and its model
        uncertainty ('crx_err'), the effective wavelength at the weighted RV
        ('RV_wave') and its model uncertainty ('RV_wave_err'), and the red. 
        Chi**2 of the fit ('crx_redchi').
    :rtype: dict
    """
    
    def func(lmfit_params, wavelengths, velocities, RV, weights):
        """The objective function for the crx model
        
        :param lmfit_params: The crx parameters.
        :type lmfit_params: :class:`lmfit.Parameters`
        :param wavelengths: The chunk wavelengths.
        :type wavelengths: ndarray[nr_chunks]
        :param velocities: The chunk velocities.
        :type velocities: ndarray[nr_chunks]
        :param RV: The weighted RV of the observation.
        :type RV: float
        :param weights: The weights of the chunks.
        :type weights: ndarray[nr_chunks] (or None)
        
        :return: The (optionally weighted) residuals between the crx model
            and the velocities.
        :rtype: ndarray[nr_chunks]
        """
        
        velocities_fit = velocity_from_chromatic_index(
                wavelengths, RV, lmfit_params['RV_wave'], lmfit_params['crx'])
        if len(np.where(np.isnan(velocities_fit))[0]) > 0:
            print(velocities_fit)
        if isinstance(weights, (list,np.ndarray)):
            return (velocities_fit - velocities) * np.sqrt(np.abs(weights))
        else:
            return velocities_fit - velocities
    
    # Take care of NaNs in the wavelengths or velocities arrays: Exclude them
    w_nan_inds = np.where(np.isnan(wavelengths))
    if len(w_nan_inds[0]) > 0:
        w_fin_inds = np.where(np.isfinite(wavelengths))
        wavelengths = wavelengths[w_fin_inds]
        velocities = velocities[w_fin_inds]
        
    v_nan_inds = np.where(np.isnan(velocities))
    if len(v_nan_inds[0]) > 0:
        v_fin_inds = np.where(np.isfinite(velocities))
        wavelengths = wavelengths[v_fin_inds]
        velocities = velocities[v_fin_inds]
    
    # Now do the parameter starting guesses
    # For the effective wavelength of the modelled RV:
    # Median of all wavelengths
    RV_wave_guess = np.median(wavelengths)
    # For the chromatic index:
    # Slope between the first (30) and last (30) chunk velocities
    vels_first = np.median(velocities[:30])
    vels_last  = np.median(velocities[-30:])
    waves_first = np.median(wavelengths[:30])
    waves_last  = np.median(wavelengths[-30:])
    
    # This follows when rearranging Equ. 21 in Zechmeister et al. (2018)
    crx_guess = (vels_last - vels_first) / np.log(waves_last / waves_first)
    
    # Set up the the lmfit parameters
    lmfit_params = lmfit.Parameters()
    lmfit_params.add('RV_wave', value=RV_wave_guess, min=np.min(wavelengths), 
                     max=np.max(wavelengths))
    lmfit_params.add('crx', value=crx_guess)
    
    # And fit
    lmfit_result = lmfit.minimize(func, lmfit_params, args=[wavelengths, velocities, RV, weights])
    
    crx_dict = {
            'crx': lmfit_result.params['crx'].value,
            'crx_err': lmfit_result.params['crx'].stderr, 
            'RV_wave': lmfit_result.params['RV_wave'].value, 
            'RV_wave_err': lmfit_result.params['RV_wave'].stderr, 
            'crx_redchi': lmfit_result.redchi
            }
    
    return crx_dict