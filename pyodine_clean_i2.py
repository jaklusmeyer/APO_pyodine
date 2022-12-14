#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 17:21:06 2020

@author: pheeren
"""

# Import packages
import pyodine

import os
#import sys
import time
import numpy as np
import h5py
import logging
#from pathos.multiprocessing import Pool
import traceback
from progressbar import ProgressBar

#import argparse
#import importlib

# Dictionary for LSF smoothing parameters
_smooth_dict = {
        'smooth_pixels': 160,
        'smooth_orders': 3,
        'order_separation': 15
        }


def clean_i2_single_observation(utilities, Pars, res_file, out_file, rv=None,
                                plot_dir=None, error_log=None, info_log=None, 
                                quiet=False, use_progressbar=False,
                                smooth_lsf=True, lsf_dict=_smooth_dict):
    """Clean a single observation spectrum from I2 features, following the
    algorithm described in Diaz +2019
    
    This routine loads a model results of an observation and aims at cleaning
    the spectrum from the I2 features. The result and analysis plots can be 
    saved to file.
    
    :param utilities: The utilities module for the instrument used in this 
        analysis.
    :type utilities: library
    :param Pars: The parameter input object to use.
    :type Pars: :class:`Parameters`
    :param res_file: The pathname of the model result to use.
    :type res_file: str
    :param out_file: The pathname where to save the resulting cleaned spectrum.
    :type out_file: str
    :param rv: Optionally the RV to use in the template correction can be 
        supplied here. If None, the chunk velocities are used (default).
    :type rv: float, int, or None
    :param plot_dir: The directory name where to save plots. If the directory 
        structure does not exist yet, it will be created in the process. If 
        None is given, no plots will be saved (default).
    :type plot_dir: str, or None
    ::param error_log: A pathname of a log-file used for error messages. If 
        None, no errors are logged.
    :type error_log: str, or None
    :param info_log: A pathname of a log-file used for info messages. If 
        None, no info is logged.
    :type info_log: str, or None
    :param quiet: Whether or not to print info messages to terminal. Defaults 
        to False (messages are printed).
    :type quiet: bool
    :param use_progressbar: Whether to show a progressbar during the algorithm. 
        Defaults to False.
    :type use_progressbar: bool
    :param smooth_lsf: Whether to smooth the LSFs over neighbouring chunks.
        Defaults to True.
    :type smooth_lsf: bool
    :param lsf_dict: A dictionary with parameters for the LSF smoothing. Should
        contain entries for 'smooth_pixels', 'smooth_orders' and 
        'order_separation'. If None is given, the default dictionary 
        _smooth_dict is used.
    :type lsf_dict: dict
    """
    
    # Check whether a logger is already setup. If no, setup a new one
    #if not logging.getLogger().hasHandlers():
    pyodine.lib.misc.setup_logging(
            config_file=Pars.log_config_file, level=Pars.log_level,
            error_log=error_log, info_log=info_log, quiet=quiet)
    
    try:
        
        # Start timer
        start_t = time.time()
        
        logging.info('')
        logging.info('---------------------------')
        logging.info('Aiming to I2-clean: {}'.format(res_file))
        
        ###########################################################################
        ## Set up the environment, and load all neccessary data and parameters
        ###########################################################################
        
        # Load result
        obs_chunks, fit_results = pyodine.fitters.results_io.restore_results_object(
                utilities, res_file)
        
        # Output pathname
        if not os.path.exists(os.path.dirname(res_file)):
            os.makedirs(os.path.dirname(res_file))
        
        # Output directory for plots (setup the directory structure if non-existent)
        if isinstance(plot_dir, str):
            if not os.path.exists(plot_dir):
                os.makedirs(plot_dir)
        
        
        ###########################################################################
        ## Now prepare the cleaning: Setup the (smoothed) LSF and the model
        ###########################################################################
        
        # Smooth the LSFs if desired
        if smooth_lsf:
            
            redchi2 = np.array([r.redchi for r in fit_results])
            
            lsf_smoothed = pyodine.lib.misc.smooth_lsf(
                    obs_chunks, lsf_dict['smooth_pixels'], lsf_dict['smooth_orders'], 
                    lsf_dict['order_separation'], fit_results,
                    redchi2=redchi2, osample=fit_results[0].model.osample_factor)
            
            LSFarr = pyodine.models.lsf.LSF_Array(
                    lsf_smoothed, np.array([ch.order for ch in obs_chunks]),
                    np.array([ch.abspix[0] for ch in obs_chunks]))
            
            
            # Build the model and fitter
            lsf_model = pyodine.models.lsf.model_index['FixedLSF']
            model = pyodine.models.spectrum.SimpleModel(
                lsf_model, fit_results[0].model.wave_model, fit_results[0].model.cont_model, 
                fit_results[0].model.iodine_atlas, stellar_template=fit_results[0].model.stellar_template, 
                lsf_array=LSFarr, osample_factor=fit_results[0].model.osample_factor, 
                conv_width=fit_results[0].model.conv_width)
        
        else:
            model = fit_results[0].model
        
        ###########################################################################
        ## Now loop over the chunks and clean them to return the I2-free spectra
        ###########################################################################
        
        # Now build the I2-free spectrum, following the receipe described in Diaz +2019
        i2_free_specs = []
        
        # Use a progressbar?
        if use_progressbar:
            bar = ProgressBar(max_value=len(obs_chunks), redirect_stdout=True)
            bar.update(0)
        
        for i, chunk in enumerate(obs_chunks):
            
            # If LSF smoothing: Setup the new parameters to use here
            if smooth_lsf:
                # First the LSF parameters
                pars = pyodine.models.base.ParameterSet(
                        {'lsf_order': chunk.order, 
                         'lsf_pixel0': chunk.abspix[0], 
                         'lsf_amplitude': 1.}
                        )
                # Now copy best-fit parameters for non-LSF parameters
                for key, value in fit_results[i].params.items():
                    if 'lsf' not in key:
                        pars[key] = value
            
            else:
                pars = fit_results[i].params
            
            # And finally compute the I2-cleaned spectrum
            i2_free_flux = model.clean_of_I2(
                    chunk, pars, require='full', chunk_ind=i, rv=rv)
            
            # Build a spectrum object out of it
            wave = model.wave_model.eval(chunk.pix, pars.filter('wave'))
            cont = model.cont_model.eval(chunk.pix, pars.filter('cont'))
            
            i2_free_spec = pyodine.components.Spectrum(i2_free_flux, wave, cont)
            
            i2_free_specs.append(i2_free_spec)
            
            # Update the progressbar
            if use_progressbar:
                bar.update(i+1)
        
        if use_progressbar:
            bar.finish()
        
        
        ###########################################################################
        ## To Do: Stitch the chunks together?!
        ##
        ## For now: Save them as HDF5
        ###########################################################################
        
        # Saving...
        
        with h5py.File(out_file, 'w') as h:
            h['flux'] = [spec.flux for spec in i2_free_specs]
            h['wave'] = [spec.wave for spec in i2_free_specs]
            h['cont'] = [spec.cont for spec in i2_free_specs]
            
            info_dict = {
                    'orig_filename': os.path.abspath(obs_chunks[0].observation.orig_filename).encode('utf8', 'replace'),
                    'res_name': res_file.encode('utf8', 'replace'),
                    'orig_header': obs_chunks[0].observation.orig_header.tostring(sep='\n').encode('utf8', 'replace'),
                    'star_name': obs_chunks[0].observation.star.name.encode('utf8', 'replace'),
                    'instrument_name': obs_chunks[0].observation.instrument.name.encode('utf8', 'replace'),
                    'smooth_lsf': smooth_lsf,
                    'smooth_dict': lsf_dict,
                    }
            
            pyodine.lib.h5quick.dict_to_group(info_dict, h, 'info')
        
        ###########################################################################
        ## And Done!
        ###########################################################################
        
        modelling_time = time.time() - start_t
        logging.info('')
        logging.info('Time to model this observation: {}'.format(modelling_time))
    
    except Exception as e:
        """
        with open(error_file, 'a') as f:
            f.write(parameters[5]+'\n')
        """
        traceback.print_exc()
        logging.info('Something went wrong with input file {}'.format(res_file), 
                     exc_info=True)
        print()
        print(e)
        
