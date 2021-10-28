#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:18:22 2021

@author: pheeren
"""

# Import packages
from pyodine.lib.misc import printLog, chauvenet_criterion
from pyodine import timeseries
from pyodine.timeseries.misc import robust_mean, robust_std

import os
import numpy as np
import matplotlib.pyplot as plt

import argparse

# Use importlib for more flexibility of which pyodine parameter file to import?!
import importlib


def combine_velocity_results(Pars, res_files=None, comb_res_in=None, 
                             diag_file=None, plot_dir=None, comb_res_out=None, 
                             vels_out=None, reject_files=None):
    """Weight and combine chunk velocities from modelling results
    
    :param Pars: The parameters to use in the routine.
    :type Pars: :class:`Timeseries_Parameters`
    :param res_files: A pathname to a text-file with pathnames of individual 
        results to load for the combination, or, alternatively, a list of 
        pathnames to individual results. If this is None, hand an existing 
        saved CombinedResults object to 'comb_res_in'!
    :type res_files: str, list, tuple, or None
    :param comb_res_in: A pathname to a saved CombinedResults object to load. 
        If this is None, hand individual results to 'res_files'!
    :type comb_res_in: str, or None
    :param diag_file: The pathname of a text-file to write diagnosis messages 
        into. If None, the messages are just printed in the terminal.
    :type diag_file: str, or None
    :param plot_dir: The directory name to save analysis plots into. If it does
        not exist, the directory is created in the process. If None, no plots
        are saved.
    :type plot_dir: str, or None
    :param comb_res_out: The pathname where to save the final CombinedResults 
        object into. If None, the results are not saved.
    :type comb_res_out: str, or None
    :param vels_out: The pathname of a text-file to write chosen timeseries 
        results into. If None, no results are written.
    :type vels_out: str, or None
    :param reject_files: A pathname to a text-file with pathnames of individual 
        results to reject from the combination, or, alternatively, a list of 
        pathnames to individual results. If None, all results are used in the 
        combination algorithm.
    :type reject_files: str, list, tuple, or None
    
    :return: The final CombinedResults object, containing the timeseries 
        results.
    :rtype: :class:`CombinedResults`
    """
    
    
    ###########################################################################
    ## Set up the environment, and load all neccessary data and parameters
    ###########################################################################
    
    # Set up the CombinedResults object, and load from file
    # Either load a list of individual fit results, handed directly through
    # res_files or in a text-file, or load a previously saved CombinedResults 
    # object if a filename has been supplied through comb_res_in
    Results = timeseries.base.CombinedResults()
    
    if isinstance(res_files, str):
        with open(res_files, 'r') as f:
            res_names = [l.strip() for l in f.readlines()]
        Results.load_individual_results(res_names)
        
    elif isinstance(res_files, (list,tuple)):
        res_names = res_files
        Results.load_individual_results(res_names)
        
    elif isinstance(comb_res_in, str):
        Results.load_combined(comb_res_in)
        
    else:
        raise ValueError('Either hand individual fit results through "res_files"' +
                         'as list or tuple or in a text-file, or an existing' +
                         'CombinedResults object through "comb_res_in"!')
    
    # Final output name for the diagnosis file (setup the directory structure 
    # if non-existent)
    if isinstance(diag_file, str):
        diag_file_dir = os.path.dirname(diag_file)
        if not os.path.exists(diag_file_dir):
            os.makedirs(diag_file_dir)
    
    # Output directory for plots (setup the directory structure if non-existent)
    if isinstance(plot_dir, str):
        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
    
    # Final output name for the CombinedResults object  (setup the directory 
    # structure if non-existent)
    if isinstance(comb_res_out, str):
        comb_res_dir = os.path.dirname(comb_res_out)
        if not os.path.exists(comb_res_dir):
            os.makedirs(comb_res_dir)
    
    # Output name for the RV text file (in .vels format or any other defined in
    # the parameter input file) (setup the directory structure if non-existent)
    if isinstance(vels_out, str):
        vels_out_dir = os.path.dirname(vels_out)
        if not os.path.exists(vels_out_dir):
            os.makedirs(vels_out_dir)
    
    # Load a list of files that should be rejected in the timeseries
    if isinstance(reject_files, (list,tuple)):
        reject_names = reject_files
    elif isinstance(reject_files, str):
        with open(reject_files, 'r') as f:
            reject_names = [l.strip() for l in f.readlines()]
    
    
    ###########################################################################
    ## Now do the velocity weighting and combination, as prescribed in the 
    ## parameter input file
    ###########################################################################
    
    if isinstance(reject_names, (list,tuple)):
        if Pars.reject_type == 'obs_names':
            Results.remove_observations(obs_names=reject_names)
        else:
            Results.remove_observations(res_names=reject_names)
    
    Results.create_timeseries(weighting_pars=Pars.weighting_pars, 
                              diag_file=diag_file, do_crx=Pars.do_crx)
    
    if isinstance(vels_out, str):
        Results.results_to_txt(vels_out, outkeys=Pars.txt_outkeys, 
                               delimiter=Pars.txt_delimiter, header=Pars.txt_header,
                               outformat=Pars.txt_outformat)
    
    ###########################################################################
    ## Possibly save the CombinedResults object and create analysis plots
    ###########################################################################
    
    if Pars.save_comb_res and isinstance(comb_res_out, str):
        print('\nSaving results to:\n\t{}'.format(comb_res_out))
        Results.save_combined(comb_res_out)
    
    if Pars.plot_analysis and isinstance(plot_dir, str):
        print('\nCreating and saving analysis plots to\n\t{}:'.format(plot_dir))
        
        # Plot velocity results
        fig = plt.figure(figsize=(10,6))
        plt.errorbar(Results.bary_date, Results.rv_bc, yerr=Results.rv_err, 
                     fmt='.', alpha=0.7, label='Weighted velocities\nstd={:.2f} m/s'.format(
                             np.nanstd(Results.rv_bc)))
        plt.plot(Results.bary_date, Results.mdvel+Results.bary_vel_corr, 
                 '.', alpha=0.5, label='Median velocities\nstd={:.2f} m/s'.format(
                         np.nanstd(Results.mdvel+Results.bary_vel_corr)))
        plt.legend()
        plt.xlabel('JD')
        plt.ylabel('RV [m/s]')
        plt.title('{}, RV time series'.format(Results.info['star_name']))
        plt.savefig(os.path.join(plot_dir, 'RV_timeseries.png'), format='png', dpi=300)
        plt.close()
        
        # Same as above, but outliers rejected
        mask_rvs, good_rvs, bad_rvs = chauvenet_criterion(Results.rv_bc)
        rv_good = Results.rv_bc[good_rvs]
        bjd_good = Results.bary_date[good_rvs]
        rv_err_good = Results.rv_err[good_rvs]
        mdvel_good = Results.mdvel[good_rvs]
        bvc_good = Results.bary_vel_corr[good_rvs]
        
        fig = plt.figure(figsize=(10,6))
        plt.errorbar(bjd_good, rv_good, yerr=rv_err_good, fmt='.', alpha=0.7,
                     label='Weighted velocities\nstd={:.2f} m/s)'.format(np.nanstd(rv_good)))
        plt.plot(bjd_good, mdvel_good+bvc_good, '.', alpha=0.5,
                 label='Median velocities\nstd={:.2f} m/s)'.format(np.nanstd(mdvel_good+bvc_good)))
        plt.legend()
        plt.xlabel('JD')
        plt.ylabel('RV [m/s]')
        plt.title('{}, RV time series, without {} outliers'.format(
                Results.info['star_name'], len(bad_rvs[0])))
        plt.savefig(os.path.join(plot_dir, 'RV_timeseries_goodobs.png'), format='png', dpi=300)
        plt.close()
        
        # Print the outliers to file, if desired
        if Pars.print_outliers:
            printLog(diag_file, '\nObservations with outlier RVs:')
            for i in range(len(bad_rvs[0])):
                printLog(diag_file, Results.res_filename[bad_rvs[0][i]])
        
        # Plot chunk-to-chunk scatter of observations
        fig = plt.figure(figsize=(10,6))
        plt.plot(Results.bary_date, Results.c2c_scatter, '.', alpha=0.7, 
                 label='std={:.2f} m/s'.format(np.nanstd(Results.c2c_scatter)))
        plt.legend()
        plt.xlabel('JD')
        plt.ylabel('Chunk scatter [m/s]')
        plt.title('{}, chunk scatter of observations'.format(Results.info['star_name']))
        plt.savefig(os.path.join(plot_dir, 'c2c_scatter.png'), format='png', dpi=300)
        plt.close()
        
        # Plot of chunk sigmas
        fig = plt.figure(figsize=(10,6))
        plt.plot(Results.auxiliary['sigma'], '.', alpha=0.7, 
                 label='Mean: {:.2f}+-{:.2f} m/s'.format(
                         robust_mean(Results.auxiliary['sigma']), 
                         robust_std(Results.auxiliary['sigma'])))
        plt.legend()
        plt.xlabel('Chunks')
        plt.ylabel('Chunk sigmas [m/s]')
        plt.title('{}, chunk sigmas'.format(Results.info['star_name']))
        plt.savefig(os.path.join(plot_dir, 'chunk_sigma.png'), format='png', dpi=300)
        plt.close()
        
        # Plot of chunk-to-chunk offsets
        fig = plt.figure(figsize=(10,6))
        plt.plot(Results.auxiliary['chunk_offsets'], '.', alpha=0.7,
                 label='Mean: {:.2f}+-{:.2f}'.format(
                         robust_mean(Results.auxiliary['chunk_offsets']), 
                         robust_std(Results.auxiliary['chunk_offsets'])))
        plt.legend()
        plt.xlabel('Chunks')
        plt.ylabel('Chunk offsets [m/s]')
        plt.title('{}, chunk offsets from observation means'.format(Results.info['star_name']))
        plt.savefig(os.path.join(plot_dir, 'chunk_offsets.png'), format='png', dpi=300)
        plt.close()
        
        # 3D plot of velocities (corrected by chunk offsets & barycentric velocities)
        vel_corrected = Results.params['velocity'] - Results.auxiliary['chunk_offsets']
        vel_corrected = vel_corrected.T + Results.bary_vel_corr
        
        fig = plt.figure(figsize=(10,10))
        plt.imshow(vel_corrected, aspect='auto')
        plt.colorbar()
        plt.xlabel('Observations')
        plt.ylabel('Chunks')
        plt.title('{}, offset-BV-corrected chunk velocities'.format(Results.info['star_name']))
        plt.savefig(os.path.join(plot_dir, 'chunk_vels_corr.png'), format='png', dpi=300)
        plt.close()
        
        # 3D plot of chunk deviations
        fig = plt.figure(figsize=(10,10))
        plt.imshow(Results.auxiliary['chunk_dev'].T, aspect='auto')
        plt.colorbar()
        plt.xlabel('Observations')
        plt.ylabel('Chunks')
        plt.title('{}, chunk deviations'.format(Results.info['star_name']))
        plt.savefig(os.path.join(plot_dir, 'chunk_devs.png'), format='png', dpi=300)
        plt.close()
        
        # Histogram of the final velocity weights
        fig = plt.figure(figsize=(10,6))
        plt.hist(Results.auxiliary['chunk_weights'].flatten(), bins=100, alpha=0.7,
                 label=r'Mean: {}+-{} (m/s)$^{-2}$'.format(
                         robust_mean(Results.auxiliary['chunk_weights']), 
                         robust_std(Results.auxiliary['chunk_weights'])))
        plt.legend()
        plt.xlabel(r'Weights [(m/s)$^{-2}$]')
        plt.title('{}, chunk weights'.format(Results.info['star_name']))
        plt.savefig(os.path.join(plot_dir, 'chunk_weights_hist.png'), format='png', dpi=300)
        plt.close()
        
        # Plot chromatic indices (if any)
        if Pars.do_crx:
            fig = plt.figure(figsize=(10,6))
            plt.errorbar(Results.bary_date, Results.crx, yerr=Results.crx_err, 
                         fmt='.', alpha=0.7, label='std={:.2f} m/s'.format(
                                 np.nanstd(Results.crx)))
            plt.legend()
            plt.xlabel('JD')
            plt.ylabel('CRX [(m/s)/Np]')
            plt.title('{}, CRX time series'.format(Results.info['star_name']))
            plt.savefig(os.path.join(plot_dir, 'CRX_timeseries.png'), format='png', dpi=300)
            plt.close()
        
    ###########################################################################
    ## Everything's done now, return the CombinedResults object
    ###########################################################################
    
    print('All done!')
    
    return Results


if __name__ == '__main__':
    
    # Set up the parser for input arguments
    parser = argparse.ArgumentParser(
            description='Weight and combine velocities from observation modelling')
    
    # Required input arguments:
    # utilities_dir, ostar_files, temp_files, temp_outname, (plot_dir=None, par_file=None)
    parser.add_argument('par_file', type=str, help='The pathname to the timeseries parameters file to use.')
    parser.add_argument('--res_files', type=str, help='A pathname to a text-file with the pathnames of modelling results.')
    parser.add_argument('--comb_res_in', type=str, help='The pathname to a saved CombinedResults object.')
    parser.add_argument('--diag_file', type=str, help='The pathname of a text-file to write diagnosis messages into.')
    parser.add_argument('--plot_dir', type=str, help='The pathname to a directory where to save analysis plots.')
    parser.add_argument('--comb_res_out', type=str, help='The pathname where to save the CombinedResults object.')
    parser.add_argument('--vels_out', type=str, help='The pathname of a text-file where to save chosen timeseries results.')
    parser.add_argument('--reject_files', type=str, help='A pathname of a text-file with the pathnames of results to reject.')
    
    # Parse the input arguments
    args = parser.parse_args()
    
    par_file = args.par_file
    res_files = args.res_files
    comb_res_in = args.comb_res_in
    diag_file = args.diag_file
    plot_dir = args.plot_dir
    comb_res_out = args.comb_res_out
    vels_out = args.vels_out
    reject_files = args.reject_files
    
    # Import and load the timeseries parameters
    par_file = os.path.splitext(par_file)[0].replace('/', '.')
    timeseries_parameters = importlib.import_module(par_file)
    Pars = timeseries_parameters.Timeseries_Parameters()
    
    # And run the velocity weighting routine
    combine_velocity_results(Pars, res_files=res_files, comb_res_in=comb_res_in, 
                             diag_file=diag_file, plot_dir=plot_dir, 
                             comb_res_out=comb_res_out, vels_out=vels_out, 
                             reject_files=reject_files)
