#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:18:22 2021

@author: pheeren
"""

# Import packages
from pyodine.lib.misc import printLog, chauvenet_criterion
from pyodine import timeseries
from pyodine.timeseries.misc import robust_mean, robust_std, reweight

import os
import numpy as np
import matplotlib.pyplot as plt



def combine_velocity_results(Pars, res_files=None, comb_res_in=None, 
                             diag_file=None, plot_dir=None, comb_res_out=None, 
                             vels_out=None, reject_files=None):
    
    
    
    
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
        Results.save_combined(comb_res_out)
    
    if Pars.plot_analysis and isinstance(plot_dir, str):
        
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
        


# Histogram of the final velocity weights
fig = plt.figure(figsize=(10,8))
plt.hist(wt1.flatten(), bins=200, alpha=0.75)
plt.legend(['Median weights (only non-0): \n{} +- {} m/s'.format(
        np.nanmedian(wt1[wt1!=0.]), np.nanstd(wt1[wt1!=0.]))])
plt.xlabel('Weights')
plt.title('Chunk weights after velocity rms weighting (0 for rejected chunks)')
plt.savefig(os.path.join(savepath, 'weights_hist1.png'), format='png', dpi=300)
plt.close()

"""
Now comes the analysis with plots
"""

# Print outliers to analysis file
printLog(diag_file, 'Outlier velocities:')
for i in range(len(bad_rvs[0])):
    printLog(diag_file, res_filename[bad_rvs[0][i]])
"""
# Load original Lick velocities
jd_orig = np.genfromtxt(orig_vst_file, unpack=True, usecols=[0])
vel_orig = np.genfromtxt(orig_vst_file, unpack=True, usecols=[1])
vel_orig -= np.median(vel_orig)
#vel_orig = np.array([-0., -0., 0., 0., 0.])
err_orig = np.genfromtxt(orig_vst_file, unpack=True, usecols=[2])


# Plot again velocity results, this time vs. original Lick results
fig = plt.figure(figsize=(10,8))
plt.errorbar(bjd, rv-np.nanmedian(rv), yerr=errvel, fmt='.', alpha=0.75)
plt.errorbar(jd_orig, vel_orig, yerr=err_orig, fmt='.', alpha=0.75)
plt.legend(['Weighted velocities\n(std={:.2f} m/s)'.format(np.nanstd(rv)), 
            'Original velocities\n(std={:.2f} m/s)'.format(np.nanstd(vel_orig))])
plt.xlabel('JD')
plt.ylabel('RV [m/s]')
plt.title('Velocity time series with original Lick results')
plt.savefig(os.path.join(savepath, 'velocity_series1.png'), format='png', dpi=300)
plt.close()


# Same as above, but outliers rejected
fig = plt.figure(figsize=(10,8))
plt.errorbar(bjd_good, rv_good-np.nanmedian(rv_good), yerr=errvel_good, fmt='.', alpha=0.75)
plt.errorbar(jd_orig, vel_orig, yerr=err_orig, fmt='.', alpha=0.75)
plt.legend(['Weighted velocities\n(std={:.2f} m/s)'.format(np.nanstd(rv_good)), 
            'Original velocities\n(std={:.2f} m/s)'.format(np.nanstd(vel_orig))])
plt.xlabel('JD')
plt.ylabel('RV [m/s]')
plt.title('Velocity time series (without outliers) with original Lick results')
plt.savefig(os.path.join(savepath, 'velocity_series1_restricted.png'), format='png', dpi=300)


# Find original velocities that correspond to ours and plot differences
ind_ours = []
ind_orig = []
for i in range(len(bjd)):
    ind = np.argmin(np.abs(bjd[i]-jd_orig))
    if abs(bjd[i] - jd_orig[ind]) < 0.01:
        ind_ours.append(i)
        ind_orig.append(ind)

jd_res = bjd[ind_ours]
full_res = rv[ind_ours] - vel_orig[ind_orig]
full_res -= np.nanmedian(full_res)
full_res_err = np.sqrt(errvel[ind_ours]**2 + err_orig[ind_orig]**2)

fig = plt.figure(figsize=(10,8))
plt.errorbar(jd_res, full_res, yerr=full_res_err, fmt='.', alpha=0.75)
plt.legend(['Residual velocities\n(std={:.2f} m/s)'.format(np.nanstd(full_res))])
plt.xlabel('JD')
plt.ylabel('Residual RV [m/s]')
plt.title('Difference to original Lick results')
plt.savefig(os.path.join(savepath, 'velocity_difference.png'), format='png', dpi=300)
plt.close()
"""
"""
# Diss plot
fig = plt.figure(figsize=(8,7))
plt.errorbar(bjd, rv-vel_orig-np.nanmean(rv-vel_orig),
             yerr=np.sqrt(errvel**2 + err_orig**2), fmt='p', alpha=0.75)
plt.plot(bjd, mdvel+bvc-vel_orig-np.nanmean(mdvel+bvc-vel_orig), 'o', alpha=0.75)
plt.legend(['Unweighted RVs\n(rms={:.2f} m/s)'.format(
            np.nanstd(mdvel+bvc-vel_orig-np.nanmean(mdvel+bvc-vel_orig))),
    'Weighted RVs\n(rms={:.2f} m/s)'.format(
    np.nanstd(rv-vel_orig-np.nanmean(rv-vel_orig)))
    ])
plt.axhline(y=0, c='g', linestyle='--', alpha=0.5)
plt.xlabel('JD')
plt.ylabel('$v_\mathrm{meas} - v_\mathrm{inj}$ [m/s]')
plt.title('Residuals measured vs. injected RVs')
plt.savefig(os.path.join(savepath, 'rv_weighted_median_diff.png'), format='png', dpi=300)
plt.close()

# Plot distribution of velocities
fig = plt.figure(figsize=(8,7))
velocity_med = np.nanmedian(vel[-1])
velocity_std = np.nanstd(vel[-1])
plt.hist(vel[-1], bins=200)#, alpha=0.75)#, range=(0,2000));
plt.xlabel('Velocity [m/s]')
#plt.xlim((0,400.))
plt.axvline(x=velocity_med, color='r', alpha=0.5)
plt.axvline(x=velocity_med-velocity_std, color='g', alpha=0.5)
plt.axvline(x=velocity_med+velocity_std, color='g', alpha=0.5)
plt.legend(['Med = {:5.1f} m/s\nStd = {:5.1f} m/s'.format(velocity_med, velocity_std)])
plt.title('Chunk velocities for $v_\mathrm{inj} = 2000 \,\mathrm{m/s}$')
#plt.xlim(-3000,500)
plt.savefig(os.path.join(savepath, 'vel_hist_4.png'),
                         format='png', dpi=300)
plt.close()
"""
"""
# Same as above but without outliers
ind_ours2 = []
ind_orig2 = []
for i in range(len(bjd_good)):
    ind = np.argmin(np.abs(bjd_good[i]-jd_orig))
    if abs(bjd_good[i] - jd_orig[ind]) < 0.2:
        ind_ours2.append(i)
        ind_orig2.append(ind)

jd_good_res = bjd_good[ind_ours2]
good_res = rv_good[ind_ours2] - vel_orig[ind_orig2]
good_res -= np.nanmedian(good_res)
good_res_errs = np.sqrt(errvel_good[ind_ours2]**2 + err_orig[ind_orig2]**2)

fig = plt.figure(figsize=(10,8))
plt.errorbar(jd_good_res, good_res, yerr=good_res_errs, fmt='.', alpha=0.75)
plt.legend(['Residual velocities\n(std={:.2f} m/s)'.format(np.nanstd(good_res))])
plt.xlabel('JD')
plt.ylabel('Residual RV [m/s]')
plt.title('Difference to original Lick results, without outliers')
plt.savefig(os.path.join(savepath, 'velocity_difference_restricted.png'), format='png', dpi=300)
plt.close()


# Compare date and bvc of this pipeline with original Lick values

# Plot difference in Julian Date between this and the original Lick data
fig = plt.figure(figsize=(10,8))
plt.plot(jd_orig[ind_orig], bjd[ind_ours]-jd_orig[ind_orig], '.', alpha=0.75)
plt.xlabel('Original JD')
plt.ylabel('JD difference [d]')
plt.title('Difference between our JDs and the original ones')
plt.savefig(os.path.join(savepath, 'JD_difference_to_lickvst.png'), format='png', dpi=300)
plt.close()


# Plot correlation of velocity residuals to barycentric velocities
fig = plt.figure(figsize=(10,8))
plt.plot(np.abs(full_res), bvc[ind_ours], '.', alpha=0.75)
plt.xlabel('Residual vels [m/s]')
plt.ylabel('BVC [m/s]')
plt.title('Correlation between RV residuals and BVCs')
plt.savefig(os.path.join(savepath, 'RVres_vs_bvcs.png'), format='png', dpi=300)
plt.close()


# Plot correlation of velocity residuals to date residuals
fig = plt.figure(figsize=(10,8))
plt.plot(np.abs(full_res), bjd[ind_ours]-jd_orig[ind_orig], '.', alpha=0.75)
plt.xlabel('Residual vels [m/s]')
plt.ylabel('Residual JDs')
plt.title('Correlation between RV residuals and JD residuals')
plt.savefig(os.path.join(savepath, 'RVres_vs_jdres.png'), format='png', dpi=300)
plt.close()
"""
"""
# Get original barycentric velocity values for this star
bcvel_all = np.genfromtxt(bvc_file, skip_header=1, unpack=True, dtype=str)

bcvels = bcvel_all[:,np.where((bcvel_all[1,:]==star_name[3:]) | (bcvel_all[1,:]==star_name))[0]]
print(bcvels.shape)

bcvels_v = bcvels[2,:].astype(np.float)
bcvels_d = bcvels[3,:].astype(np.float) + 2440000.


# Show dates from this routine and the original one
fig = plt.figure(figsize=(10,8))
plt.vlines(bcvels_d, ymin=0., ymax=0.5, label='lick', color='r', linewidth=0.7)
plt.vlines(bjd, ymin=0.5, ymax=1., label='pyodine', color='b', linewidth=0.7)
plt.legend()
plt.xlabel('JD')
plt.title('Date signatures of pyodine and orig. Lick')
plt.savefig(os.path.join(savepath, 'JD_comparison_lickbvc.png'), format='png', dpi=300)
plt.close()


# Compare barycentric dates from this routine to original ones
ind_pyodine = []
ind_lickbvc = []
for i in range(len(bjd)):
    ind = np.argmin(np.abs(bjd[i]-bcvels_d))
    if abs(bjd[i] - bcvels_d[ind]) < .5:
        ind_pyodine.append(i)
        ind_lickbvc.append(ind)

fig = plt.figure(figsize=(10,8))
plt.plot(bcvels_d[ind_lickbvc], bcvels_d[ind_lickbvc]-bjd[ind_pyodine], '.', alpha=0.75)
plt.xlabel('Lick BVC JD')
plt.ylabel('JD difference [d]')
plt.title('Difference between pyodine JDs and Lick BVC dates')
plt.savefig(os.path.join(savepath, 'JD_difference_to_lickbvc.png'), format='png', dpi=300)
plt.close()


# Residuals from barycentric velocities from this routine to original ones
fig = plt.figure(figsize=(10,8))
plt.plot(bjd[ind_pyodine], bvc[ind_pyodine]-bcvels_v[ind_lickbvc], '.', alpha=0.75)
plt.xlabel('pyodine JD')
plt.ylabel('BVC difference [m/s]')
plt.title('Difference between pyodine BVCs and Lick BVCs')
plt.savefig(os.path.join(savepath, 'BVC_difference_to_lickbvc.png'), format='png', dpi=300)
plt.close()
"""
