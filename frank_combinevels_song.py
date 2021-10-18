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

# This is where the results sit
#resdir = '../data_lick/data_res/hip88048/multigausslick_ch_norecenter/hip88048all_4o_flat_nonorm_tempwc_nobpm/'
resdir = '../data_song/data_res/sun/multigauss_song_new/sun_flat_lsfc_lsfnew/collected_results.h5'
#resdir = '../Documents/song/arcturus/collected_test.h5'

# Savepath of the original Lick velocity file
#orig_vst_file = '../data_lick/data_res/vsthip88048.vels'
orig_vst_file = None #'../data_song/data_res/vsthip102488_song_5.vels'

# Savepath of the Lick bvc file
#bvc_file = '../data_lick/bcvel.ascii'
#bvc_file = '../data_song/HIP102488_2020-09-01T085227.txt'

# This is where the analysis plots go
save_dir_name = 'plots_song_analysis_bvc0'
if os.path.isdir(resdir):
    savepath = os.path.join(resdir, save_dir_name)
else:
    savepath = os.path.join(os.path.dirname(resdir), save_dir_name)
if not os.path.isdir(savepath):
    os.mkdir(savepath)

reject_file_name = 'sun_reject.vels'
if os.path.isdir(resdir) and reject_file_name is not None:
    reject_file = os.path.join(resdir, reject_file_name)
elif reject_file_name is not None:
    reject_file = os.path.join(os.path.dirname(resdir), reject_file_name)
else:
    reject_file = None

alpha = 1.8
beta = 8.0
sigma = 2.0

"""
If a directory path is supplied, collect all results for one star, write into one data cube,
and load it.
If the name of the data cube is supplied, load it.
"""

if os.path.isdir(resdir):
    print('Collecting individual results...')
    dirs = [d for d in os.listdir(resdir) if os.path.isdir(os.path.join(resdir,d))]
    files = []
    for d in dirs:
        fs = os.listdir(os.path.join(resdir,d))
        for f in fs:
            if 'results1.h5' in f:
                files.append(os.path.join(resdir, d, f))
    files.sort()
    
    Res_comb = timeseries.base.CombinedResults(files)
    Res_comb.save_combined(os.path.join(resdir, 'collected_results.h5'))
    print('Combined results file created.')
    
    Res_comb = timeseries.base.CombinedResults(os.path.join(resdir, 'collected_results.h5'))

else: 
    try:
        print('Loading combined results file...')
        Res_comb = timeseries.base.CombinedResults(resdir)
    except:
        raise IOError('File not found, or could not be read!')

"""
Get the required parameters from the combined results cube.
"""

star_name = Res_comb.info['star_name']
vel = Res_comb.params['velocity']
vel_err = Res_comb.errors['velocity']
redchi2 = Res_comb.redchi2
medcnts = Res_comb.medcnts
bjd = Res_comb.observation['bary_date']
bvc = Res_comb.observation['bary_vel_corr']
res_filename = Res_comb.res_filenames

if reject_file:
    reject = np.genfromtxt(reject_file, unpack=True, dtype=str)
    if reject.shape == ():
        reject = [str(reject)]
else:
    reject = []

if len(reject) > 0:
    print('Using info from reject file...')
    ind = []
    for i in range(len(res_filename)):
        if res_filename[i] not in reject:
            ind.append(i)
    star_name = Res_comb.info['star_name']
    vel = Res_comb.params['velocity'][ind]
    vel_err = Res_comb.errors['velocity'][ind]
    redchi2 = Res_comb.redchi2[ind]
    medcnts = Res_comb.medcnts[ind]
    bjd = Res_comb.observation['bary_date'][ind]
    bvc = Res_comb.observation['bary_vel_corr'][ind]
    #bvc = np.zeros(bjd.shape)
    res_filename = [res_filename[i] for i in ind]

"""
date_ind = np.where(np.logical_and(bjd>2457522.3, bjd<2457522.38))
vel = vel[date_ind]
vel_err = vel_err[date_ind]
bjd = bjd[date_ind]
bvc = bvc[date_ind]
res_filename = [res_filename[i] for i in date_ind[0]]
"""

"""
# Without the edge velocities
ind = np.where((np.arange(704) % 44 < 10) | (np.arange(704) % 44 > 34))
vel = vel[:,ind[0]]
vel_err = vel_err[:,ind[0]]
redchi2 = redchi2[:,ind[0]]
medcnts = medcnts[:,ind[0]]
"""
# Name of the diagnosis file
diag_file = os.path.join(savepath, 'py{}diag.txt'.format(star_name))

printLog(diag_file, '--------------------------------------------------')
printLog(diag_file, '- Pyodine chunk combination (based on SONG code) -')
printLog(diag_file, '--------------------------------------------------')
printLog(diag_file, '')

nr_obs_all, nr_chunks = vel.shape
printLog(diag_file, 'Nr. of obs, chunks per obs: {}, {}'.format(nr_obs_all, nr_chunks))

# Give an overview of counts and redchi2 values
printLog(diag_file, 'Med. counts: {} +- {}'.format(np.nanmedian(medcnts), np.nanstd(medcnts)))
printLog(diag_file, 'Med. red. Chi^2: {} +- {}'.format(np.nanmedian(redchi2), np.nanstd(redchi2)))

# Where are chunk velocities or barycentric velocities nan?
ind_nan = np.where(np.isnan(vel))
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

nr_obs = nr_obs_all

"""
Here the actual velocity weighting starts.
"""

# Set up the barycentric corrected velocities
vel_bc = np.transpose(np.transpose(vel) + bvc)

# Calculate d2 array
d2 = np.zeros((nr_obs_all, nr_chunks))
for i in range(nr_obs):
    d2[i,:] = vel[i,:] - robust_mean(vel[i,150:350])

# Calculate relative chunk-to-chunk offsets
offsets_chunk = robust_mean(d2, axis=0)
offsets_chunk -= robust_mean(d2)

printLog(diag_file, '')
printLog(diag_file, 'Chunk-to-chunk offsets from observation mean:')
printLog(diag_file, 'Median: {:.2f} +- {:.2f}'.format(
        np.nanmedian(offsets_chunk), np.nanstd(offsets_chunk)))

# Sigma and deviation arrays
sig = np.zeros(nr_chunks)
dev = np.zeros((nr_obs_all, nr_chunks))

for j in range(nr_chunks):
    sig[j] = robust_std(d2[:,j])
    dev[:,j] = (d2[:,j] - robust_mean(d2[:,j])) / sig[j]

for i in range(nr_obs):
    dev[i,:] = dev[i,:] - robust_mean(dev[i,:])

ind = np.where(np.logical_or(sig<=4., sig>=1000.))
sig[ind] = 1000.

printLog(diag_file, '')
printLog(diag_file, 'Chunk sigmas:')
printLog(diag_file, 'Median: {:.2f} +- {:.2f}'.format(
        np.nanmedian(sig), np.nanstd(sig)))

# Loop over observations
rv = np.zeros(nr_obs)  # Weighted velocity for each observation
rv_bc = np.zeros(nr_obs)
c2c_scatter = np.zeros(nr_obs)
ww = np.zeros(vel_bc.shape)
mdvel = np.zeros(nr_obs)
errvel = np.zeros(nr_obs)

# The weights for the chunks based on the scatter in their time-series.
wt0 = 1./sig**2.   # The weight from the 'sigmas' for each chunk.
wt1 = np.zeros(vel_bc.shape)

for i in range(nr_obs):
    
    # [From iSONG:]
    # We multiply the 'sigma' for the chunks by the re-weight function.
    # This will give the points with small deviations from the mean curve for 
    # the chunk higher weight. In this way we combine the global (sig) sigmas
    # with down-weighting for this particular image, if needed.
    
    wd = reweight(dev[i,:], alpha, beta, sigma)
    # TODO: Check for NaN's and correct to 0.01 m/s
    ind = np.where(np.logical_or(wd == 0., np.isnan(wd)))
    
    wd[ind] = 0.01
    
    sig1 = sig / wd
    
    wt = 1.0 / sig1**2
    wt1[i] = wt
    
    # Overwrite weights
    #wt = 1.0 / sig**2
    
    
    velocity = vel[i,:] - offsets_chunk
    mdvel[i] = np.nanmedian(velocity)
    rv[i] = np.nansum(velocity * wt) / np.nansum(wt)
    rv_bc[i] = rv[i] + bvc[i]
    c2c_scatter[i] = robust_std(velocity)
    ww[i,:] = wt / np.nansum(wt)
    errvel[i] = 1./np.sqrt(np.nansum(wt))

print('Weights 1:', np.sqrt(1./np.nansum(wt0)))

print('Weights 2:', np.sqrt(1./np.nansum(np.median(wt1, axis=0))))


"""
Save RVs and other information
"""

# Write RV results to a file in vst format
with open(os.path.join(savepath, 'pyvst{}.vels'.format(star_name)), 'w') as f:
    for i in range(nr_obs):
        outstring = '{:.5f}\t{:.4f}\t{:.4f}\n'.format(
                bjd[i], rv_bc[i], errvel[i])
        f.write(outstring)

# Append the results from the weighting algorithm to the CombinedResult data cube
#if date_ind:
#    Res_comb = timeseries.base.CombinedResults(res_filename)

Res_comb.tseries = {}

Res_comb.tseries['c2c_scatter'] = c2c_scatter
Res_comb.tseries['rv'] = rv_bc
Res_comb.tseries['wt0'] = wt0
Res_comb.tseries['wt1'] = wt1
Res_comb.tseries['drv'] = errvel

Res_comb.save_combined(os.path.join(savepath, 'collected_results.h5'))

"""
Plots about velocities, sigmas, devs and weights.
"""

# 3D plot of velocities
fig = plt.figure(figsize=(10,8))
plt.imshow(vel_bc, aspect='auto')
plt.colorbar()
plt.xlabel('Chunks')
plt.ylabel('Observations')
plt.title('Velocities, BV-corrected')
plt.savefig(os.path.join(savepath, 'vel_cube0.png'), format='png', dpi=300)
plt.close()

# 3D plot of velocities corrected by chunk offsets
fig = plt.figure(figsize=(10,8))
plt.imshow(d2, aspect='auto')
plt.colorbar()
plt.xlabel('Chunks')
plt.ylabel('Observations')
plt.title('d2')
plt.savefig(os.path.join(savepath, 'vel_cube1.png'), format='png', dpi=300)
plt.close()

# Plot of chunk-to-chunk offsets
fig = plt.figure(figsize=(10,8))
plt.plot(offsets_chunk, '.', alpha=0.75)
plt.legend(['Chunk offset median: {:.2f} +- {:.2f}'.format(
        np.nanmedian(offsets_chunk), np.nanstd(offsets_chunk))])
plt.xlabel('Chunks')
plt.ylabel('Mean offset [m/s]')
plt.title('Chunk-to-chunk offsets from obs. means')
plt.savefig(os.path.join(savepath, 'chunk_offsets.png'), format='png', dpi=300)
plt.close()

# Plot of chunk sigmas
fig = plt.figure(figsize=(10,8))
plt.plot(sig, '.', alpha=0.75)
plt.legend(['Sigma median: {:.2f} +- {:.2f}'.format(
        np.nanmedian(sig), np.nanstd(sig))])
plt.xlabel('Chunks')
plt.ylabel('Chunk sigma [m/s]')
plt.title('Chunk series deviation from obs. medians')
plt.savefig(os.path.join(savepath, 'chunk_sigma.png'), format='png', dpi=300)
plt.close()

# 3D plot of chunk deviations
fig = plt.figure(figsize=(10,8))
plt.imshow(dev, aspect='auto')
plt.colorbar()
plt.xlabel('Chunks')
plt.ylabel('Observations')
plt.title('Chunk deviations')
plt.savefig(os.path.join(savepath, 'vel_cube2.png'), format='png', dpi=300)
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

# Plot velocity results
fig = plt.figure(figsize=(10,8))
plt.errorbar(bjd, rv_bc-np.median(rv_bc), yerr=errvel, fmt='.', alpha=0.75, label='Weighted velocities')
plt.plot(bjd, mdvel+bvc-np.median(mdvel+bvc), '.', alpha=0.75, label='Median velocities')
plt.legend(['Median velocities\n(std={:.2f} m/s)'.format(np.nanstd(mdvel+bvc)),
           'Weighted velocities\n(std={:.2f} m/s)'.format(np.nanstd(rv_bc))])
plt.xlabel('JD')
plt.ylabel('RV [m/s]')
plt.title('Velocity time series')
plt.savefig(os.path.join(savepath, 'velocity_series0.png'), format='png', dpi=300)
plt.close()


# Now subtract straight line
linear_model = np.polyfit(bjd, rv_bc-np.median(rv_bc), 1)
linear_model_fn = np.poly1d(linear_model)

fig = plt.figure(figsize=(10,8))
plt.errorbar(bjd, rv_bc-np.median(rv_bc)-linear_model_fn(bjd), yerr=errvel, fmt='.', alpha=0.75, label='Weighted velocities')
plt.plot(bjd, mdvel+bvc-np.median(mdvel+bvc), '.', alpha=0.75, label='Median velocities')
plt.legend(['Median velocities\n(std={:.2f} m/s)'.format(np.nanstd(mdvel+bvc)),
           'Weighted velocities\n(std={:.2f} m/s)'.format(np.nanstd(rv_bc-linear_model_fn(bjd)))])
plt.xlabel('JD')
plt.ylabel('RV [m/s]')
plt.title('Velocity time series')
plt.savefig(os.path.join(savepath, 'velocity_series0_corr.png'), format='png', dpi=300)
plt.close()


# Same as above, but outliers rejected
mask_rvs, good_rvs, bad_rvs = chauvenet_criterion(rv_bc)
rv_good = rv_bc[good_rvs]
bjd_good = bjd[good_rvs]
errvel_good = errvel[good_rvs]

fig = plt.figure(figsize=(10,8))
plt.errorbar(bjd_good, rv_good-np.median(rv_good), yerr=errvel_good, fmt='.', alpha=0.75),
plt.plot(bjd_good, mdvel[good_rvs]+bvc[good_rvs]-np.median(mdvel[good_rvs]+bvc[good_rvs]), '.', alpha=0.75)
plt.legend(['Median velocities\n(std={:.2f} m/s)'.format(np.nanstd(mdvel[good_rvs]+bvc[good_rvs])),
           'Weighted velocities\n(std={:.2f} m/s)'.format(np.nanstd(rv_good))])
plt.xlabel('JD')
plt.ylabel('RV [m/s]')
plt.title('Velocity time series, without {} outliers'.format(len(bad_rvs[0])))
plt.savefig(os.path.join(savepath, 'velocity_series0_restricted.png'), format='png', dpi=300)
plt.close()

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
