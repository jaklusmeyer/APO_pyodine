#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 10:15:46 2020

@author: pheeren
"""

import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from os import path
from .template.base import StellarTemplate_Chunked
from .timeseries.misc import robust_std

def plot_chunkmodel(fit_results, chunk_array, chunk_nr, template=True, tellurics=None, 
                    weight=None, title='', savename=None, dpi=300, show_plot=False):
    """Create a plot of the template, iodine, observed and modeled spectrum
    and residuals between the two for a given chunk. If it is for an
    O-star observation, set template to False. If tellurics are given,
    shade the regions affected by them. If savename is given,
    the plot is saved.
    
    Args:
        fit_results (list): A list containing instances of :class:'LmfitResult'
            for each chunk.
        chunk_array (:class:'ChunkArray'): The chunks of the observation.
        chunk_nr (int): The index of the chunk to plot.
        template (Optional[bool]): True if the template should be plotted,
            False if not (e.g. in O-star models). Defaults to True.
        tellurics (Optional[:class:'SimpleTellurics']): An instance of tellurics.
            If None, they are not included (default).
        weight (Optional[ndarray]): An array of pixel weights, to display which
            pixels were excluded in the modelling. Defaults to None.
        title (Optional[str]): A title for the plot. If None, a default title
            is used.
        savename (Optional[str]): If a pathname is given, the plot is saved
            there. Defaults to None.
        dpi (Optional[int]): DPI of the saved plot. Defaults to 300.
        show_plot (Optional[bool]): If True, the plot is showed during execution.
            Defaults to False.
    """
    # Evaluate model
    try:
        sp = fit_results[chunk_nr].fitted_spectrum
        iod = fit_results[chunk_nr].model.iodine_atlas
        
        if template:
            fig = plt.figure(figsize=(10,10))
            nr_axes = 4
            gs = fig.add_gridspec(4, 1,  height_ratios=(1, 1, 1, 0.5))
        else:
            fig = plt.figure(figsize=(10,8.5))
            gs = fig.add_gridspec(3, 1,  height_ratios=(1, 1, 0.5))
            nr_axes = 3
        
        ax = []
        
        if template:
            # Compute Doppler shift of this chunk
            beta = fit_results[chunk_nr].params['velocity'] / 299792458.
            doppler = np.sqrt((1. + beta) / (1. - beta))
            if isinstance(fit_results[chunk_nr].model.stellar_template, StellarTemplate_Chunked):
                temp = fit_results[chunk_nr].model.stellar_template[chunk_nr]
                temp_shifted = fit_results[chunk_nr].model.stellar_template[chunk_nr]
            else:
                # Get template range for this chunk, and doppler shifted template
                temp = fit_results[chunk_nr].model.stellar_template.get_wavelength_range(sp.wave[0], sp.wave[-1])
                temp_shifted = fit_results[chunk_nr].model.stellar_template.get_wavelength_range(sp.wave[0]/doppler, sp.wave[-1]/doppler)
            #ax.append(fig.add_subplot(nr_axes,1,1))
            ax.append(fig.add_subplot(gs[0]))
            # Plot modeled template
            ax[-1].plot(temp.wave, temp.flux, drawstyle='steps-mid')
            ax[-1].plot(temp_shifted.wave*doppler, temp_shifted.flux, drawstyle='steps-mid')
            ax[-1].set_xlim((sp.wave[0], sp.wave[-1]))
            plt.legend(['Template', 'Shifted template (v={})'.format(int(fit_results[chunk_nr].params['velocity']))])
            ax[-1].tick_params(axis='both', top=True, right=True, direction='in', which='both', labelbottom=False)
        
        #ax.append(fig.add_subplot(nr_axes,1,nr_axes-2))
        ax.append(fig.add_subplot(gs[nr_axes-3]))
        # Plot iodine
        iod_chunk = iod.get_wavelength_range(sp.wave[0], sp.wave[-1])
        ax[-1].plot(iod_chunk.wave, iod_chunk.flux)
        ax[-1].set_xlim((sp.wave[0], sp.wave[-1]))
        ax[-1].tick_params(axis='both', top=True, right=True, direction='in', which='both', labelbottom=False)
        
        #ax.append(fig.add_subplot(nr_axes,1,nr_axes-1))
        ax.append(fig.add_subplot(gs[nr_axes-2]))
        # Plot observation
        ax[-1].plot(sp.wave, chunk_array[chunk_nr].flux, drawstyle='steps-mid')
        # Plot model (spectrum and continuum)
        ax[-1].plot(sp.wave, sp.flux, drawstyle='steps-mid')
        ax[-1].plot(sp.wave, sp.cont)
        ax[-1].set_xlim((sp.wave[0], sp.wave[-1]))
        plt.legend(['Observation', 'Fitted model', 'Continuum'])
        ax[-1].tick_params(axis='both', top=True, right=True, direction='in', which='both', labelbottom=False)
        
        #ax.append(fig.add_subplot(nr_axes,1,nr_axes))
        ax.append(fig.add_subplot(gs[nr_axes-1]))
        # Residual plot (normalize with mean flux)
        ax[-1].plot(sp.wave, fit_results[chunk_nr].residuals / sp.flux, drawstyle='steps-mid')
        ax[-1].set_xlim((sp.wave[0], sp.wave[-1]))
        ax[-1].tick_params(axis='both', top=True, right=True, direction='in', which='both')
        
        # If tellurics are given, add them as shaded regions
        if tellurics is not None:
            ind = np.where((tellurics.dict_tellurics['wave_stop']>sp.wave[0]) & 
                           (tellurics.dict_tellurics['wave_start']<sp.wave[-1]))
            for tell_line in ind[0]:
                for j in range(len(ax)):
                    ax[j].axvspan(tellurics.dict_tellurics['wave_start'][tell_line],
                                tellurics.dict_tellurics['wave_stop'][tell_line], alpha=0.1, color='k')
        if weight is not None:
            ind2 = np.where(weight != 0.)
            ind3 = np.where(weight == 0.)
            
            if len(ind3[0]) > 0:
                ax[-1].plot(sp.wave[ind3], np.zeros(len(ind3[0])), 'P', color='r', alpha=0.5)
                plt.legend(['rms={:.3f}%'.format((robust_std(fit_results[chunk_nr].residuals/sp.flux)*1e2)),
                            'Weights = 0',
                            'rms_c={:.3f}%'.format(
                                    (robust_std(fit_results[chunk_nr].residuals[ind2]/ \
                                               sp.flux[ind2])*1e2))])
            else:
                plt.legend(['rms={:.3f}%'.format((robust_std(fit_results[chunk_nr].residuals/sp.flux)*1e2))])
        else:
            plt.legend(['rms={:.3f}%'.format((robust_std(fit_results[chunk_nr].residuals/sp.flux)*1e2))])
        
        fig.subplots_adjust(hspace=0)
        ax[-1].set_xlabel('Wavelength [$\AA$]')
        
        if title == '':
            title = 'Chunk: {}, Order: {}, Pixels: {} - {}'.format(
                    chunk_nr, chunk_array[chunk_nr].order,
                    chunk_array[chunk_nr].abspix[0], chunk_array[chunk_nr].abspix[-1])
        ax[0].set_title(title)
        
        if savename is not None:
            fmt = path.splitext(savename)[1][1:]
            plt.savefig(savename, format=fmt, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()
    except Exception as e:
        print(e)


def plot_residual_hist(fit_results, residual_arr=None, tellurics=None, title='', 
                    savename=None, dpi=300, show_plot=False):
    """Create a histogram of all chunk residuals (in percent)
    
    If tellurics are given, both the complete residuals and those outside of 
    the telluric regions are plotted. 
    
    Args:
        fit_results (list): A list containing instances of :class:'LmfitResult'
            for each chunk.
        residual_arr (Optional[tuple, list, ndarray]): If supplied, just use 
            this data for the plot (if tuple: all residuals and sub residuals, 
            otherwise just one array of residuals). If left free, calculate 
            the data from fit_results.
        tellurics (Optional[:class:'SimpleTellurics']): An instance of 
            tellurics, to disregard affected pixels in the computation of 
            the sub residuals array. Defaults to None.
        title (Optional[str]): A title for the plot. If None, a default title
            is used.
        savename (Optional[str]): If a pathname is given, the plot is saved
            there. If None, then no saving.
        dpi (Optional[int]): DPI of the saved plot. Defaults to 300.
        show_plot (Optional[bool]): If True, the plot is shown during execution.
            Defaults to False.
       
       Return:
           tuple or ndarray: The calculated residuals (either both all and sub,
                or only all).
    """
    if residual_arr is None:
        all_res = np.array([robust_std(r.residuals/r.fitted_spectrum.flux*1e2) for r in fit_results \
                            if r.lmfit_result is not None])
        if tellurics is not None:
            sub_res = []
            for r in fit_results:
                if r.lmfit_result is not None:
                    ind = np.array(tellurics.is_affected(r.fitted_spectrum.wave))
                    sub_res.append(robust_std(r.residuals[np.where(ind==0)]/r.fitted_spectrum.flux[np.where(ind==0)])*1e2)
            sub_res = np.array(sub_res)
        else:
            sub_res = None
    else:
        if type(residual_arr) is tuple:
            if len(residual_arr) == 2:
                all_res, sub_res = residual_arr
            else:
                all_res = residual_arr[0]
                sub_res = None
        else:
            all_res = residual_arr
            sub_res = None
    
    fig, ax = plt.subplots(figsize=(14,8))
    ax.hist(all_res, bins=np.arange(0,10,0.1))#np.nanmax(all_res)*0.8,0.1))
    ax.set_xlabel('Residuals [%]')
    
    axins = ax.inset_axes([0.48, 0.4, 0.5, 0.58])
    axins.hist(all_res, bins=np.arange(0,2,0.02))
    
    if sub_res is not None:
        ax.hist(sub_res, bins=np.arange(0,10,0.1), alpha=0.5)#np.nanmax(sub_res),0.1), alpha=0.5)
        axins.hist(sub_res, bins=np.arange(0,2,0.02), alpha=0.5)
        axins.legend(['All: med={:.3f}$\pm${:.3f}%'.format(np.nanmedian(all_res), np.nanstd(all_res)), 
                      'Outside tell.: med={:.3f}$\pm${:.3f}%'.format(np.nanmedian(sub_res), np.nanstd(sub_res))])
    else:
        axins.legend(['All: med={:.3f}$\pm${:.3f}%'.format(np.nanmedian(all_res), np.nanstd(all_res))])
    
    if title == '':
        title = 'Distribution of chunk residuals'
    ax.set_title(title)
    
    if savename is not None:
        fmt = path.splitext(savename)[1][1:]
        plt.savefig(savename, format=fmt, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()
    
    if sub_res is not None:
        return all_res, sub_res
    else:
        return all_res


def plot_chunk_scatter(scatter=None, scatter_fmt='o', scatter_alpha=1., 
                       scatter_label=None,
                       curve=None, curve_fmt='-', curve_alpha=1., 
                       curve_label=None,
                       errorbar=None, errorbar_yerr=None, errorbar_fmt='.', 
                       errorbar_alpha=1., errorbar_label=None,
                       hlines=None, hlines_fmt='-', hlines_alpha=1.,
                       hlines_label=None, hlines_color=None,
                       ylabel=None, ylog=False, yrange=None, grid=True, 
                       nr_chunks_order=None, nr_orders=None, title='', 
                       savename=None, dpi=300, show_plot=False):
    """Create any chunk scatter plots
    
    Plot any data for all chunks via 'scatter' and/or 'curve'.
    
    Args:
        scatter (Optional[list, ndarray, tuple]): Data that should be plotted
            as scatter. Can be 2D, that is multiple datasets for all chunks.
            All other scatter keywords only come in place if this is provided.
        scatter_fmt (Optional[str]): The format of the scatter points. Defaults
            to 'o'.
        scatter_alpha (Optional[float]): The alpha-value of the scatter points 
            (should be between 0. and 1.). Defaults to 1.
        scatter_label (Optional[list, ndarray, tuple, str]): Provide labels for
            the scatter data. Can be a list with labels for each scatter
            dataset, if multiple are given.
        curve (Optional[list, ndarray, tuple]): Data that should be plotted
            as curve. Can be 2D, that is multiple datasets for all chunks.
            All other curve keywords only come in place if this is provided.
        curve_fmt (Optional[str]): The linestyle format. Defaults to '-'.
        curve_alpha (Optional[float]): The alpha-value of the curves 
            (should be between 0. and 1.). Defaults to 1.
        curve_label (Optional[list, ndarray, tuple, str]): Provide labels for
            the curve data. Can be a list with labels for each curve
            dataset, if multiple are given.
        errorbar (Optional[list, ndarray, tuple]): Data that should be plotted
            as errorbars. Can be 2D, that is multiple datasets for all chunks.
            All other errorbar keywords only come in place if this is provided.
        errorbar_yerr (Optional[list, ndarray, tuple]): The y-errors of the
            errorbar data. Make sure the format works! If None, no errorbars
            are plotted.
        errorbar_fmt (Optional[str]): The linestyle format. Defaults to '.'.
        errorbar_alpha (Optional[float]): The alpha-value of the errorbars 
            (should be between 0. and 1.). Defaults to 1.
        errorbar_label (Optional[list, ndarray, tuple, str]): Provide labels 
            for the errorbar data. Can be a list with labels for each errorbar
            dataset, if multiple are given.
        hlines (Optional[int, float, list, ndarray, tuple]): Data that should 
            be plotted as horizontal line. If not a single number, then 
            multiple lines are plotted. All other hlines keywords only come in 
            place if this is provided.
        hlines_fmt (Optional[str]): The linestyle format. Defaults to '.'.
        hlines_alpha (Optional[float]): The alpha-value of the hlines 
            (should be between 0. and 1.). Defaults to 1.
        hlines_label (Optional[list, ndarray, tuple, str]): Provide labels 
            for the hlines. Can be a list with labels for each hlines
            dataset, if multiple are given.
        hlines_color (Optional[list, ndarray, tuple, str, int, float]): The
            color(s) of the line(s).
        ylabel (Optional[str]): The label of the y-axis. Defaults to None.
        ylog (Optional[bool]): If True, logscale y-axis. Defaults to False.
        yrange(Optional[list,tuple]): The range of the y-axis. Defaults to 
            None.
        grid (Optional[bool]): Whether or not to plot a grid. Defaults to True.
        nr_chunks_order (Optional[int]): Number of chunks per order. If this
            and nr_orders is given, the order borders are plotted.
        nr_orders (Optional[int]): Number of orders. If this and 
            nr_chunks_orders is given, the order borders are plotted.
        title (Optional[str]): A title for the plot. If None, a default title
            is used.
        savename (Optional[str]): If a pathname is given, the plot is saved
            there. If None, then no saving.
        dpi (Optional[int]): DPI of the saved plot. Defaults to 300.
        show_plot (Optional[bool]): If True, the plot is shown during execution.
            Defaults to False.
    """
    # Make sure that if some input data was wrong, the whole program does not
    # crash
    try:
        fig = plt.figure(figsize=(12,6))
        
        # If scatter data is provided, plot it
        if isinstance(scatter, (np.ndarray, list, tuple)):
            if not isinstance(scatter[0], (np.ndarray, list, tuple)):
                scatter = [scatter]
            if not isinstance(scatter_label, (np.ndarray, list, tuple, str)):
                scatter_label = [None] * len(scatter)
            elif isinstance(scatter_label, str):
                scatter_label = [scatter_label] * len(scatter)
            for i in range(len(scatter)):
                if ylog:
                    plt.semilogy(scatter[i], scatter_fmt, alpha=scatter_alpha,
                                 label=scatter_label[i])
                else:
                    plt.plot(scatter[i], scatter_fmt, alpha=scatter_alpha,
                             label=scatter_label[i])
        
        # If curve data is provided, plot it
        if isinstance(curve, (np.ndarray, list, tuple)):
            if not isinstance(curve[0], (np.ndarray, list, tuple)):
                curve = [curve]
            if not isinstance(curve_label, (np.ndarray, list, tuple, str)):
                curve_label = [None] * len(curve)
            elif isinstance(curve_label, str):
                curve_label = [curve_label] * len(curve)
            for i in range(len(curve)):
                if ylog:
                    plt.semilogy(curve[i], linestyle=curve_fmt, alpha=curve_alpha,
                                 label=curve_label[i])
                else:
                    plt.plot(curve[i], linestyle=curve_fmt, alpha=curve_alpha,
                             label=curve_label[i])
        
        # If errorbar data is provided, plot it
        if isinstance(errorbar, (np.ndarray, list, tuple)):
            if not isinstance(errorbar[0], (np.ndarray, list, tuple)):
                errorbar = [errorbar]
                # Assume that the errorbars are also just given for one dataset
                # (if at all)
                errorbar_yerr = [errorbar_yerr]
            else:
                # If errorbars are None, make sure that the plot still works
                if not isinstance(errorbar_yerr, (np.ndarray, list, tuple)):
                    errorbar_yerr = [None] * len(errorbar)
            if not isinstance(errorbar_label, (np.ndarray, list, tuple, str)):
                errorbar_label = [None] * len(errorbar)
            elif isinstance(errorbar_label, str):
                errorbar_label = [errorbar_label] * len(errorbar)
            for i in range(len(errorbar)):
                plt.errorbar(np.arange(len(errorbar[i])), errorbar[i], 
                             yerr=errorbar_yerr[i], fmt=errorbar_fmt, 
                             alpha=errorbar_alpha, label=errorbar_label[i])
            if ylog:
                plt.yscale('log')
        
        # If hlines data is provided, plot it
        if isinstance(hlines, (np.ndarray, list, tuple, int, float)):
            if isinstance(hlines, (int, float)):
                hlines = [hlines]
            if not isinstance(hlines_label, (np.ndarray, list, tuple, str)):
                hlines_label = [None] * len(hlines)
            elif isinstance(hlines_label, str):
                hlines_label = [hlines_label] * len(hlines)
            if not isinstance(hlines_color, (np.ndarray, list, tuple, str, int, float)):
                hlines_color = ['k'] * len(hlines)
            elif isinstance(hlines_color, (str, int, float)):
                hlines_color = [hlines_color] * len(hlines)
            for i in range(len(hlines)):
                plt.axhline(hlines[i], linestyle=hlines_fmt, alpha=hlines_alpha, 
                            label=hlines_label[i], color=hlines_color[i])
            if ylog:
                plt.yscale('log')
        
        plt.grid(grid)
        
        if nr_chunks_order and nr_orders:
            for o in range(nr_orders):
                plt.axvline(nr_chunks_order*o, color='k', alpha=0.5, linestyle=':')
        
        if ((isinstance(scatter_label, list) and scatter_label[0] != None) or
            (isinstance(curve_label, list) and curve_label[0] != None) or
            (isinstance(errorbar_label, list) and errorbar_label[0] != None) or
            (isinstance(hlines_label, list) and hlines_label[0] != None)):
            plt.legend()
        
        if isinstance(yrange, (list,tuple)) and len(yrange) == 2:
            plt.ylim(yrange)
        
        plt.xlabel('Chunk #')
        plt.ylabel(ylabel)
        
        if title == '':
            title = 'Chunk scatter'
        plt.title(title)
        
        if savename is not None:
            fmt = path.splitext(savename)[1][1:]
            plt.savefig(savename, format=fmt, dpi=dpi)
        if show_plot:
            plt.show()
        plt.close()
    
    except Exception as e:
        print(e)


def plot_lsfs_grid(lsf_array, chunks, x_nr=3, y_nr=3, alpha=1.0, xlim=None, 
                   grid=True, savename=None, dpi=300, show_plot=False):
    """Plot a grid of evaluated LSFs
    
    A grid of 'x_nr x y_nr' LSFs is plotted.
    
    Args:
        lsf_array (ndarray[nr_chunks, nr_pix]): An array of evaluated LSFs for
            all chunks of an observation.
        chunks (:class:'ChunkArray'): A list of chunks of an observation.
        x_nr (Optional[int]): Number of LSFs to plot in x-direction
            (along the orders, dispersion direction). Defaults to 3.
        x_nr (Optional[int]): Number of LSFs to plot in y-direction
            (across the orders, cross-dispersion direction). Defaults to 3.
        alpha (Optional[float]): The alpha-value of the curves (should be 
              between 0. and 1.). Defaults to 1.
        xlim (Optional[tuple, list, ndarray]): The x-limits in LSF pixels of 
            each subplot (left and right). If None, the whole evaluated LSF is
            plotted.
        grid (Optional[bool]): Whether or not to plot a grid. Defaults to True.
        savename (Optional[str]): If a pathname is given, the plot is saved
            there. If None, then no saving.
        dpi (Optional[int]): DPI of the saved plot. Defaults to 300.
        show_plot (Optional[bool]): If True, the plot is shown during execution.
            Defaults to False.
    """    
    
    nr_chunks_order = len(chunks[chunks[0].order])
    nr_orders = len(chunks.orders)
    
    fig = plt.figure(figsize=(14,14))
    gs = gridspec.GridSpec(x_nr, y_nr, height_ratios=[1]*y_nr, width_ratios=[1]*x_nr)
    ax = []
    for i in range(y_nr):
        for j in range(x_nr):
            ax.append(plt.subplot(gs[i, j]))
            #ind = int(nr_chunks_total/18.) + int(nr_chunks_total/9.) * (i*3+j)
            ind = (int(nr_orders/y_nr/2.) + int(nr_orders/y_nr) * i) * nr_chunks_order + \
                    (int(nr_chunks_order/x_nr/2.) + int(nr_chunks_order/x_nr) * j)
            ax[-1].plot(lsf_array[ind], '-o', alpha=alpha)
            
            ax[-1].set_title('Chunk {}, order {}, pixel {}'.format(
                    ind, chunks[ind].order, chunks[ind].abspix[int(len(chunks[ind])/2.)]))
            #ax[-1].set_title('Chunk {}'.format(ind))
            ax[-1].set_ylim(np.nanmin(lsf_array)-0.002, np.nanmax(lsf_array)+0.002)
            if isinstance(xlim, (list, tuple, np.ndarray)) and len(xlim) == 2:
                ax[-1].set_xlim(xlim[0], xlim[1])
            plt.grid(grid)
        
    if savename is not None:
        fmt = path.splitext(savename)[1][1:]
        plt.savefig(savename, format=fmt, dpi=dpi)
    if show_plot:
        plt.show()
    plt.close()


