"""
    Set here all the important parameters for the I2 reduction pipeline, both
    modeling of individual observations as well as the template creation.
    
    Paul Heeren, 3/02/2021
"""

from pyodine import models

class Parameters:
    
    def __init__(self):
        # General parameters:
        self.osample_obs = 4                # Oversample factor for the observation modeling
        self.ref_spectrum = 'arcturus'      # Reference spectrum to use in normalizer (arcturus or sun)
        self.telluric_mask = 'carmenes'     # Telluric mask to use (carmenes, uves or hitran); 
                                            # if None, tellurics are not taken care of
        self.tell_wave_range = (None,6500)  # Load tellurics only within this wavelength range
        self.tell_dispersion = 0.002        # Dispersion (i.e. wavelength grid) of telluric mask
        self.order_range = (None,None)      # Order range (min,max) to use in observation modeling
        self.chunk_width = 40               # Width of chunks in pixels in observation modeling
        self.chunk_padding = 6              # Padding (left and right) of the chunks in pixels
        self.chunks_per_order = None        # Maximum number of chunks per order (optional)
        self.number_cores = 4               # Number of processor cores for multiprocessing
        
        # Chunked: Normalize chunks in the beginning?
        self.normalize_chunks = False
        
        # Weighting of pixels
        self.bad_pixel_mask = False          # Whether to run the bad pixel mask
        self.bad_pixel_cutoff = 0.22        # Cutoff parameter for the bad pixel mask
        self.correct_obs = False            # Whether to correct the observation in regions of weight = 0.
        self.weight_type = 'flat'           # Type of weights (flat or lick, as implemented in pyodine.components.Spectrum)
        
        # Order used for velocity guess (should be outside I2 region)
        self.velgues_order_range = (2,15)
        
        # I2
        self.i2_to_use = 2                  # Index of I2 FTS to use (see archive/conf.py)
        self.wavelength_scale = 'air'       # Which wavelength scale to use (air or vacuum - should always be the first)
        
        # LSF model
        self.lsf_conv_width = 6.           # LSF is evaluated over this many pixels (times 2)
        
        # Now to the run info: For each modelling run, define a new entry in the following dictionary
        # with all the neccessary information needed
        # (except fitting parameters, those are defined further in constrain_parameters() below)
        self.model_runs = {
                0:
                {# First define the LSF
                 'lsf_model': models.lsf.SingleGaussian,    # LSF model to use (this is absolutely neccessary)
                 # Before the chunks are modeled, you can fit the wavelength guesses from the chunks
                 # over the orders (in order to use the smoothed version as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs?!):
                 'pre_wave_slope_deg': 3,                   # Polynomial degree of dispersion fitting (default: 0, no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept
                 # Fitting keywords
                 'use_chauvenet_pixels': True,             # Chauvenet criterion for pixel outliers? (default: False)
                 # Plotting keywords
                 'plot_success': True,                      # Create plot of fitting success (default: False)
                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (default: False)
                 # Save the result?
                 'save_result': True,                      # Save the result of this run (default: True)
                 # After the chunks have been modeled, you can fit the wavelength results
                 # over the orders (in order to use the smoothed version as input for next run):
                 'wave_slope_deg': 3,                       # Polynomial degree of dispersion fitting (default: 0, no fitting)
                 'wave_intercept_deg': 3,                   # Same as above, for wavelength intercept
                 # Plot some chunks with models and residuals?
                 'plot_chunks': [414],                      # A list with indices of chunks that will be plotted
                 # Calculate, save (and plot) median parameter results (default: False)
                 'median_pars': True,                       #FIXME: Plot?!
                 },
                
                1:
                {# First define the LSF
                 'lsf_model': models.lsf.MultiGaussian_Lick,
                 # Fitting keywords
                 'use_chauvenet_pixels': True,
                 # Plotting keywords
                 'plot_success': True,
                 'plot_analysis': True,
                 # Save the result?
                 'save_result': True,
                 # After the chunks have been modeled, you can fit the wavelength results
                 # over the orders (in order to use the smoothed version as input for next run):
                 'wave_slope_deg': 3,
                 'wave_intercept_deg': 3,
                 # Plot some chunks with models and residuals?
                 'plot_chunks': [414],
                 # Calculate, save (and plot) median parameter results (default: False)
                 'median_pars': True,                       #FIXME: Plot?!
                 }}
        """,
                
                2:
                {# First define the LSF
                 'lsf_model': models.lsf.MultiGaussian_Lick,
                 # Fitting keywords
                 'use_chauvenet_pixels': True,
                 # Plotting keywords
                 'plot_success': True,
                 'plot_analysis': True,
                 # Save the result?
                 'save_result': True,
                 # After the chunks have been modeled, you can fit the wavelength results
                 # over the orders (in order to use the smoothed version as input for next run):
                 'wave_slope_deg': 3,
                 'wave_intercept_deg': 3,
                 # Plot some chunks with models and residuals?
                 'plot_chunks': [414],
                 # Calculate, save (and plot) median parameter results (default: False)
                 'median_pars': True,                       #FIXME: Plot?!
                 }
        }"""
        """,
                
                2:
                {# First define the LSF
                 'lsf_model': models.lsf.FixedLSF,
                 # For fixed lsf, choose run from which to smooth lsf
                 'smooth_lsf_run': 1,                       # If not given, use last run
                 'smooth_pixels': 160,                      # Pixels (in dispersion direction) to smooth over
                 'smooth_orders': 3,                        # Orders (in cross-disp direction) to smooth over
                 'order_separation': 15,                    # Avg. pixels between orders in raw spectrum
                 'smooth_manual_redchi': False,             # If true, calculate smooth weights from manual redchi2
                 'smooth_osample': 0,                       # If larger zero: Oversampling to use in smoothing
                 # Fitting keywords
                 'use_chauvenet_pixels': True,
                 # Plotting keywords
                 'plot_success': True,
                 'plot_analysis': True,
                 # Save the result?
                 'save_result': True,
                 # After the chunks have been modeled, you can fit the wavelength results
                 # over the orders (in order to use the smoothed version as input for next run):
                 'wave_slope_deg': 3,
                 'wave_intercept_deg': 3,
                 # Plot some chunks with models and residuals?
                 'plot_chunks': [414],
                 }
                
        }"""
        

    def constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        """
           Constrain the lmfit_parameters for the models, however you wish!
           Input:
               lmfit_params: The model parameters as lmfit parameter set
               run_id: Which run is this?
               run_results: Dictionary with important observation info and
                            results from previous modelling runs
        """
        ###########################################################################
        # RUN 0
        # Mainly there for first wavelength solution to feed into the next runs.
        # Use a single-Gaussian model.
        ###########################################################################
        if run_id == 0:
            for i in range(len(lmfit_params)):
                # Velocity to velocity guess
                lmfit_params[i]['velocity'].set(
                        value=run_results[run_id]['velocity_guess'])
                
                # SingleGaussian model - just constrain the lsf_fwhm
                lmfit_params[i]['lsf_fwhm'].set(
                        value=2.2, min=0.5, max=4.0)
                
                # Chunked normalized: Fix continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
        
        ###########################################################################
        # RUN 1
        # A full model now, with the full LSF description (constrained?!).
        # Use (smoothed) wavelength results from run 0 as starting values.
        ###########################################################################
        elif run_id == 1:
            # Dictionary of median lsf parameters from previous run
            median_lsf_pars = {p[4:]: run_results[0]['median_pars'][p] for p in run_results[0]['median_pars'] if 'lsf' in p}
            # Fit the lsf from last run to get good starting parameters
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf_pars)
            print(lsf_fit_pars)
            
            for i in range(len(lmfit_params)):
                # Velocity to median velocity from last run
                lmfit_params[i]['velocity'].set(
                        value=run_results[0]['median_pars']['velocity'])#value=run_results[run_id]['velocity_guess'])#
                # Constrain the iodine and template
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[0]['median_pars']['iod_depth'])
                lmfit_params[i]['tem_depth'].set(
                        value=run_results[0]['median_pars']['tem_depth'])
                
                # Wavelength dispersion: use results from before
                if run_results[0]['results'][i].lmfit_result is not None:
                    lmfit_params[i]['wave_slope'].set(
                        value=run_results[0]['wave_slope_fit'][i])#value=run_results[0]['results'][i].params['wave_slope'])##), #,
                        #min=run_results[0]['wave_slope_fit'][i]*0.99,
                        #max=run_results[0]['wave_slope_fit'][i]*1.01)
                    lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['wave_intercept_fit'][i])#run_results[0]['results'][i].params['wave_intercept'])##,
                        #min=run_results[0]['wave_intercept_fit'][i]*0.99,
                        #max=run_results[0]['wave_intercept_fit'][i]*1.01)
                    # Chunked normalized:
                    # Continuum parameters: use results from before (and fix)
                    lmfit_params[i]['cont_intercept'].set(
                        value=run_results[0]['results'][i].params['cont_intercept'])
                    lmfit_params[i]['cont_slope'].set(
                        value=run_results[0]['results'][i].params['cont_slope'])
                    #        value=0., vary=False)
                    #median_lsf_pars = {p[4:]: run_results[0]['results'][i].params[p] \
                    #                   for p in run_results[0]['results'][i].params.keys() if 'lsf' in p}
                
                else:
                    lmfit_params[i]['wave_slope'].set(
                        value=run_results[0]['wave_slope_fit'][i])#), #,
                    lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['wave_intercept_fit'][i])#,
                        #min=run_results[0]['wave_intercept_fit'][i]*0.99,
                        #max=run_results[0]['wave_intercept_fit'][i]*1.01)
                    #median_lsf_pars = {p[4:]: run_results[0]['median_pars'][p] \
                    #                   for p in run_results[0]['median_pars'] if 'lsf' in p}
                    
                #lsf_fit_pars = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf_pars)
                
                # LSF parameters:
                for p in lsf_fit_pars.keys():
                    lmfit_params[i]['lsf_'+p].set(
                        value=lsf_fit_pars[p],
                        min=lsf_fit_pars[p]-0.3,#0.3,
                        max=lsf_fit_pars[p]+0.3)#0.3)
                
                """
                # LSF parameters: limit them
                lsf_parnames = ['lsf_left_1', 'lsf_left_2', 'lsf_left_3', 'lsf_left_4', 'lsf_left_5',
                                'lsf_right1', 'lsf_right2', 'lsf_right3', 'lsf_right4', 'lsf_right5']
                lsf_limits = [[0.,1.4], [0.,1.0], [-0.1,0.8], [-0.2,0.6], [-0.3,0.4],
                              [0.,1.4], [0.,1.0], [-0.1,0.8], [-0.2,0.6], [-0.3,0.4]]
                for j, p in enumerate(lsf_parnames):
                    lmfit_params[i][p].set(
                            min=lsf_limits[j][0],
                            max=lsf_limits[j][1])
                """
                """
                # LSF parameters: limit them
                lmfit_params[i]['lsf_exponent'].set(
                        min=1.5, max=2.5)
                lmfit_params[i]['lsf_sigma'].set(
                        value=run_results[0]['results'][i].params['lsf_fwhm']/2.355,
                        min=run_results[0]['results'][i].params['lsf_fwhm']/2.355 - 0.3,
                        max=run_results[0]['results'][i].params['lsf_fwhm']/2.355 + 0.3)
                lmfit_params[i]['lsf_left'].set(
                        min=-0.5, max=0.5)
                lmfit_params[i]['lsf_right'].set(
                        min=-0.5, max=0.5)
                """
        elif run_id == 2:
            # Dictionary of median lsf parameters from previous run
            median_lsf_pars = {p[4:]: run_results[1]['median_pars'][p] for p in run_results[1]['median_pars'] if 'lsf' in p}
            # Fit the lsf from last run to get good starting parameters
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[1]['lsf_model'], median_lsf_pars)
            print(lsf_fit_pars)
            
            for i in range(len(lmfit_params)):
                # Velocity to median velocity from last run
                lmfit_params[i]['velocity'].set(
                        value=run_results[1]['median_pars']['velocity'])#run_results[run_id]['velocity_guess'])
                # Constrain the iodine and template
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[1]['median_pars']['iod_depth'])
                lmfit_params[i]['tem_depth'].set(
                        value=run_results[1]['median_pars']['tem_depth'])
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[1]['wave_slope_fit'][i])#), #run_results[0]['results'][i].params['wave_slope'],
                        #min=run_results[0]['wave_slope_fit'][i]*0.99,
                        #max=run_results[0]['wave_slope_fit'][i]*1.01)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[1]['wave_intercept_fit'][i])#,
                        #min=run_results[0]['wave_intercept_fit'][i]*0.99,
                        #max=run_results[0]['wave_intercept_fit'][i]*1.01)
                
                # Chunked normalized:
                # Continuum parameters: use results from before (and fix)
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[1]['results'][i].params['cont_intercept'])
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[1]['results'][i].params['cont_slope'])
                #        value=0., vary=False)
                
                # LSF parameters:
                for p in lsf_fit_pars.keys():
                    lmfit_params[i]['lsf_'+p].set(
                        value=lsf_fit_pars[p],
                        min=lsf_fit_pars[p]-0.3,
                        max=lsf_fit_pars[p]+0.3)
        """
        ###########################################################################
        # RUN 2
        # Final run. Use smoothed LSF results from run 1 and keep them fixed.
        # Only vary other parameters (use results from run 1 as starting values).
        ###########################################################################
        elif run_id == 2:
            for i in range(len(lmfit_params)):
                # Velocity to median velocity from last run
                lmfit_params[i]['velocity'].set(
                        value=run_results[1]['median_pars']['velocity'])#run_results[run_id]['velocity_guess'])
                # Constrain the iodine and template
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[1]['median_pars']['iod_depth'])
                lmfit_params[i]['tem_depth'].set(
                        value=run_results[1]['median_pars']['tem_depth'])
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[1]['wave_slope_fit'][i]),
                        #min=run_results[1]['wave_slope_fit'][i]*0.99,
                        #max=run_results[1]['wave_slope_fit'][i]*1.01)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[1]['wave_intercept_fit'][i])
                        #min=run_results[1]['wave_intercept_fit'][i]*0.96,
                        #max=run_results[1]['wave_intercept_fit'][i]*1.04)
                        
                # Chunked normalized:
                # Continuum parameters: use results from before (and fix)
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[1]['results'][i].params['cont_intercept'])
                lmfit_params[i]['cont_slope'].set(
                        value=0., vary=False)
                
                # For fixed, smoothed LSF: Don't vary order and pixel0 values!
                # (and better also not the amplitude)
                lmfit_params[i]['lsf_order'].vary = False
                lmfit_params[i]['lsf_pixel0'].vary = False
                lmfit_params[i]['lsf_amplitude'].vary = False
        """
        return lmfit_params
    
    
    
class Template_Parameters:
    
    def __init__(self):
        # General parameters:
        self.osample_obs = 4                # Oversample factor for the observation modeling
        self.ref_spectrum = 'arcturus'      # Reference spectrum to use in normalizer (arcturus or sun)
        self.telluric_mask = 'carmenes'     # Telluric mask to use (carmenes, uves or hitran); 
                                            # if None, tellurics are not taken care of
        self.tell_wave_range = (None,6500)  # Load tellurics only within this wavelength range
        self.tell_dispersion = 0.002        # Dispersion (i.e. wavelength grid) of telluric mask
        self.chunk_width = 40               # Width of chunks in pixels in observation modeling
        self.chunk_padding = 12             # Padding (left and right) of the chunks in pixels
        self.chunks_per_order = 44          # Maximum number of chunks per order (optional)
        
        # Chunked: Normalize chunks in the beginning?
        self.normalize_chunks = False
        
        # LSF model
        self.lsf_conv_width = 6.           # LSF is evaluated over this many pixels (times 2)
        
        # Weighting of pixels
        self.bad_pixel_mask = True          # Whether to run the bad pixel mask
        self.bad_pixel_cutoff = 0.22        # Cutoff parameter for the bad pixel mask
        self.correct_obs = True             # Whether to correct the observation in regions of weight = 0.
        self.weight_type = 'flat'           # Type of weights (flat or lick, as implemented in pyodine.components.Spectrum)
        
        # Order used for velocity guess (should be outside I2 region)
        self.velgues_order_range = (10,22)
        
        # I2
        self.i2_to_use = 2                  # Index of I2 FTS to use (see archive/conf.py)
        self.wavelength_scale = 'air'       # Which wavelength scale to use (air or vacuum - should always be the first)
        
        # Template creation specific
        self.temp_order_range = (15,30)     # Order range (min,max) to use in template creation
        # Jansson deconvolution LSF smoothing?
        self.jansson_lsf_smoothing = {
                'do_smoothing': False,           # If False, do no LSF smoothing
                'smooth_lsf_run': 1,            # If not given, use last run
                'smooth_pixels': 160,           # Pixels (in dispersion direction) to smooth over
                'smooth_orders': 3,             # Orders (in cross-disp direction) to smooth over
                'order_separation': 15,         # Avg. pixels between orders in raw spectrum
                'smooth_manual_redchi': False,  # If true, calculate smooth weights from manual redchi2
                #'smooth_osample': 0,           # Not needed, same as osample_temp in deconvolution_pars!!!
                }
        self.jansson_run_model = 1              # Model (LSF, wave, cont) from this run used in deconvolution
        self.chunk_weights_redchi = False       # Use fitting red.Chi2 as chunk weights for template? (otherwise analytic)
        # Deconvolution parameters
        self.deconvolution_pars = {
                'osample_temp': 10,             # Oversampling of template
                'jansson_niter': 1200,          # Max. number of iterations in Jansson deconvolution
                'jansson_zerolevel': 0.00,      # Spectrum zero-level in Jansson deconvolution
                'jansson_contlevel': 1.02,      # Spectrum continuum-level in Jansson deconvolution
                'jansson_conver': 0.1,          # Convergence parameter in Jansson deconvolution (careful with that!)
                'jansson_chi_change': 1e-6,     # Minimum change of red.Chi2 in Jansson deconvolution before iterations 
                                                # are stopped
                'lsf_conv_width': self.lsf_conv_width,  # LSF is evaluated over this many pixels (times 2)
                }
        
        # Now to the run info: For each modelling run, define a new entry in the following dictionary
        # with all the neccessary information needed
        # (except fitting parameters, those are defined further in constrain_parameters() below)
        self.model_runs = {
                0:
                {# First define the LSF
                 'lsf_model': models.lsf.SingleGaussian,    # LSF model to use (this is absolutely neccessary)
                 # Before the chunks are modeled, you can fit the wavelength guesses from the chunks
                 # over the orders (in order to use the smoothed version as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs?!):
                 'pre_wave_slope_deg': 3,                   # Polynomial degree of dispersion fitting (default: 0, no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept
                 # Fitting keywords
                 'use_chauvenet_pixels': False,             # Chauvenet criterion for pixel outliers? (default: False)
                 # Plotting keywords
                 'plot_success': True,                      # Create plot of fitting success (default: False)
                 'plot_analysis': True,                     # Create analysis plots (residuals etc.) (default: False)
                 # Save the result?
                 'save_result': True,                      # Save the result of this run (default: True)
                 # After the chunks have been modeled, you can fit the wavelength results
                 # over the orders (in order to use the smoothed version as input for next run):
                 'wave_slope_deg': 3,                       # Polynomial degree of dispersion fitting (default: 0, no fitting)
                 'wave_intercept_deg': 3,                   # Same as above, for wavelength intercept
                 # Plot some chunks with models and residuals?
                 'plot_chunks': [414],                      # A list with indices of chunks that will be plotted
                 # Calculate, save (and plot) median parameter results (default: False)
                 'median_pars': True,                       #FIXME: Plot?!
                 },
                
                1:
                {# First define the LSF
                 'lsf_model': models.lsf.MultiGaussian_Lick,
                 # Before the chunks are modeled, you can fit the wavelength guesses from the chunks
                 # over the orders (in order to use the smoothed version as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs?!):
                 'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (default: 0, no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept
                 # Fitting keywords
                 'use_chauvenet_pixels': True,
                 # Plotting keywords
                 'plot_success': True,
                 'plot_analysis': True,
                 # Save the result?
                 'save_result': True,
                 # After the chunks have been modeled, you can fit the wavelength results
                 # over the orders (in order to use the smoothed version as input for next run):
                 'wave_slope_deg': 3,
                 'wave_intercept_deg': 3,
                 # Plot some chunks with models and residuals?
                 'plot_chunks': [414],
                 # Calculate, save (and plot) median parameter results (default: False)
                 'median_pars': True,                       #FIXME: Plot?!
                 },
                
                2:
                {# First define the LSF
                 'lsf_model': models.lsf.MultiGaussian_Lick,
                 # Before the chunks are modeled, you can fit the wavelength guesses from the chunks
                 # over the orders (in order to use the smoothed version as input for the run)
                 # (probably only makes sense before first run, later use smoothed results from previous runs?!):
                 'pre_wave_slope_deg': 0,                   # Polynomial degree of dispersion fitting (default: 0, no fitting)
                 'pre_wave_intercept_deg': 0,               # Same as above, for wavelength intercept
                 # Fitting keywords
                 'use_chauvenet_pixels': True,
                 # Plotting keywords
                 'plot_success': True,
                 'plot_analysis': True,
                 # Save the result?
                 'save_result': True,
                 # After the chunks have been modeled, you can fit the wavelength results
                 # over the orders (in order to use the smoothed version as input for next run):
                 'wave_slope_deg': 3,
                 'wave_intercept_deg': 3,
                 # Plot some chunks with models and residuals?
                 'plot_chunks': [414],
                 # Calculate, save (and plot) median parameter results (default: False)
                 'median_pars': True,                       #FIXME: Plot?!
                 }
        }
        """,
                
                2:
                {# First define the LSF
                 'lsf_model': models.lsf.FixedLSF,
                 # For fixed lsf, choose run from which to smooth lsf
                 'smooth_lsf_run': 1,                       # If not given, use last run
                 'smooth_pixels': 160,                      # Pixels (in dispersion direction) to smooth over
                 'smooth_orders': 3,                        # Orders (in cross-disp direction) to smooth over
                 'order_separation': 15,                    # Avg. pixels between orders in raw spectrum
                 'smooth_manual_redchi': False,             # If true, calculate smooth weights from manual redchi2
                 'smooth_osample': 0,                       # Oversampling to use in smoothing (if 0: same as models)
                 # Fitting keywords
                 'use_chauvenet_pixels': True,
                 # Plotting keywords
                 'plot_success': True,
                 'plot_analysis': True,
                 # Save the result?
                 'save_result': True,
                 # After the chunks have been modeled, you can fit the wavelength results
                 # over the orders (in order to use the smoothed version as input for next run):
                 'wave_slope_deg': 2,
                 'wave_intercept_deg': 2,
                 # Plot some chunks with models and residuals?
                 'plot_chunks': [414],
                 }
                
        }"""
        

    def temp_constrain_parameters(self, lmfit_params, run_id, run_results, fitter):
        """
           Constrain the lmfit_parameters for the models, however you wish!
           Input:
               lmfit_params: The model parameters as lmfit parameter set
               run_id: Which run is this?
               run_results: Dictionary with important observation info and
                            results from previous modelling runs
        """
        ###########################################################################
        # RUN 0
        # Mainly there for first wavelength solution to feed into the next runs.
        # Use a single-Gaussian model.
        ###########################################################################
        if run_id == 0:
            for i in range(len(lmfit_params)):
                # Fix velocity and template depth (these are 0. and constant 1 then)
                lmfit_params[i]['velocity'].set(
                        vary=False)
                lmfit_params[i]['tem_depth'].set(
                        vary=False)
                # Constrain the iodine to not become negative
                lmfit_params[i]['iod_depth'].set(min=0.1)
                
                # SingleGaussian model - just constrain the lsf_fwhm
                lmfit_params[i]['lsf_fwhm'].set(
                        value=2.2, min=0.5, max=4.0)
                
                # Chunked normalized: Fix continuum slope to 0
                #lmfit_params[i]['cont_slope'].set(
                #        value=0., vary=False)
        
        ###########################################################################
        # RUN 1
        # A full model now, with the full LSF description (constrained?!).
        # Use (smoothed) wavelength results from run 0 as starting values.
        ###########################################################################
        elif run_id == 1:
            # Dictionary of median lsf parameters from previous run
            median_lsf_pars = {p[4:]: run_results[0]['median_pars'][p] for p in run_results[0]['median_pars'] if 'lsf' in p}
            # Fit the lsf from last run to get good starting parameters
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[0]['lsf_model'], median_lsf_pars)
            print(lsf_fit_pars)
            
            for i in range(len(lmfit_params)):
                # Fix velocity and template depth (these are 0. and constant 1 then)
                lmfit_params[i]['velocity'].set(
                        vary=False)
                lmfit_params[i]['tem_depth'].set(
                        vary=False)
                # Constrain the iodine to not become negative
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[0]['median_pars']['iod_depth'])
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[0]['wave_slope_fit'][i])#, #['results'][i].params['wave_slope'],
                        #min=run_results[0]['wave_slope_fit'][i]*0.99, #['results'][i].params['wave_slope']*0.96,
                        #max=run_results[0]['wave_slope_fit'][i]*1.01) #['results'][i].params['wave_slope']*1.04)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['wave_intercept_fit'][i])#, #['results'][i].params['wave_intercept'])
                        #min=run_results[0]['wave_intercept_fit'][i]*0.999999,
                        #max=run_results[0]['wave_intercept_fit'][i]*1.000001)
                    
                # Continuum parameters: use results from before
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[0]['results'][i].params['cont_slope'])
                #       min=run_results[0]['results'][i].params['cont_slope']*0.96,
                #       max=run_results[0]['results'][i].params['cont_slope']*1.04)
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[0]['results'][i].params['cont_intercept'])
                #       min=run_results[0]['results'][i].params['cont_intercept']*0.96,
                #       max=run_results[0]['results'][i].params['cont_intercept']*1.04)
                
                # LSF parameters:
                for p in lsf_fit_pars.keys():
                    lmfit_params[i]['lsf_'+p].set(
                        value=lsf_fit_pars[p],
                        min=lsf_fit_pars[p]-0.3,
                        max=lsf_fit_pars[p]+0.3)
            
        elif run_id == 2:
            # Dictionary of median lsf parameters from previous run
            median_lsf_pars = {p[4:]: run_results[1]['median_pars'][p] for p in run_results[1]['median_pars'] if 'lsf' in p}
            # Fit the lsf from last run to get good starting parameters
            lsf_fit_pars = fitter.fit_lsfs(self.model_runs[1]['lsf_model'], median_lsf_pars)
            print(lsf_fit_pars)
            
            for i in range(len(lmfit_params)):
                # Fix velocity and template depth (these are 0. and constant 1 then)
                lmfit_params[i]['velocity'].set(
                        vary=False)
                lmfit_params[i]['tem_depth'].set(
                        vary=False)
                # Constrain the iodine to not become negative
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[1]['median_pars']['iod_depth'])
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[1]['wave_slope_fit'][i], #['results'][i].params['wave_slope'],
                        min=run_results[0]['wave_slope_fit'][i]*0.99, #['results'][i].params['wave_slope']*0.96,
                        max=run_results[0]['wave_slope_fit'][i]*1.01) #['results'][i].params['wave_slope']*1.04)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[0]['wave_intercept_fit'][i], #['results'][i].params['wave_intercept'])
                        min=run_results[0]['wave_intercept_fit'][i]*0.999999,
                        max=run_results[0]['wave_intercept_fit'][i]*1.000001)
                    
                # Continuum parameters: use results from before
                lmfit_params[i]['cont_slope'].set(
                        value=run_results[0]['results'][i].params['cont_slope'])
                #       min=run_results[0]['results'][i].params['cont_slope']*0.96,
                #       max=run_results[0]['results'][i].params['cont_slope']*1.04)
                lmfit_params[i]['cont_intercept'].set(
                        value=run_results[0]['results'][i].params['cont_intercept'])
                #       min=run_results[0]['results'][i].params['cont_intercept']*0.96,
                #       max=run_results[0]['results'][i].params['cont_intercept']*1.04)
                
                # LSF parameters:
                for p in lsf_fit_pars.keys():
                    lmfit_params[i]['lsf_'+p].set(
                        value=lsf_fit_pars[p],
                        min=lsf_fit_pars[p]-0.3,
                        max=lsf_fit_pars[p]+0.3)
                
                """
                # LSF parameters: limit them
                lsf_parnames = ['lsf_left_1', 'lsf_left_2', 'lsf_left_3', 'lsf_left_4', 'lsf_left_5',
                                'lsf_right1', 'lsf_right2', 'lsf_right3', 'lsf_right4', 'lsf_right5']
                lsf_limits = [[0.,1.4], [0.,1.0], [-0.1,0.8], [-0.2,0.6], [-0.3,0.4],
                              [0.,1.4], [0.,1.0], [-0.1,0.8], [-0.2,0.6], [-0.3,0.4]]
                for j, p in enumerate(lsf_parnames):
                    lmfit_params[i][p].set(
                            min=lsf_limits[j][0],
                            max=lsf_limits[j][1])
                """
                """
                # LSF parameters: limit them
                lmfit_params[i]['lsf_exponent'].set(
                        min=1.5, max=2.5)
                lmfit_params[i]['lsf_sigma'].set(
                        value=run_results[0]['results'][i].params['lsf_fwhm']/2.355,
                        min=run_results[0]['results'][i].params['lsf_fwhm']/2.355 - 0.3,
                        max=run_results[0]['results'][i].params['lsf_fwhm']/2.355 + 0.3)
                lmfit_params[i]['lsf_left'].set(
                        min=-0.5, max=0.5)
                lmfit_params[i]['lsf_right'].set(
                        min=-0.5, max=0.5)
                """
        """
        ###########################################################################
        # RUN 2
        # Final run. Use smoothed LSF results from run 1 and keep them fixed.
        # Only vary other parameters (use results from run 1 as starting values).
        ###########################################################################
        elif run_id == 2:
            for i in range(len(lmfit_params)):
                # Fix velocity and template depth (these are 0. and constant 1 then)
                lmfit_params[i]['velocity'].set(
                        vary=False)
                lmfit_params[i]['tem_depth'].set(
                        vary=False)
                # Constrain the iodine to not become negative
                lmfit_params[i]['iod_depth'].set(
                        value=run_results[1]['median_pars']['iod_depth'])
                
                # Wavelength dispersion: use results from before
                lmfit_params[i]['wave_slope'].set(
                        value=run_results[1]['wave_slope_fit'][i])#, #['results'][i].params['wave_slope'])
                        #min=run_results[1]['wave_slope_fit'][i]*0.99,
                        #max=run_results[1]['wave_slope_fit'][i]*1.01)
                
                # Wavelength intercept: use results from before
                lmfit_params[i]['wave_intercept'].set(
                        value=run_results[1]['wave_intercept_fit'][i])#, #['results'][i].params['wave_intercept'])
                        #min=run_results[1]['wave_intercept_fit'][i]*0.99,
                        #max=run_results[1]['wave_intercept_fit'][i]*1.01)
                     
                # Continuum parameters: use results from before
                lmfit_params[i]['cont_slope'].set(
                       value=run_results[1]['results'][i].params['cont_slope'])
                #       min=run_results[1]['results'][i].params['cont_slope']*0.96,
                #       max=run_results[1]['results'][i].params['cont_slope']*1.04)
                lmfit_params[i]['cont_intercept'].set(
                       value=run_results[1]['results'][i].params['cont_intercept'])
                #       min=run_results[1]['results'][i].params['cont_intercept']*0.96,
                #       max=run_results[1]['results'][i].params['cont_intercept']*1.04)
                
                # For fixed, smoothed LSF: Don't vary order and pixel0 values!
                # (and better also not the amplitude)
                lmfit_params[i]['lsf_order'].vary = False
                lmfit_params[i]['lsf_pixel0'].vary = False
                lmfit_params[i]['lsf_amplitude'].vary = False
        """
        return lmfit_params