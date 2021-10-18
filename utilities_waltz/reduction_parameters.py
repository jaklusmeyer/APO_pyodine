"""
    Set here all the important parameters for tracing, extracting and
    (possibly) further analyzing the spectra.
    
    Paul Heeren, 3/10/2019
"""
        

class Parameters_Waltz:
    
    def __init__(self):
        
        # General information:
        self.input_directory = '../data/data_raw/u11_x/'        # Directory to work on
        self.avoid_plot      = True        # Whether you want to avoid plots
        self.npools          = 2#4            # Number of pools to use in extraction
        self.output_directory = '../data/data_ext/u11_red2'   # Results directory
        self.order_dir       = 'CERES/wavcals_waltz_01042021'    # ThAr template specs for each order
        self.n_useful        = 46           # Max. nr. of useful orders
        self.ro0             = 79           # Physical order number of lowest extracted order
        #self.rotate          = False        # Whether raw spectra need to be rotated
        #self.flip            = True         # Whether raw spectra need to be flipped in cross-disp. direction
        self.object2do       = 'all'        # What Objects to do ('all' for all in directory, 'new' for not yet reduced ones)
        self.force_pre_process = False       # Forcing to always pre-process spectra
        self.force_flat_extract = False      # Forcing to always extract flats
        self.force_thar_extract = False     # Forcing to always extract thars
        self.force_thar_wavcal = False       # If True, force wavelength calibration even if files already exist
        self.force_sci_extract = False       # Forcing to always extract Science spectra
        self.force_spectral_file_build = True # Forcing to always build the complete wavelength calibrated Science output
        self.bary_wave_scale = False        # If True, the wavelength scale of the Science spectra will already be corrected for barycentric velocity
        self.wavelength_type = None         # Conver wavelengths from 'air_to_vac' or 'vac_to_air' (or None)
        self.output_science  = 'both'       # One of 'SONG', 'CERES', or 'both'
        
        # For scattered light removal:
        self.scatter_type = 'median'        # Either 'median' or 'min' (which values between the orders?)
        self.scatter_neg = False            # Negative values in scattered light allowed?
        self.scatter_option = 1             # If 1, force scatter estimate close to CCD edge
        
        # For order tracing and extraction:
        self.order_multiplicity = 1         # How many fibers/slits are dispersed at once?
        self.order_width     = 8            # Width of the echelle orders in pixels
        self.order_separation = 30          # Separation of the echelle orders in pixels (never used at the moment)
        self.nsigmas         = 4.           # Min. nr. standard dev. of orders above background
        self.min_extract_col = 0            # Min. pixel to start order extraction (in dispersion direction)
        self.max_extract_col = -1           # Max. pixel for order extraction (-1 if no boundary)
        self.trace_degree    = 5            # Nr. coefficients of polynomial to fit the order traces
        self.Marsh_alg       = 0            # Use Marsh (0) or Horne (1) optimal extraction algorithm
        self.NSigma_Marsh    = 5            # Nr. standard dev. to ignore cosmics in weight calculation
        self.NCosmic_Marsh   = 5            # Same as above, but for optimal extraction algorithm
        self.S_Marsh         = 0.4          # Fraction of a pixel considered for the interpolation
        self.N_Marsh         = 3            # Degree of polynomial
        
        # Parameters for wavelength solution:
        self.code_order_corr = 34           # Code order correlating to phys_order_corr (first guess it)
        self.d_order_corr    = 6            # Range around order_corr to find conversion factor
        self.order_conv_fac  = 79           # Conversion factor between code and physical orders
        self.phys_order_corr = 113          # Physical order correlating to code_order_corr
        self.max_pix_shift_disp = 50       # Max. pixel shift of lines in disp. direction
        self.shift_increment = 0.1          # Cross-correlation step in pixels
        self.binning         = 1            # Binning of the ThAr spectrum (line atlas pixels must be corrected)
        
        self.Inverse_m       = False        # Inverse spectrum in cross-dispersion direction
        self.use_cheby       = True         # Use Chebyshev polynomials (recommended)
        self.Dump_Argon      = False        # Dump Argon lines in the analysis?
        self.do_xc           = True         # Do a pre-CCF for each order to better determine pixel shift?!
        self.del_width       = 5.0          # If pre-CCF, reject lines farther away than this
        self.line_width      = 5            # Half-width of zones around lines where the Gaussians are fitted
        self.position_fact   = 1            # Factor used to shift the line pixel positions (should always be 1!)
        self.rms_max         = 150          # Max. rms in m/s of initial wav. solution
        self.minlines_init   = 10           # Minimum lines per order for initial wav. solution
        self.line_sigma      = 2.0          # Approximate sigma of ThAr lines
        self.MRMS            = 150          # Max. rms in m/s of global wav. solution
        self.minlines_global = 500          # Minimum lines for global wav. solution
        self.ncoef_x         = 3
        self.ncoef_m         = 8
        self.npar_wsol = int((min(self.ncoef_x,self.ncoef_m) + 1) * (2*max(self.ncoef_x,self.ncoef_m) - \
                          min(self.ncoef_x,self.ncoef_m) + 2) / 2)
        
        # Parameters for continuum normalization
        self.cont_deg        = 3            # Degree of the fitted polynomial
        self.lower_lim       = 1.           # Lower limit (in terms of sigma) for flux values
        self.upper_lim       = 5.           # Upper limit (in terms of sigma) for flux values
        self.min_frac_bins   = 0.1          # Minimum fraction of original bins to be used
