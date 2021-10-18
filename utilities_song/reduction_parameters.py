"""
    Set here all the important parameters for tracing, extracting and
    (possibly) further analyzing the spectra.
    
    Paul Heeren, 3/10/2019
"""

class Parameters_Lick:
    
    def __init__(self):
        
        # General information:
        self.input_directory = '../data/data_raw/hip96459_iodspec'    # Directory to work on
        self.avoid_plot      = True         # Whether you want to avoid plots
        self.JustExtract     = True         # Whether you just want to extract (or more)
        self.npools          = 4            # Number of pools to use in extraction
        self.output_directory = '../data/data_ext/hip96459_iodspec'   # Results directory
        self.order_dir       = 'CERES/wavcals_lick3'    # ThAr template specs for each order
        self.n_useful        = 38           # Max. nr. of useful orders (39?)
        self.ro0             = 83           # Physical order number of lowest extracted order
        
        # For order tracing and extraction:
        self.order_multiplicity = 1         # How many fibers/slits are dispersed at once?
        self.order_width     = 6            # Width of the echelle orders in pixels
        self.order_separation = 1           # Separation of the echelle orders in pixels
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
        self.code_order_corr = 52           # Code order correlating to phys_order_corr (first guess it)
        self.d_order_corr    = 1            # Range around order_corr to find conversion factor
        self.phys_order_corr = 112          # Physical order correlating to code_order_corr
        self.order_conv_fac  = 60           # Conversion factor between code and physical orders
        self.Inverse_m       = False        # Inverse spectrum in cross-dispersion direction
        self.use_cheby       = True         # Use Chebyshev polynomials (recommended)
        self.Dump_Argon      = False        # Dump Argon lines in the analysis?
        self.rms_max         = 100          # Max. rms in m/s of initial wav. solution
        self.minlines_init   = 10           # Minimum lines for initial wav. solution
        self.line_sigma      = 6.2          # Approximate sigma of ThAr lines
        self.MRMS            = 1000         # Max. rms in m/s of global wav. solution
        self.minlines_global = 400          # Minimum lines for global wav. solution
        self.ncoef_x         = 2
        self.ncoef_m         = 5
        self.npar_wsol = int((min(self.ncoef_x,self.ncoef_m) + 1) * (2*max(self.ncoef_x,self.ncoef_m) - \
                          min(self.ncoef_x,self.ncoef_m) + 2) / 2)
        
        # self.models_path = base+"data/COELHO_MODELS/R_40000b/"  # for atmospheric modelling
        