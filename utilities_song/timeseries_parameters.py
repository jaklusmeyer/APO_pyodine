"""
    Set here all the important parameters for the timeseries combination
    of individual chunk results.
    
    Paul Heeren, 8/03/2021
"""

class Timeseries_Parameters:
    
    def __init__(self):        
        
        # If you hand a list of filenames to reject, are these the names of the
        # individual modelling results ('res_files') or of the original 
        # observations ('obs_files')?
        self.reject_type = 'res_files'
        
        # This dictionary defines the parameters used in the weighting
        # algorithm: 
        # - 'good_chunks' and 'good_orders' define which chunks to use in the
        #   computation of observation means (to offset-correct all chunks)
        # - 'sig_limits' defines the lower and upper limit in m/s that chunk
        #   timeseries are allowed to have - outliers are corrected to 
        #   'sig_correct';
        # - 'reweight_alpha', 'reweight_beta' and 'reweight_sigma' are
        #   used in the reweight function to create the chunk weights;
        # - 'weight_correct' is the value that weights of 0 or NaN are 
        #   corrected to.
        self.weighting_pars = {
                'good_chunks': None, #(3, 15), #(150, 350)
                'good_orders': None, #(6,14)
                'sig_limits': [4., 1000.],
                'sig_correct': 1000.,
                'reweight_alpha': 1.8,
                'reweight_beta': 8.0,
                'reweight_sigma': 2.0,
                'weight_correct': 0.01,
                }
        
        # Do chromatic index computation?
        self.do_crx = True
        
        # For writing timeseries results to a text-file:
        self.txt_outkeys = ['bary_date', 'rv', 'rv_err']    # Write these results
        self.txt_delimiter = '\t'                           # Delimiter to use
        self.txt_header = None                              # Header line
        self.txt_outformat = ['%10.5f', '6.4f', '3.4f']     # Output format (make sure
                                                            # this matches the keys!)
        
        # Save the final results to file?
        self.save_comb_res = True
        
        # Create and save analysis plots?
        self.plot_analysis = True
        
        # Print observation names of RV outliers? (Only if plot_analysis is True!)
        self.print_outliers = True
        