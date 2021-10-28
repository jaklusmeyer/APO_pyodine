"""
    Set here all the important parameters for the timeseries combination
    of individual chunk results.
    
    Paul Heeren, 8/03/2021
"""

class TimeseriesParameters:
    
    def __init__(self):
        # General parameters:
        self.diagnosis_out = True           # Whether to print and plot diagnosis to file
        
        # Chunk requirements:
        self.min_counts = 1000.             # Minimum counts that a chunk should have!
        self.max_counts = 1000000.           # Maximum counts that a chunk should have!
        self.max_redchi2 = 100000000.           # Maximum red. Chi^2 a chunk should have! #58000
        
        # Second step to reject chunks (first step is with chauvenet criterion):
        # Discard poor red. Chi2 chunks and those with large photon errors
        # -> only use best 2 or 3 sigma (95.5 or 99.7 %)
        self.rejection_percentile = 0.997   # Rejection percentile limit
        
        # Weight chunk series by their relative scatter around observation medians
        self.default_sigma = 1000.          # Default sigma of chunk series with less than 4 good ones ([m/s])
        
        # Finally compute weighted observation velocities
        # -> only use best 2 or 3 sigma (95.5 or 99.7 %)
        self.useage_percentile = 0.997      # Useage percentile limit of chunks
        
        
        
        
        
        self.reject_type = 'res_files'
        self.weighting_pars = {
                'reweight_alpha': 1.8,
                'reweight_beta': 8.0,
                'reweight_sigma': 2.0,
                'weight_correct': 0.01,
                'sig_limits': [4., 1000.],
                'sig_correct': 1000.,
                'good_chunks': None, #(3, 15), #(150, 350)
                'good_orders': None #(6,14)
                }
        self.do_crx = True
        self.txt_outkeys = ['bary_date', 'rv', 'rv_err']
        self.txt_delimiter = '\t'
        self.txt_header = None
        self.txt_outformat = ['%10.5f', '6.4f', '3.4f']
        self.plot_analysis = True
        