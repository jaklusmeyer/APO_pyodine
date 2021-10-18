"""
    Set here all the important parameters for the timeseries combination
    of individual chunk results.
    
    Paul Heeren, 8/03/2021
"""

class Parameters_SONG:
    
    def __init__(self):
        # General parameters:
        self.diagnosis_out = True           # Whether to print and plot diagnosis to file
        
        # Chunk requirements:
        self.min_counts = 1000.             # Minimum counts that a chunk should have!
        self.max_counts = 100000.           # Maximum counts that a chunk should have!
        self.max_redchi2 = 58000.           # Maximum red. Chi^2 a chunk should have!
        
        # Second step to reject chunks (first step is with chauvenet criterion):
        # Discard poor red. Chi2 chunks and those with large photon errors
        # -> only use best 2 or 3 sigma (95.5 or 99.7 %)
        self.rejection_percentile = 0.997   # Rejection percentile limit
        
        # Weight chunk series by their relative scatter around observation medians
        self.default_sigma = 1000.          # Default sigma of chunk series with less than 4 good ones ([m/s])
        
        # Finally compute weighted observation velocities
        # -> only use best 2 or 3 sigma (95.5 or 99.7 %)
        self.useage_percentile = 0.997      # Useage percentile limit of chunks