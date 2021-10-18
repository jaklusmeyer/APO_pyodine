import os
import h5py
import numpy as np
from .. import fitters
from ..lib import h5quick



class CombinedResults():
    """Container for timeseries fitting results
    
    This object class is a container for all the individual fit results from
    all observations of a star, and can be used as input in the final velocity 
    weighting to receive RVs. If that has already been performed, it can 
    additionally store the results from the weighting algorithm and RVs.
    
    Args:
        filename (Optional[str, list, tuple]): Either a string to the path from
            which to load a :class:'CombinedResults' object, or a list or tuple
            of pathnames to individual fit results to load.
    """
    
    def __init__(self, filename=None):
        # Check whether one combined result should be loaded,
        # or a list or tuple of individual results
        if isinstance(filename, str):
            self.load_combined(filename)
        elif isinstance(filename, (list,tuple)):
            self.load_individual_results(filename)
        
    
    
    def load_individual_results(self, filenames):
        """Load individual fit results
        
        Args:
            filenames (list or tuple): The pathnames of the files.
        """
        self.nr_files = len(filenames)
        self.res_filenames = [os.path.abspath(f) for f in filenames]
        
        # Get param names and general info from first file
        res_dict = fitters.load_results(filenames[0])
        self.param_names = [k for k in res_dict['params'].keys()]
        self.chunk_names = [k for k in res_dict['chunks'].keys()]
        self.info = {
                'star_name': res_dict['observation']['star_name'].decode(),
                'instrument_name': res_dict['observation']['instrument_name'].decode()
                }
        if 'model' in res_dict.keys() and res_dict['model'] != None:
            self.info['lsf_model'] = res_dict['model']['lsf_model'].decode()
            self.info['stellar_template'] = res_dict['model']['stellar_template'].decode()
            self.info['iodine_file'] = res_dict['model']['iodine_file'].decode()
            self.info['osample_factor'] = res_dict['model']['osample_factor']
        
        self.nr_chunks = len(res_dict['chunks'][self.chunk_names[0]])
        
        # Allocate arrays
        self.observation = {
            #'time_start': np.zeros(nfiles, dtype='S23'),
            'bary_date': np.zeros(self.nr_files),
            'bary_vel_corr': np.zeros(self.nr_files),
            'orig_filename': [''] * self.nr_files
        }
        self.params = {k: np.zeros((self.nr_files, self.nr_chunks)) for k in self.param_names}
        self.errors = {k: np.zeros((self.nr_files, self.nr_chunks)) for k in self.param_names}
        self.chunks = {k: np.zeros((self.nr_files, self.nr_chunks)) for k in self.chunk_names}
        self.redchi2 = np.zeros((self.nr_files, self.nr_chunks))
        self.medcnts = np.zeros((self.nr_files, self.nr_chunks))
        
        # Load the results from all files
        for i, file in enumerate(filenames):
            res_dict = fitters.load_results(file)
            
            for k in self.observation.keys():
                self.observation[k][i] = res_dict['observation'][str(k)]
            for k in self.param_names:
                self.params[k][i] = res_dict['params'][k]
                self.errors[k][i] = res_dict['errors'][k]
            for k in self.chunk_names:
                self.chunks[k][i] = res_dict['chunks'][k]
            self.redchi2[i] = res_dict['redchi2']
            self.medcnts[i] = res_dict['medcounts']
        
        # Initiate an empty tseries property
        self.tseries = {}
        
    def save_combined(self, filename):
        """Save the combined fit results to file
        
        Args:
            filename (str): The filename.
        """
        
        with h5py.File(filename, 'w') as h:
            h5quick.dict_to_group(self.observation, h, 'observation')
            h5quick.dict_to_group(self.params, h, 'params')
            h5quick.dict_to_group(self.errors, h, 'errors')
            h5quick.dict_to_group(self.chunks, h, 'chunks')
            h5quick.dict_to_group(self.info, h, 'info')
            h5quick.dict_to_group(self.tseries, h, 'tseries')
            h['redchi2'] = self.redchi2
            h['medcounts'] = self.medcnts
            h['res_filenames'] = [f.encode('utf8', 'replace') for f in self.res_filenames]
    
    def load_combined(self, filename):
        with h5py.File(filename, 'r') as h:
            self.observation = h5quick.h5data(h['observation'])
            self.params = h5quick.h5data(h['params'])
            self.errors = h5quick.h5data(h['errors'])
            self.chunks = h5quick.h5data(h['chunks'])
            self.info = h5quick.h5data(h['info'])
            try:
                self.tseries = h5quick.h5data(h['tseries'])
            except:
                self.tseries = {}
            self.redchi2 = h5quick.h5data(h['redchi2'])
            self.medcnts = h5quick.h5data(h['medcounts'])
            self.res_filenames = [f.decode() for f in h5quick.h5data(h['res_filenames'])]
            
        self.nr_chunks = self.redchi2.shape[1]
        self.nr_files = len(self.res_filenames)
        self.param_names = [k for k in self.params.keys()]
        self.chunk_names = [k for k in self.chunks.keys()]
        
            