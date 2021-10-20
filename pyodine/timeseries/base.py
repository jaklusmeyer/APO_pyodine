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
    
    :param filename: Either a string to the path from which to load an 
        existing :class:`CombinedResults` object, or a list or tuple of 
        pathnames to individual fit results to load.
    :type filename: str, list, tuple, or None
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
        
        :param filenames: The pathnames of the files.
        :type filenames: list or tuple 
        """
        self.nr_files = len(filenames)
        self.res_filenames = [os.path.abspath(f) for f in filenames]
        
        # Check whether the results are saved as 'dill' or 'h5py'. 
        # We assume that all files are of same type,
        # so we should be able to tell from the first one
        
        
        # Get param names and general info from first file. It can be either
        # 'h5py' or 'dill', so take care of that.
        filetype = fitters.filetype_from_ext(filenames[0])
        result = fitters.load_results(filenames[0], filetype=filetype)
        # If it was a 'dill' file, transform the recovered object structure
        # to a dictionary
        if not isinstance(result, dict):
            result = fitters.create_results_dict(result)
        
        self.param_names = [k for k in result['params'].keys()]
        self.chunk_names = [k for k in result['chunks'].keys()]
        self.info = {
                'star_name': result['observation']['star_name'].decode(),
                'instrument_name': result['observation']['instrument_name'].decode()
                }
        if 'model' in result.keys() and result['model'] != None:
            self.info['lsf_model'] = result['model']['lsf_model'].decode()
            self.info['stellar_template'] = result['model']['stellar_template'].decode()
            self.info['iodine_file'] = result['model']['iodine_file'].decode()
            self.info['osample_factor'] = result['model']['osample_factor']
        
        self.nr_chunks = len(result['chunks'][self.chunk_names[0]])
        
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
        self.residuals = np.zeros((self.nr_files, self.nr_chunks))
        self.medcnts = np.zeros((self.nr_files, self.nr_chunks))
        
        # Now load the results from all files and fill up the object properties,
        # again making sure about the file formats
        for i, file in enumerate(filenames):
            filetype = fitters.filetype_from_ext(file)
            result = fitters.load_results(file, filetype=filetype)
            # If it was a 'dill' file, transform the recovered object structure
            # to a dictionary
            if not isinstance(result, dict):
                result = fitters.create_results_dict(result)
            
            for k in self.observation.keys():
                self.observation[k][i] = result['observation'][str(k)]
            for k in self.param_names:
                self.params[k][i] = result['params'][k]
                self.errors[k][i] = result['errors'][k]
            for k in self.chunk_names:
                self.chunks[k][i] = result['chunks'][k]
            self.redchi2[i] = result['redchi2']
            self.residuals[i] = result['residuals']
            self.medcnts[i] = result['medcounts']
        
        # Initiate an empty timeseries property
        self.tseries = {}
    
        
    def save_combined(self, filename):
        """Save the combined fit results to file
        
        :param filename: Save under this filename.
        :type filename: str
        """
        
        # Make sure that the file extension matches the h5py format, and
        # correct if this is not the case
        match, new_filename = fitters.check_filename_format(
                filename, 'h5py', correct=True)
        
        with h5py.File(new_filename, 'w') as h:
            h5quick.dict_to_group(self.observation, h, 'observation')
            h5quick.dict_to_group(self.params, h, 'params')
            h5quick.dict_to_group(self.errors, h, 'errors')
            h5quick.dict_to_group(self.chunks, h, 'chunks')
            h5quick.dict_to_group(self.info, h, 'info')
            h5quick.dict_to_group(self.tseries, h, 'tseries')
            h['redchi2'] = self.redchi2
            h['residuals'] = self.residuals
            h['medcounts'] = self.medcnts
            h['res_filenames'] = [f.encode('utf8', 'replace') for f in self.res_filenames]
    
    
    def load_combined(self, filename):
        """Load a combined fit results object from file
        
        :param filename: The pathname of the file.
        :type filename: str
        """
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
            self.residuals = h5quick.h5data(h['residuals'])
            self.medcnts = h5quick.h5data(h['medcounts'])
            self.res_filenames = [f.decode() for f in h5quick.h5data(h['res_filenames'])]
            
        self.nr_chunks = self.redchi2.shape[1]
        self.nr_files = len(self.res_filenames)
        self.param_names = [k for k in self.params.keys()]
        self.chunk_names = [k for k in self.chunks.keys()]
        
            