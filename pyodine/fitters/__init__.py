import h5py
import numpy as np
from ..lib import h5quick
import dill

from . import lmfit_wrapper
# TODO: This should not be imported by default (maybe in a try/except block?)

_result_keys = ['observation', 'chunks', 'params', 'errors', 'reports', 
                'redchi2', 'medcounts', 'model']

def save_results(filename, fit_results, filetype='h5py'):
    """Preliminary function to save a set of results
    
    Second argument is a list of :class:'LmfitResult' objects.
    If 'filename' exists, it will be overwritten.
    
    Args:
        filename (str): Output path for the results file.
        fit_results (list): A list of :class:'LmfitWrapper.LmfitResult' objects.
        filetype (Optional[str]): In which format should the results be written?
            Default is 'h5py', which saves the most important data in a
            dictionary format to a compact hdf5-file. If 'dill' is specified
            instead, the full object structure is saved and can be recovered
            later for in-depth analysis. Note that this requires a lot more 
            memory!
    """
    
    # If savetype is dill, save the whole object structure
    if filetype == 'dill':
        with open(filename, 'wb') as f:
            dill.dump(fit_results, f)
    
    elif filetype == 'h5py':
        # Collect observation info
        # Note: Unicode strings are currently not supported by h5py, so we need to
        # convert to bytestring (ascii). Special characters are replaced with
        # question marks
        obs = fit_results[0].chunk.observation
        obsinfo = {
            'instrument_name': obs.instrument.name.encode('utf8', 'replace'),
            'star_name': obs.star.name.encode('utf8', 'replace'),
            'orig_header': obs.orig_header.tostring(sep='\n').encode('utf8', 'replace'),
            'time_start': obs.time_start.isot.encode('utf8', 'replace'),
            'bary_date': obs.bary_date,
            'bary_vel_corr': obs.bary_vel_corr,
        }
        
        # Collect modelling info
        modinfo = {
                'lsf_model': fit_results[0].model.lsf_model.name().encode('utf8', 'replace'),
                'iodine_file': fit_results[0].model.iodine_atlas.orig_filename.encode('utf8', 'replace'),
                'osample_factor': fit_results[0].model.osample_factor
                }
        # If not a fit result from O-star modelling
        if fit_results[0].model.stellar_template is not None:
            # Include the template info and the original filename of the observation
            modinfo['stellar_template'] = fit_results[0].model.stellar_template.orig_filename.encode('utf8', 'replace')
            obsinfo['orig_filename'] = obs.orig_filename.encode('utf8', 'replace')
        else:
            # Include the original filenames of all modelled O-star observations (if more than one)
            if hasattr(obs, 'all_filenames'):
                obsinfo['orig_filename'] = [f.encode('utf8', 'replace') for f in obs.all_filenames]
            else:
                obsinfo['orig_filename'] = obs.orig_filename.encode('utf8', 'replace')
    
        # Assume parameter names are the same in all elements of input array
        param_names = list(fit_results[0].params.keys())
        nchunk = len(fit_results)
    
        # Collect info from all chunks
        chunks = {k: np.zeros(nchunk, dtype='int')
                  for k in ('order', 'firstpix', 'lastpix', 'padding')}
        params = {k: np.zeros(nchunk, dtype='float64')
                  for k in param_names}
        errors = {k: np.zeros(nchunk, dtype='float64')
                  for k in param_names}
        reports = []
        redchi2 = []
        medcnts = []
        for i in range(nchunk):
            res = fit_results[i]
            # Get chunk info
            chunks['order'][i] = res.chunk.order
            chunks['firstpix'][i] = res.chunk.abspix[0]
            chunks['lastpix'][i] = res.chunk.abspix[-1]
            chunks['padding'][i] = res.chunk.padding
            # Get parameter values and errors
            for p in param_names:
                params[p][i] = res.params[p]
                errors[p][i] = res.errors[p]
            reports += [res.report]
            redchi2 += [res.redchi]
            medcnts += [res.medcounts]
        reports = np.array(reports, dtype='S')
        redchi2 = np.array(redchi2)
        medcnts = np.array(medcnts)
    
        # Save to HDF5
        with h5py.File(filename, 'w') as h:
            h5quick.dict_to_group(obsinfo, h, 'observation')
            h5quick.dict_to_group(chunks, h, 'chunks')
            h5quick.dict_to_group(params, h, 'params')
            h5quick.dict_to_group(errors, h, 'errors')
            h['reports'] = reports
            h['redchi2'] = redchi2
            h['medcounts'] = medcnts
            h5quick.dict_to_group(modinfo, h, 'model')
    
    else:
        raise KeyError('The savetype must be either of "h5py" or "dill",' + \
                       'but {} was supplied!'.format(filetype))
    

def load_results(filename, filetype='h5py'):
    """Function to load a set of results
    
    Returns them either as dictionary, if the filetype is 'h5py', or as the
    original object structure if filetype is 'dill'.
    
    Args:
        filename (str): Path to the file.
        filetype (Optional[str]): If 'h5py', the results are returned as
            dictionary (default). If 'dill', the original object structure
            is recovered. Note: Obviously the filetype must match the actual
            type of the file to be loaded!
    
    Return:
        dict or list: The fit results.
    """
    
    if filetype == dill:
        try:
            with open(filename, 'rb') as f:
                fit_results = dill.load(f)
            return fit_results
        except Exception as e:
            raise(e)
        
    
    elif filetype == 'h5py':
        try:
            fit_results = {}
            with h5py.File(filename, 'r') as h:
                for key in _result_keys:
                    try:
                        fit_results[key] = h5quick.h5data(h[key])
                    except:
                        fit_results[key] = None
                        print('Key {} not in result file!'.format(key))
            return fit_results
        except Exception as e:
            raise(e)
    
    else:
        raise KeyError('The savetype must be either of "h5py" or "dill",' + \
                       'but {} was supplied!'.format(filetype))
    

def save_template_results(filename, fit_results):
    """Preliminary function to save a set of results from the template creation
    
    This is essentially the same as save_results(), except for the missing
    template path name (obviously). So pull both together?! -> included now
    up there. Let's test it.
    
    If 'filename' exists, it will be overwritten.
    
    Args:
        filename (str): Output path for the results file.
        fit_results (list): A list of :class:'LmfitWrapper.LmfitResult' objects.
    """
    
    # Collect observation info
    # Note: Unicode strings are currently not supported by h5py, so we need to
    # convert to bytestring (ascii). Special characters are replaced with
    # question marks
    obs = fit_results[0].chunk.observation
    obsinfo = {
        'instrument_name': obs.instrument.name.encode('utf8', 'replace'),
        'star_name': obs.star.name.encode('utf8', 'replace'),
        'orig_filename': obs.orig_filename.encode('utf8', 'replace'),
        'orig_header': obs.orig_header.tostring(sep='\n').encode('utf8', 'replace'),
        'time_start': obs.time_start.isot.encode('utf8', 'replace'),
        'bary_date': obs.bary_date,
        'bary_vel_corr': obs.bary_vel_corr,
    }
    
    # Collect modelling info
    modinfo = {
            'lsf_model': fit_results[0].model.lsf_model.name().encode('utf8', 'replace'),
            'iodine_file': fit_results[0].model.iodine_atlas.orig_filename.encode('utf8', 'replace'),
            'osample_factor': fit_results[0].model.osample_factor
            }

    # Assume parameter names are the same in all elements of input array
    param_names = list(fit_results[0].params.keys())
    nchunk = len(fit_results)

    # Collect info from all chunks
    chunks = {k: np.zeros(nchunk, dtype='int')
              for k in ('order', 'firstpix', 'lastpix', 'padding')}
    params = {k: np.zeros(nchunk, dtype='float64')
              for k in param_names}
    errors = {k: np.zeros(nchunk, dtype='float64')
              for k in param_names}
    reports = []
    redchi2 = []
    medcnts = []
    for i in range(nchunk):
        res = fit_results[i]
        # Get chunk info
        chunks['order'][i] = res.chunk.order
        chunks['firstpix'][i] = res.chunk.abspix[0]
        chunks['lastpix'][i] = res.chunk.abspix[-1]
        chunks['padding'][i] = res.chunk.padding
        # Get parameter values and errors
        for p in param_names:
            params[p][i] = res.params[p]
            errors[p][i] = res.errors[p]
        reports += [res.report]
        redchi2 += [res.redchi]
        medcnts += [res.medcounts]
    reports = np.array(reports, dtype='S')
    redchi2 = np.array(redchi2)
    medcnts = np.array(medcnts)

    # Save to HDF5
    with h5py.File(filename, 'w') as h:
        h5quick.dict_to_group(obsinfo, h, 'observation')
        h5quick.dict_to_group(chunks, h, 'chunks')
        h5quick.dict_to_group(params, h, 'params')
        h5quick.dict_to_group(errors, h, 'errors')
        h['reports'] = reports
        h['redchi2'] = redchi2
        h['medcounts'] = medcnts
        h5quick.dict_to_group(modinfo, h, 'model')
        
