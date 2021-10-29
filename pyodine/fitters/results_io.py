import h5py
import numpy as np
from ..lib import h5quick
import dill
import os

_group_keys = ('observation', 'chunks', 'params', 'errors', 'model')
_array_keys = ('reports', 'redchi2', 'residuals', 'medcnts')

_fileformats = {
        'h5py': ('.h5',),
        'dill': ('.pkl',)
        }


def create_results_dict(fit_results):
    """Pack most important fit results of an observation into a dictionary
    
    :param fit_results: The results of an observation.
    :type fit_results: list[:class:`LmfitWrapper.LmfitResult`]
    
    :return: The dictionary with the results.
    :rtype: dict
    """
    
    res_dict = {k: None for k in _group_keys}
    res_dict.update( {k: None for k in _array_keys} )
    
    # Collect observation info
    # Note: Unicode strings are currently not supported by h5py, so we need to
    # convert to bytestring (ascii). Special characters are replaced with
    # question marks
    obs = fit_results[0].chunk.observation
    res_dict['observation'] = {
            'instrument_name': obs.instrument.name.encode('utf8', 'replace'),
            'star_name': obs.star.name.encode('utf8', 'replace'),
            'orig_header': obs.orig_header.tostring(sep='\n').encode('utf8', 'replace'),
            'time_start': obs.time_start.isot.encode('utf8', 'replace'),
            'bary_date': obs.bary_date,
            'bary_vel_corr': obs.bary_vel_corr
            }
    
    # Collect modelling info, and the original filename(s) of the modelled 
    # observation(s)
    res_dict['model'] = {
            'lsf_model': fit_results[0].model.lsf_model.name().encode('utf8', 'replace'),
            'iodine_file': os.path.abspath(fit_results[0].model.iodine_atlas.orig_filename).encode('utf8', 'replace'),
            'osample_factor': fit_results[0].model.osample_factor
            }
    # If not a fit result from O-star modelling
    if fit_results[0].model.stellar_template is not None:
        # Include the template info and the original filename of the observation
        res_dict['model']['stellar_template'] = os.path.abspath(fit_results[0].model.stellar_template.orig_filename).encode('utf8', 'replace')
        res_dict['observation']['orig_filename'] = os.path.abspath(obs.orig_filename).encode('utf8', 'replace')
    else:
        # Include the original filenames of all modelled O-star observations (if more than one)
        if hasattr(obs, 'all_filenames'):
            res_dict['observation']['orig_filename'] = [os.path.abspath(f).encode('utf8', 'replace') for f in obs.all_filenames]
        else:
            res_dict['observation']['orig_filename'] = os.path.abspath(obs.orig_filename).encode('utf8', 'replace')

    # Assume parameter names are the same in all elements of input array
    param_names = list(fit_results[0].params.keys())
    nchunk = len(fit_results)

    # Collect info from all chunks
    res_dict['reports']   = np.array([res.report for res in fit_results], dtype='S')
    res_dict['redchi2']   = np.array([res.redchi for res in fit_results])
    res_dict['residuals'] = np.array([res.residuals for res in fit_results])
    res_dict['medcnts']   = np.array([res.medcnts for res in fit_results])
    
    res_dict['chunks'] = {k: np.zeros(nchunk, dtype='int') for k in ('order', 'firstpix', 'lastpix', 'padding')}
    res_dict['params'] = {k: np.zeros(nchunk, dtype='float64') for k in param_names}
    res_dict['errors'] = {k: np.zeros(nchunk, dtype='float64') for k in param_names}
    
    for i in range(nchunk):
        res = fit_results[i]
        # Get chunk info
        res_dict['chunks']['order'][i]    = res.chunk.order
        res_dict['chunks']['firstpix'][i] = res.chunk.abspix[0]
        res_dict['chunks']['lastpix'][i]  = res.chunk.abspix[-1]
        res_dict['chunks']['padding'][i]  = res.chunk.padding
        # Get parameter values and errors
        for p in param_names:
            res_dict['params'][p][i] = res.params[p]
            res_dict['errors'][p][i] = res.errors[p]
        #res_dict['reports'] += [res.report]
        #res_dict['redchi2'] += [res.redchi]
        #res_dict['residuals'] += [res.rel_res_mean()]
        #res_dict['medcnts'] += [res.medcounts]
    #res_dict['reports'] = np.array(reports, dtype='S')
    #res_dict['redchi2'] = np.array(redchi2)
    #res_dict['residuals'] = np.array(residuals)
    #res_dict['medcnts'] = np.array(medcnts)
    
    return res_dict


def filetype_from_ext(filename):
    """Determine the filetype from the filename extension
    
    :param filename: The pathname of the file.
    :type filename: str
    
    :return: The filetype matching the filename extension.
    :rtype: str
    """
    
    # Split the filename and check the extension
    ext = os.path.splitext(filename)[1]
    
    file_format = None
    for key in _fileformats.keys():
        if ext in _fileformats[key]:
            file_format = key
    
    if isinstance(file_format, str):
        return file_format
    else:
        raise KeyError('File extension {} does not correspond to any known format!'.format(ext))
    
    
def check_filename_format(filename, filetype, correct=True):
    """Check whether the extension of the filename matches the chosen filetype
    
    If the filename does not match and keyword 'correct' is True, a corrected 
    filename is returned. Otherwise just the old one.
    
    :param filename: The chosen filename.
    :type filename: str
    :param filetype: The chosen filetype.
    :type filetype: str
    :param correct: Whether to correct the filename if it does not match the
        filetype. Defaults to True.
    :type correct: bool
    
    :return: Whether the filename matches the filetype or not.
    :rtype: bool
    :return: The filename (corrected if it does not match and 'correct==True').
    :rtype: str
    """
    # Check whether the filetype is known
    if filetype in _fileformats.keys():
        # Split the filename and check the extension
        file_ext = os.path.splitext(filename)
        if file_ext[1] not in _fileformats[filetype]:
            print('The extension {} does not match the filetype {}.'.format(
                    file_ext[1], filetype))
            print('It should be one of: ', _fileformats[filetype])
            
            # Possibly correct
            if correct:
                print('Correcting it to: {}'.format(_fileformats[filetype][0]))
                new_filename = file_ext[0] + _fileformats[filetype][0]
            else:
                new_filename = filename
                
            return False, new_filename
        else:
            return True, filename
    
    else:
        raise KeyError('Filetype {} is not known! Must be either of: '.format(
                filetype), _fileformats.keys())
    

def save_results(filename, fit_results, filetype='h5py'):
    """Preliminary function to save a set of fit results
    
    The results are either saved in 'h5py' format, which is basically a
    dictionary of the most important fit results, or in 'dill' format, where
    the whole object structure is saved and can be recovered later for in-depth 
    analysis. Note that this requires a lot more memory!
    If 'filename' exists, it will be overwritten.
    
    :param filename: Output path for the results file.
    :type filename: str
    :param fit_results: The list of :class:`LmfitWrapper.LmfitResult` objects 
        to save.
    :type fit_results: list[:class:`LmfitWrapper.LmfitResult`]
    :param filetype: In which format should the results be written? Default is 
        'h5py', which saves the most important data in a dictionary format to a 
        compact hdf5-file. If 'dill' is specified instead, the full object 
        structure is saved.
    :type filetype: str
    """
    
    # First check whether the filename matches the chosen type, and correct if
    # it does not
    match, new_filename = check_filename_format(filename, filetype, correct=True)
    
    # If savetype is dill, save the whole object structure
    if filetype == 'dill':
        with open(new_filename, 'wb') as f:
            dill.dump(fit_results, f)
    
    # Else create a results dictionary and save that as hdf5
    elif filetype == 'h5py':
        res_dict = create_results_dict(fit_results)
        
        with h5py.File(new_filename, 'w') as h:
            for key in _group_keys:
                h5quick.dict_to_group(res_dict[key], h, key)
            
            for key in _array_keys:
                #print(key)
                #print(res_dict[key])
                #h5quick.dict_to_group(res_dict[key], h, key)
                h[key] = res_dict[key]
            """
            h5quick.dict_to_group(res_dict['observation'], h, 'observation')
            h5quick.dict_to_group(res_dict['chunks'], h, 'chunks')
            h5quick.dict_to_group(res_dict['params'], h, 'params')
            h5quick.dict_to_group(res_dict['errors'], h, 'errors')
            h['reports'] = reports
            h['redchi2'] = redchi2
            h['residuals'] = residuals
            h['medcounts'] = medcnts
            h5quick.dict_to_group(modinfo, h, 'model')
            """
    
    else:
        raise KeyError('The savetype must be either of "h5py" or "dill",' + \
                       'but {} was supplied!'.format(filetype))
        
        """
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
        residuals = []
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
            residuals += [res.rel_res_mean()]
            medcnts += [res.medcounts]
        reports = np.array(reports, dtype='S')
        redchi2 = np.array(redchi2)
        residuals = np.array(residuals)
        medcnts = np.array(medcnts)
        """
        
    

def load_results(filename, filetype='h5py', force=True):
    """Function to load a set of results
    
    Returns them either as dictionary, if the filetype is 'h5py', or as the
    original object structure if filetype is 'dill'. If 'filetype' does not 
    match the extension of the filename and keyword 'force' is True, then 
    the routine attempts to still load the results as corresponding to the
    filename extension.
    
    :param filename: Path to the file.
    :type filename: str
    :param filetype: If 'h5py', the results are returned as dictionary 
        (default). If 'dill', the original object structure is recovered.
    :type filetype: str
    :param force: Whether or not to force the loading of the results, even if
        the filename does not match the filetype.
    :type force: bool
    
    :return: The fit results.
    :rtype: dict or list[:class:`LmfitResult`]
    """
    
    # First check whether the filename matches the chosen type
    match, filename = check_filename_format(filename, filetype, correct=False)
    
    # If it does not match but loading is forced, attempt to determine the
    # format from the filename
    if not match and force:
        filetype = filetype_from_ext(filename)
        match = True
    
    if match:
        # For dill: recover the whole object structure
        if filetype == 'dill':
            try:
                with open(filename, 'rb') as f:
                    fit_results = dill.load(f)
                return fit_results
            except Exception as e:
                raise(e)
            
        # For h5py: load the results as dictionary
        elif filetype == 'h5py':
            try:
                fit_results = {}
                with h5py.File(filename, 'r') as h:
                    for key in _group_keys + _array_keys:
                        try:
                            fit_results[key] = h5quick.h5data(h[key])
                        except:
                            fit_results[key] = None
                            print('Key {} not in result file!'.format(key))
                return fit_results
            except Exception as e:
                raise(e)        
    
    else:
        raise KeyError('The filename {} does not match the fileformat {}!'.format(
                filename, filetype))
    

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
        h['medcnts'] = medcnts
        h5quick.dict_to_group(modinfo, h, 'model')
        
