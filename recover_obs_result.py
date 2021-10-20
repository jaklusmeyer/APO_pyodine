#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:01:45 2021

@author: paul
"""

import pyodine
import utilities_song as utilities
import numpy as np


def recover_obs_results(filename):
    
    # First check the filetype corresponding to the filename
    filetype = pyodine.fitters.filetype_from_ext(filename)
    result = pyodine.fitters.load_results(filename, filetype)
    
    # If it was a 'dill' file, the result should already be recovered as
    # object structure
    if not isinstance(result, dict):
        chunks = pyodine.components.ChunkArray()
        for r in result:
            chunks.append(r.chunk)
        
        return chunks, result
    
    # For a 'h5py' file in contrast, we need to try and manually rebuild the
    # object structure from the information in the saved dictionary
    else:
        #print(result)
        
        # First the observation name(s)
        if isinstance(result['observation']['orig_filename'], list):
            obs_path = [f.decode() for f in result['observation']['orig_filename']]
        else:
            obs_path = [result['observation']['orig_filename'].decode()]
        
        # Now the stellar template name (if any)
        if 'stellar_template' in result['model'].keys():
            temp_path = result['model']['stellar_template'].decode()
        else:
            temp_path = None
        
        iod_path = result['model']['iodine_file'].decode()
        osample = result['model']['osample_factor']
        lsf_name = result['model']['lsf_model'].decode()
        res_chunks = result['chunks']
        res_params = result['params']
        
        #iod_id = list(utilities.conf.my_iodine_atlases.keys())[list(utilities.conf.my_iodine_atlases.values()).index(iod_path)]
        
        orders = np.unique(np.array([o for o in res_chunks['order']]))
        
        # Load the data
        obs = utilities.load_pyodine.ObservationWrapper(obs_path)
        if temp_path:
            temp = pyodine.template.base.StellarTemplate_Chunked(temp_path)
        else:
            temp = None
        iod = utilities.load_pyodine.IodineAtlas(iod_path)
        
        # Build the model and fitter
        lsf_model = pyodine.models.lsf.model_index[lsf_name]
        wave_model = pyodine.models.wave.LinearWaveModel
        cont_model = pyodine.models.cont.LinearContinuumModel
        model = pyodine.models.spectrum.SimpleModel(
                            lsf_model, wave_model, cont_model, iod, stellar_template=temp, 
                            osample_factor=osample, conv_width=6.)
        
        fitter = pyodine.fitters.lmfit_wrapper.LmfitWrapper(model)
        
        # Build the chunk array
        if temp:
            chunks = pyodine.chunks.wave_defined(obs, temp, 
                                                     width=abs(res_chunks['lastpix'][1]-res_chunks['firstpix'][1])+1,
                                                     orders=orders, padding=res_chunks['padding'][2], order_correction=0)
            print(res_chunks['firstpix'][0])
        else:
            chunks = pyodine.chunks.user_defined(obs, width=abs(res_chunks['lastpix'][1]-res_chunks['firstpix'][1])+1,
                                                     orders=orders, padding=res_chunks['padding'][2],
                                                     chunks_per_order=22, pix_offset0=res_chunks['firstpix'][0])
            print(res_chunks['firstpix'][0])
        nr_chunks_total = len(chunks)
        nr_chunks = len(chunks.get_order(orders[0]))
        nr_orders = len(orders)
        
        print('Total number of created chunks: {} (in result file: {})'.format(nr_chunks_total, len(res_chunks['lastpix'])))
        print('Number of created chunks per order: {}'.format(nr_chunks))
        
        pars = model.guess_params(chunks[0])
        for key in res_params.keys():
            pars[key] = res_params[key][0]
        print(pars)
        
        # Loop over the chunks to build the fit_result object
        fit_result = []
        for i, chunk in enumerate(chunks):
            pars = model.guess_params(chunk)
            for key in res_params.keys():
                pars[key] = res_params[key][i]
            lmfit_pars = fitter.convert_params(pars, to_lmfit=True)
            for key in lmfit_pars.keys():
                lmfit_pars[key].set(vary=False)
            # "Fit"
            fit_result.append(fitter.fit(chunk, lmfit_pars, chunk_ind=i))
        
        return chunks, fit_result