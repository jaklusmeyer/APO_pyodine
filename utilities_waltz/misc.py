#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 18:19:21 2021

@author: Paul Heeren
"""

from astropy.io import fits
from os import path, listdir, walk
import glob

# Here the different available spectrum types are defined.
# These arrays are used globally and also when this module is imported elsewhere.
_cal_spec_types = ['i2', 'dark', 'bias', 'thar', 'flat']
_obj_spec_types = ['science_i2', 'science_noi2']
_spec_types = _cal_spec_types + _obj_spec_types
# Respective differentiating arguments (depending on header keys/values)
# need to be set in function classify_by_header()


def list_files(directory, constraint='*'):
    """
    Convenience function to list all files in a given directory (sorted by name).
    No sub-directories are shown. If a constraint is given, only files
    full-filling that constraint are returned.
    """
    all_files = sorted(glob.glob(path.join(directory, constraint)))
    return [path.abspath(file) for file in all_files if path.isfile(file)]


def search_fits_by_header(key, value, directory, recursive=True):
    """
    Search through a directory (and its subdirectories if recursive is True)
    and return a list with names of all fits-files containing the specified
    combination of key and value in their headers.
    Also return the type of the spectrum.
    """
    
    out_files = []
    out_types = []
    
    if recursive is False:
        all_files = list_files(directory, '*.fits')
        for file in all_files:
            with fits.open(file) as hdul:
                if key.upper() in list(hdul[0].header.keys()):
                    if hdul[0].header[key] == value:
                        out_files.append(path.abspath(file))
                        out_types.append(classify_by_header(hdul[0].header))
    
    else:
        for root, dirs, files in walk(directory):
            for file in files:
                if '.fits' in file:
                    with fits.open(path.join(root, file)) as hdul:
                        if key.upper() in list(hdul[0].header.keys()):
                            if hdul[0].header[key] == value:
                                out_files.append(path.abspath(path.join(root, file)))
                                out_types.append(classify_by_header(hdul[0].header))
    return out_files, out_types


def classify_by_header(header):
    """
    Classify by given header if the associated spectrum is a calibration spec
    (i2, thar, flat, bias, dark) or science spec (with or without i2). 
    Return identifier.
    """
    # Bias:
    if header['SHUTTER'] == 2 and header['EXPTIME'] <= 1.2E-5:
        spec_type = 'bias'
    # Dark:
    elif header['SHUTTER'] == 2 and header['EXPTIME'] > 1.0E-4:
        spec_type = 'dark'
    # Flat:
    elif header['SHUTTER'] == 0 and header['PCUPSTAT'].lower().strip() == 'out' and \
    header['LEDPOWER'].lower().strip() == 'on' and header['I2CESTAT'].lower().strip() == 'out':
        spec_type = 'flat'
    # I2:
    elif header['SHUTTER'] == 0 and header['PCUPSTAT'].lower().strip() == 'out' and \
    header['LEDPOWER'].lower().strip() == 'on' and header['I2CESTAT'].lower().strip() == 'in':
        spec_type = 'i2'
    # ThAr:
    elif header['SHUTTER'] == 0 and header['PCUPSTAT'].lower().strip() == 'out' and \
    header['THAPOWER'].lower().strip() == 'on':
        spec_type = 'thar'
    # Science with I2:
    elif header['SHUTTER'] == 0 and header['PCUPSTAT'].lower().strip() == 'in' and \
    header['I2CESTAT'].lower().strip() == 'in':
        spec_type = 'science_i2'
    # Science without I2:
    elif header['SHUTTER'] == 0 and header['PCUPSTAT'].lower().strip() == 'in' and \
    header['I2CESTAT'].lower().strip() == 'out':
        spec_type = 'science_noi2'
    # Unidentified
    else:
        spec_type = 'unknown'
    
    return spec_type


def classify_fits_file(filename):
    """
    Convenience function to wrap 'classify_by_header' and supply it with the
    header from the fits file. Returns the type of the fits spectrum.
    """
    with fits.open(filename) as hdul:
        spec_type = classify_by_header(hdul[0].header)
    
    return spec_type
    

def search_calibrations_from_filename(filename, spec_types=None):
    """
    Search for all calibration spectra that correspond to the file given.
    They need to be stored in the same directory!
    Returns one list with filenames, one list with types of the spectra.
    """
    directory = path.split(filename)[0]
    
    all_files = list_files(directory, '*.fits')
    
    all_types = []
    for file in all_files:
        all_types.append(classify_fits_file(file))
    
    if spec_types is None:
        spec_types = _cal_spec_types
    
    all_files, all_types = filter_for_type(all_files, all_types, spec_types)
    
    return all_files, all_types
    
    
def filter_for_type(all_files, all_types, spec_types):
    """
    Convenience function to filter a list of files (and respective types)
    for the desired types (as defined by spec_types).
    """
    all_files_updated = []
    all_types_updated = []
    for i, tp in enumerate(all_types):
        if tp in spec_types:
            all_files_updated.append(all_files[i])
            all_types_updated.append(all_types[i])
    return all_files_updated, all_types_updated
    
    


