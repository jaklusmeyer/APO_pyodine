#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 18:21:08 2021

@author: Paul Heeren
"""

#from os import path, mkdir, remove, listdir
from . import misc


def create_filelist_from_hip(hip, savename='reduce_files.txt', 
                             raw_directory='../data/data_raw', ext_directory='../data/data_ext', 
                             res_directory='', spec_types=None):
    """
    Create a filelist, containing input filenames and desired result filenames
    for the reduction pipeline, from a given HIP number: This routine will
    search through the given directory (default is '~/data/data_raw'), and collect 
    all science spectra names as well as the names of the associated calibration
    specs.
    Use spec_types to chose only certain types (as defined in 
    misc._spec_types).
    Needs to be given as tuple or list.
    """
    # Prepare arrays with desired spec_types
    if spec_types is None:
        cal_spec_types = misc._cal_spec_types
        obj_spec_types = misc._obj_spec_types
    else:
        cal_spec_types = [tp for tp in spec_types if tp in misc._cal_spec_types]
        obj_spec_types = [tp for tp in spec_types if tp in misc._obj_spec_types]
    
    print('Object to search for: ', hip)
    print('Directory to search: ', raw_directory)
    
    hip_files, hip_types = misc.search_fits_by_header('OBJECT', hip, raw_directory)
    
    # Only keep those with type as in obj_spec_types
    hip_files, hip_types = misc.filter_for_type(hip_files, hip_types, obj_spec_types)
    
    print('Found {} science spectra.'.format(len(hip_files)))
    
    # Fetch calibration spectra
    all_files = []
    all_types = []
    cal_nr = 0
    for hip_type, hip_file in zip(hip_types, hip_files):
        cal_files, cal_types = misc.search_calibrations_from_filename(hip_file, cal_spec_types)
        cal_nr += len(cal_files)
        all_files.append([hip_file] + cal_files)
        all_types.append([hip_type] + cal_types)
    
    print('Found {} calibration spectra.'.format(cal_nr))
    
    with open(savename, 'w') as outfile:
        for i in range(len(all_files)):
            outstring = ''
            for file, tp in zip(all_files[i], all_types[i]):
                outstring += file + '\t\t' + tp + '\n'
            outstring += '\n'
            outfile.write(outstring)
    
    print('Results written to {}!'.format(savename))


def create_filelist_from_directory(raw_directory, savename='reduce_files.txt', 
                                   ext_directory=None, res_directory=None, 
                                   spec_types=None):
    """
    Create a filelist, containing input filenames and desired result filenames
    for the reduction pipeline, for exactly one directory: This routine will
    just collect all spectra in this directory and classify them.
    Use spec_types to chose only certain types (as defined in 
    misc._spec_types).
    Needs to be given as tuple or list.
    """
    # Prepare arrays with desired spec_types
    if spec_types is None:
        cal_spec_types = misc._cal_spec_types
        obj_spec_types = misc._obj_spec_types
        spec_types = misc._spec_types
    else:
        cal_spec_types = [tp for tp in spec_types if tp in misc._cal_spec_types]
        obj_spec_types = [tp for tp in spec_types if tp in misc._obj_spec_types]
        spec_types = [tp for tp in spec_types if tp in misc._spec_types]
    
    print('Directory to search: ', raw_directory)
    
    all_files = misc.list_files(raw_directory, constraint='*.fits')
    
    all_types = []
    for file in all_files:
        all_types.append(misc.classify_fits_file(file))
    
    all_files, all_types = misc.filter_for_type(all_files, all_types, spec_types)
    
    print('Found {} spectra.'.format(len(all_files)))
    
    with open(savename, 'w') as outfile:
        outstring = ''
        for file, tp in zip(all_files, all_types):
            outstring += file + '\t\t' + tp + '\n'
        outstring += '\n'
        outfile.write(outstring)
    
    print('Results written to {}!'.format(savename))
        
    
    