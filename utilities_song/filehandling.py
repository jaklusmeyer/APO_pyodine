#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 18:21:08 2021

@author: Paul Heeren
"""

#from os import path, mkdir, remove, listdir
from . import misc


def extract_create_filelist_from_time():


def extract_create_filelist_from_hip(hip, savename='extract_files.txt', 
                             raw_directory='../data/data_raw', ext_directory='../data/data_ext', 
                             meta_directory='../data/meta_ext', spec_types=None):
    """
    Create a filelist, containing input filenames for the reduction pipeline, 
    from a given HIP number. This routine will also search for calibration
    spectra that go with each science spectrum (these calibration spectra
    should sit in the same sub-directory as the science spectrum).
    Input:
        hip: HIP number of star to extract (format: 'HIPxxxxx')
        savename: filename of filelist containing the spectrum names
        raw_directory: directory to search in for spectra
        ext_directory: directory to safe the extracted spectra in
        meta_directory: directory to safe the extraction meta-data in
        spec_types: desired spectrum types for extraction; leave free if all
                    are wished for; otherwise hand list or tuple with desired
                    types as defined in misc._spec_types
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
    
    hip_files, hip_infos = misc.search_fits_by_header('OBJECT', hip, raw_directory)
    
    # Only keep those with type as in obj_spec_types
    hip_files, hip_infos = misc.filter_for_type(hip_files, hip_infos, obj_spec_types)
    
    print('Found {} science spectra.'.format(len(hip_files)))
    
    # Fetch calibration spectra
    all_files = []
    all_infos = []
    cal_nr = 0
    for hip_file, hip_info in zip(hip_files, hip_infos):
        cal_files, cal_infos = misc.search_calibrations_from_filename(hip_file, cal_spec_types)
        cal_nr += len(cal_files)
        all_files.append([hip_file] + cal_files)
        all_infos.append([hip_info] + cal_infos)
    
    print('Found {} calibration spectra.'.format(cal_nr))
    """
    with open(savename, 'w') as outfile:
        for i in range(len(all_files)):
            outstring = ''
            for file, tp in zip(all_files[i], all_types[i]):
                outstring += file + '\t\t' + tp + '\n'
            outstring += '\n'
            outfile.write(outstring)
    
    print('Results written to {}!'.format(savename))
    """


def create_filelist_from_directory(raw_directory, savename='reduce_files.txt', 
                                   ext_directory='../data/data_ext', 
                                   meta_directory='../data/meta_ext', 
                                   spec_types=None):
    """
    Create a filelist, containing input filenames for the reduction pipeline, 
    from all files in a given directory. Spectra will be split up for different
    nights (for optimal extraction).
    Input:
        raw_directory: directory to search in for spectra
        savename: filename of filelist containing the spectrum names
        ext_directory: directory to safe the extracted spectra in
        meta_directory: directory to safe the extraction meta-data in
        spec_types: desired spectrum types for extraction; leave free if all
                    are wished for; otherwise hand list or tuple with desired
                    types as defined in misc._spec_types
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
    # First list all files in the directory
    all_files = misc.list_files(raw_directory, constraint='*.fits')
    
    # Now get information from headers
    all_infos = []
    for file in all_files:
        all_infos.append(misc.info_from_fits(file))
    # Filter out only those types that are required
    all_files, all_infos = misc.filter_for_type(all_files, all_infos, spec_types)
    
    # Sort the files by nights into groups (in a dictionary)
    night_dic = misc.sort_by_night(all_files, all_infos)
    
    print('Found overall {} spectra.'.format(len(all_files)))
    
    # Construct input file for extraction
    misc.create_extraction_input_file(savename, night_dic, raw_directory,
                                      ext_directory, meta_directory)
    """
    with open(savename, 'w') as outfile:
        outstring = ''
        for file, tp in zip(all_files, all_types):
            outstring += file + '\t\t' + tp + '\n'
        outstring += '\n'
        outfile.write(outstring)
    
    print('Results written to {}!'.format(savename))
    """


def read_extraction_filelist(filename):
    """
    Read an extraction input filelist, created with e.g. create_filelist_from_directory.
    """
    with open(filename, 'r') as f:
        # Check if this file is of correct type
        line = f.readline()
        if 'Extraction input file' not in line:
            raise ValueError('File is not an Extraction input file: ', filename)
        # Now loop over all other lines and find out how many input directories there are
        # (plus respective extraction and meta directories), as well as how many nights
        # within each directory there are (search for the timestamps)
        raw_directories = []
        ext_directories = []
        meta_directories = []
        night_dics = []
        for i, line in enumerate(f.readlines()):
            pass
    
    
    