#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 18:19:21 2021

@author: Paul Heeren
"""

from astropy.io import fits
from astropy.time import Time, TimeDelta
from os import path, listdir, walk
import glob
import numpy as np

# Here the different available spectrum types are defined.
# These arrays are used globally and also when this module is imported elsewhere.
_cal_spec_types = ['i2', 'dark', 'bias', 'thar', 'flat']
_obj_spec_types = ['science']
_spec_types = _cal_spec_types + _obj_spec_types
# Respective differentiating arguments (depending on header keys/values)
# need to be set in function classify_by_header()

# Timezone of observatory (to convert utc to local time)
_timezone = -8


def list_files(directory, constraint='*'):
    """
    Convenience function to list all files in a given directory (sorted by name).
    No sub-directories are shown. If a constraint is given, only files
    full-filling that constraint are returned.
    """
    all_files = sorted(glob.glob(path.join(directory, constraint)))
    return [path.abspath(file) for file in all_files if path.isfile(file)]


def info_from_header(header):
    """
    Deduce information from fits header:
        - type
          (calibration type: i2, thar, flat, bias, dark;
          or science type)
        - date and times (start, midpoint, exposure time)
        - object
    Returns as dictionary.
    """
    # Bias:
    if header['OBJECT'].lower().strip() == 'bias':
        spec_type = 'bias'
    # Dark:
    if header['OBJECT'].lower().strip() == 'dark':
        spec_type = 'dark'
    # Flat:
    elif header['OBJECT'].lower().strip() == 'wideflat':
        spec_type = 'flat'
    # I2:
    elif header['OBJECT'].lower().strip() == 'i2':
        spec_type = 'i2'
    # ThAr:
    elif header['OBJECT'].lower().strip() == 'thar':
        spec_type = 'thar'
    # Science:
    else:
        spec_type = 'science'
    
    try:
        exp_time = float(header['EXPOSURE'])
        if 'millisec' in header.comments['EXPOSURE'].lower():
            exp_time = exp_time / 1000
    except:
        exp_time = float(header['EXPTIME'])
        if 'millisec' in header.comments['EXPTIME'].lower():
            exp_time = exp_time / 1000
    
    # Set up info dict
    info_dict = {
            'spec_type': spec_type,
            'start_date': header['DATE-OBS'].strip(),
            'start_jd': Time(header['DATE-OBS'].strip(), format='isot', scale='utc').jd,
            'mp_mid': header['MP-MID'].strip(),
            'exp_time': exp_time,
            'object': header['OBJECT']
            }
    
    return info_dict


def info_from_fits(filename):
    """
    Convenience function to wrap 'info_from_header' and supply it with the
    header from the fits file. Returns the info dictionary.
    """
    with fits.open(filename) as hdul:
        info_dict = info_from_header(hdul[0].header)
        info_dict['filename'] = filename
    
    return info_dict


def search_fits_by_header(key, value, directory, recursive=True):
    """
    Search through a directory (and its subdirectories if recursive is True)
    and return a list with names of all fits-files containing the specified
    combination of key and value in their headers.
    Also return the type of the spectrum, and the start date (as JD).
    """
    
    out_files = []
    out_infos = []
    
    if recursive is False:
        all_files = list_files(directory, '*.fits')
        for file in all_files:
            with fits.open(file) as hdul:
                if key.upper() in list(hdul[0].header.keys()):
                    if hdul[0].header[key] == value:
                        out_files.append(path.abspath(file))
                        info_dict = info_from_header(hdul[0].header)
                        info_dict['filename'] = file
                        out_infos.append(info_dict)
    
    else:
        for root, dirs, files in walk(directory):
            for file in files:
                if '.fits' in file:
                    with fits.open(path.join(root, file)) as hdul:
                        if key.upper() in list(hdul[0].header.keys()):
                            if hdul[0].header[key] == value:
                                out_files.append(path.abspath(path.join(root, file)))
                                info_dict = info_from_header(hdul[0].header)
                                info_dict['filename'] = file
                                out_infos.append(info_dict)
    # Sort the resulting files
    # No, should already be sorted (from all_files)
    #out_files_sort = sorted(out_files)
    #out_infos_sort = [d for (d,file) in sorted(zip(out_infos, out_files))]
    
    return out_files, out_infos



def search_calibrations_from_filename(filename, spec_types=None):
    """
    Search for all calibration spectra that correspond to the file given.
    They need to be stored in the same directory!
    Returns one list with filenames, one list with infos.
    """
    """
    directory = path.split(filename)[0]
    
    all_files = list_files(directory, '*.fits')
    
    all_infos = []
    for file in all_files:
        tp, dt = type_date_from_fits(file)
        all_types.append(tp)
        all_dates.append(dt)
    
    if spec_types is None:
        spec_types = _cal_spec_types
    
    all_files, all_types, all_dates = filter_for_type(all_files, all_types, all_dates,
                                                      spec_types)
    
    # Sort the resulting files
    all_files_sort = sorted(all_files)
    all_types_sort = [tp for (tp,file) in sorted(zip(all_types, all_files))]
    all_dates_sort = [dt for (dt,file) in sorted(zip(all_dates, all_files))]
    
    return all_files_sort, all_types_sort, all_dates_sort
    """
    pass


def search_calibrations_from_filename_date(filename, date, spec_types=None):
    """
    Search for all calibration spectra that correspond to the file given.
    They need to be stored in the same directory, and their dates need to 
    match the handed date (i.e. they should be from the same night)!
    Returns one list with filenames, one list with types of the spectra, one
    list with dates.
    """
    if spec_types is None:
        spec_types = _cal_spec_types
    
    date_t = Time(date, format='jd', scale='utc')
    
    directory = path.split(filename)[0]
    
    all_files = list_files(directory, '*.fits')
    
    all_types = []
    all_dates = []
    for file in all_files:
        tp, dt = type_date_from_fits(file)
        if tp in spec_types:
            pass
    
    
def filter_for_type(all_files, all_infos, spec_types):
    """
    Convenience function to filter a list of files (and respective infos)
    for the desired types (as defined by spec_types).
    """
    all_files_updated = []
    all_infos_updated = []
    for i, d in enumerate(all_infos):
        if d['spec_type'] in spec_types:
            all_files_updated.append(all_files[i])
            all_infos_updated.append(all_infos[i])
    return all_files_updated, all_infos_updated
    
    
def sort_by_night(all_files, all_infos):
    """
    Sort input files by their respective dates, night by night.
    Lick Observatory is in timezone UTC-8, so we use this to define nights
    to start at 12 (noon) one day and end one day later at 12 (noon), in that timezone.
    """
    night_start = '12:00:00.0'  # Local time when night starts (so early!)
    night_end   = '12:00:00.0'  # Local time when night ends (seems more astronomers-like)
    
    night_dic = {}
    
    for i in range(len(all_files)):
        dt = Time(all_infos[i]['start_jd'], format='jd', scale='utc')
        # Convert to time zone
        dt_tz = dt + TimeDelta(_timezone*3600, format='sec')
        # This is the night_start on that day
        reduced_date = Time(
                dt_tz.isot[:11]+night_start, format='isot', scale='utc')
        
        if dt_tz >= reduced_date:
            # in this case reduced_date is actually the start_date
            # and one day later is the end date
            start_date = reduced_date
            end_date = reduced_date + TimeDelta(1, format='jd')
        else:
            # in this case reduced_date is the end_date
            # and one day earlier is the start date
            end_date = reduced_date
            start_date = reduced_date - TimeDelta(1, format='jd')
        
        # Append local start time to info dictionary
        all_infos[i]['start_local'] = dt_tz.isot
        
        # Now check if these dates are already in the dictionary
        if start_date.isot in night_dic.keys():
            night_dic[start_date.isot].append(all_infos[i])
        else:
            night_dic[start_date.isot] = [all_infos[i]]
    
    return night_dic
    
    
def create_extraction_input_file(savename, night_dic, raw_directory, ext_directory, 
                                 meta_directory):
    """
    Create the input file for the extraction pipeline, holding the spectra to be
    extracted along with the respective calibration spectra for each night, where to 
    save them (ext_directory), and where to save extraction meta data (meta_directory).
    If night_dic, raw_directory, ext_directory and meta_directory are lists or tuples,
    an input file for several directories is created.
    """
    # First check if one or several input directories (and corresponding output
    # directories and night_dics)
    if isinstance(raw_directory, tuple) or isinstance(raw_directory, list):
        nr_dirs = len(raw_directory)
        # The other arguments given need to correspond (be list or tuple, and same length)
        for var, var_str in zip([night_dic, ext_directory, meta_directory],
                                ['night_dic', 'ext_directory', 'meta_directory']):
            if isinstance(var, tuple) or isinstance(var, list):
                if len(var) != nr_dirs:
                    raise TypeError(
                            'Input variable {} not of same length as raw_directory: {}'.format(
                                    var_str, len(var)))
            else:
                raise TypeError('Input variable {} not list or tuple.'.format(var_str))
    elif isinstance(raw_directory, str):
        # Pack into list in order to make handling below possible
        raw_directory = [raw_directory]
        for var, var_str in zip([ext_directory, meta_directory],
                                ['ext_directory', 'meta_directory']):
            if isinstance(var, str):
                var = [var]
            else:
                raise TypeError('Input variable {} is not a string!'.format(var_str))
        if isinstance(night_dic, dict):
            night_dic = [night_dic]
        else:
            raise TypeError('Input variable night_dic is not a dictionary!')
    
    # Get momentary date/time from system (in utc)
    t_now = Time.now()
    
    # Now write to file
    with open(savename, 'w') as f:
        f.write('Extraction input file, created {} (UTC).\n'.format(t_now.isot))
        f.write('File-listings:\n')
        f.write('filename\t\tobject\tspec_type\tstart_date\t\tstart_jd\tmp_mid\t\texp_time\n')
        for i in range(len(raw_directory)):
            f.write('-----------------------------------------------------------\n')
            f.write('raw_directory: {}\n'.format(raw_directory[i]))
            f.write('ext_directory: {}\n'.format(ext_directory[i]))
            f.write('meta_directory: {}\n\n'.format(meta_directory[i]))
            
            for key, value in night_dic[i].items():
                f.write('local start: {}\n'.format(key))
                f.write('(utc: {})\n'.format(
                        Time(key, format='isot', scale='utc') - \
                                TimeDelta(_timezone*3600, format='sec')
                        ))
                
                for j in range(len(value)):
                    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            path.split(value[j]['filename'])[1],
                            value[j]['object'],
                            value[j]['spec_type'],
                            value[j]['start_date'],
                            value[j]['start_jd'],
                            value[j]['mp_mid'],
                            value[j]['exp_time']
                            ))
                f.write('\n')
    print('Done!')
        
    
    
