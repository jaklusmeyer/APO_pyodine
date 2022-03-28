"""
Created on Fri Oct  2 15:52:56 2020

Wrapper for the package barycorrpy by Shubham Kanodia and Jason Wright,
a Python version of Jason Eastman and Jason Wright’s IDL code [BaryCorr]
(http://astroutils.astronomy.ohio-state.edu/exofast/pro/exofast/bary/zbarycorr.pro) 
based on [Wright and Eastman (2014)](https://arxiv.org/pdf/1409.4774.pdfby).

This wrapper makes some functions easily accessible for the pyodine package.

@author: Paul Heeren
"""

import logging
import sys
import numpy as np

from barycorrpy.barycorrpy import get_BC_vel
from barycorrpy.barycorrpy import utc_tdb


def bvc_wrapper(bvc_dict, timeseries_dict, use_hip=True, z_meas=None):
    """A simple function to get the barycentric velocities for given 
    observation times, for a star and instrument, all defined in the 
    dictionaries of the CombinedResults object, as well as the correct time 
    (barycentric julian date in the barycentric dynamical time standard).
    If usehip==True: Use the built-in hip catalogue to find star's coords.
    If an array of absolute measured redshifts is handed to z_meas, then the
    precise (multiplicative) algorithm is used (non-predictive).
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    if not isinstance(z_meas, (list, np.ndarray, tuple, float, int)):
        z_meas = 0.0
    
    if use_hip is True and ('star_name' in bvc_dict and 'hip' in 
                            bvc_dict['star_name'].lower()):
        
        hip_nr = hip_from_name(bvc_dict['star_name'])
        
        logging.info('BVC through HIP number: {}'.format(hip_nr))
        
        # Calculate barycentric correction
        bcvel, warn0, stat0 = get_BC_vel(
                JDUTC = timeseries_dict['bary_date'],
                hip_id = hip_nr,
                lat = bvc_dict['instrument_lat'],
                longi = bvc_dict['instrument_long'],
                alt = bvc_dict['instrument_alt'],
                ephemeris = 'de430',
                zmeas = z_meas
                )
        
        # JDUTC to BJDTDB time converter
        bjdtdb, warn1, stat1 = utc_tdb.JDUTC_to_BJDTDB(
                JDUTC = timeseries_dict['bary_date'],
                hip_id = hip_nr,
                lat = bvc_dict['instrument_lat'],
                longi = bvc_dict['instrument_long'],
                alt = bvc_dict['instrument_alt']
                )
    
    else:
        
        ra    = bvc_dict['star_ra']
        dec   = bvc_dict['star_dec']
        pmra  = bvc_dict['star_pmra']
        pmdec = bvc_dict['star_pmdec']
        rv0   = bvc_dict['star_rv0']
        if not np.isfinite(pmra):
            pmra = 0.
        if not np.isfinite(pmdec):
            pmdec = 0.
        
        logging.info('BVC through coordinates:')
        logging.info('RA, DEC:       {}, {} (deg)'.format(ra, dec))
        logging.info('PM_RA, PM_DEC: {}, {} (mas/yr)'.format(pmra, pmdec))
        
        # Calculate barycentric correction
        bcvel, warn0, stat0 = get_BC_vel(
                JDUTC = timeseries_dict['bary_date'],
                ra = ra,
                dec = dec,
                pmra = pmra,
                pmdec = pmdec,
                rv = rv0,
                lat = bvc_dict['instrument_lat'],
                longi = bvc_dict['instrument_long'],
                alt = bvc_dict['instrument_alt'],
                ephemeris = 'de430',
                zmeas = z_meas
                )
        
        # JDUTC to BJDTDB time converter
        bjdtdb, warn1, stat1 = utc_tdb.JDUTC_to_BJDTDB(
                JDUTC = timeseries_dict['bary_date'],
                ra = ra,
                dec = dec,
                pmra = pmra,
                pmdec = pmdec,
                rv = rv0,
                lat = bvc_dict['instrument_lat'],
                longi = bvc_dict['instrument_long'],
                alt = bvc_dict['instrument_alt']
                )
        
    # JDUTC to JDTDB time converter
    #jdtdb = utc_tdb.JDUTC_to_JDTDB(utctime = Obstime,
    #                               leap_update=True)
    
    return bcvel, bjdtdb


def hip_from_name(star_name):
    
    if 'hip' in star_name.lower():
        hip = star_name.lower().replace('hip', '').strip()
        if hip.isdecimal():
            hip = int(hip)
        else:
            hip = None
    else:
        hip = None
    
    return hip
