from os.path import splitext, abspath

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits as pyfits
from astropy.time import Time, TimeDelta
import h5py
from pyodine import components

from utilities_song_apo import conf

__all__ = ["IodineTemplate", "ObservationWrapper"]


class IodineTemplate(components.IodineAtlas):
    """The iodine template class to be used in the modelling
    :param iodine_cell_id: The iodine cell ID to identify the I2 template
        spectrum by in the :ref:`overview_utilities_conf`, or the direct pathname to the I2
        template spectrum.
    :type iodine_cell_id: int or str
    """
    def __init__(self, iodine_cell):
        if not isinstance(iodine_cell, (int,str)):
            raise KeyError('Argument "iodine_cell" must be either int or string!')
        elif isinstance(iodine_cell, int):
            if iodine_cell in conf.my_iodine_atlases.keys():
                self.orig_filename = conf.my_iodine_atlases[iodine_cell]
            else:
                raise ValueError('Unknown iodine_cell ID!')
        elif isinstance(iodine_cell, str):
            self.orig_filename = iodine_cell
        
        with h5py.File(self.orig_filename, 'r') as h:
            flux = h['flux_normalized'][()]
            wave = h['wavelength_air'][()]    # originally: wavelength_air
        # Ensure increasing wavelength
        if 'APO_I2_cleaned_normalized.h5' in self.orig_filename or 'aposong' in self.orig_filename.lower():
                print('This is APO SONG iodine data ')
        if wave[0] > wave[-1]:
            print('Reversing iodine template wavelength array to ensure increasing order.')
            wave = wave[::-1]
            flux = flux[::-1]
        super().__init__(flux, wave)


class ObservationWrapper(components.Observation):
    """A wrapper for the representation of SONG observation spectra, based
    on the parent class :class:`pyodine.components.Observation`
    
    :param filename: The filename of the observation to load.
    :type filename: str
    :param instrument: The instrument used to obtain the observation. If None,
        the information is drawn from the Fits-header (default).
    :type instrument: :class:`components.Instrument`
    :param star: The star of the observation. If None, the information is 
        drawn from the Fits-header (default).
    :type star: :class:`components.Star`
    """

    # Custom properties
    _spec = None    # Internal storage of spectral flux
    _wave = None    # Internal storage of wavelength solution
    _cont = None    # Internal storage of extracted continuum

    def __init__(self, filename, instrument=None, star=None):
        flux, wave, cont, header = load_file(filename)

        self._flux = flux
        self._wave = wave
        self._cont = cont
        
        """ Weights added. Using this formula from dop code for now
            (the value of 0.008 is the flatfield noise - should be changed)
        if weight is None:# or len(weight) is not len(self.flux):
            #self._weight = (1./self._flux) / (1. + self._flux * 0.008**2)
            self._weight = np.ones(self._flux.shape)
        else:
            self._weight = weight"""

        self.nord = flux.shape[0]
        self.npix = flux.shape[1]

        self.orig_header = header
        self.orig_filename = abspath(filename)

        self.instrument = instrument or get_instrument(header)
        self.star = star or get_star(header)
        self.iodine_in_spectrum, self.iodine_cell_id = check_iodine_cell(header)

        # Camera details
        self.exp_time = get_exposuretime(header, self.instrument)  # or_none(header, 'EXPOSURE')
        self.flux_level = None      # FIXME: Define a flux measure
        self.gain = None            # FIXME: Not in header
        self.readout_noise = None   # FIXME: Not in header
        self.dark_current = None    # FIXME: Not in header

        # Timing
        self.time_start = Time(header['DATE-OBS'].strip(), format='isot', scale='utc')
        self.time_weighted = None

        self.bary_date = get_barytime(header, self.instrument)

        self.bary_vel_corr = or_none(header, 'BVC') * 1000.     # km/s in SONG header
        #self.topo_bary_factor = or_none(header, 'BVCFACT')
        #self.mjd_corr = or_none(header, 'MID-JD')#'MBJD')
        #self.moon_vel = or_none(header, 'MOONVEL') * 1000.  # convert to m/s
        if self.bary_date is None:
            print(f"\033[91mWarning: bary_date is still None for {filename}\033[0m")

        # TODO: Implement flux check
        # TODO: Re-calculate BVC

    def __getitem__(self, order) -> components.Spectrum:
        """Return one or more spectral orders
        
        :param order: The order(s) of the spectrum to return.
        :type order: int, list, ndarray, slice
        
        :return: The desired order(s).
        :rtype: :class:`Spectrum` or list[:class:`Spectrum`]
        """
        # Return one order
        if type(order) is int or hasattr(order, '__int__'):
            flux = self._flux[order]
            wave = self._wave[order]
            cont = self._cont[order]
            #weight = self._weight[order]
            return components.Spectrum(flux, wave=wave, cont=cont)#, weight=weight)
        elif isinstance(order, (list, np.ndarray)):
            return [self.__getitem__(int(i)) for i in order]  # Return MultiOrderSpectrum instead?
        elif type(order) is slice:
            return self.__getitem__([int(i) for i in np.arange(self.nord)[order]])
        else:
            raise IndexError(type(order))

def load_file(filename) -> components.Observation:
    """
    A convenience function to load observation data from file.
    Supports both multi-extension (.fits) and PyODINE-style single-HDU files.

    :param filename: The filename of the observation to load.
    :type filename: str

    :return: flux, wavelength, continuum, header
    """
    from os.path import splitext
    ext = splitext(filename)[1]

    if ext != '.fits':
        raise TypeError(f"Unsupported file format: {ext}")

    try:
        print(f"Trying to load file: {filename}")
        h = pyfits.open(filename)
        header = h[0].header

        # check structure of apo fits file apo song has uncertainities from jons pyvista routine
        if h[0].data is not None and h[0].data.ndim == 3:
            data = h[0].data
            if data.shape[0] < 3:
                raise ValueError("Expected 3 layers (flux, wave, cont) in data cube.")
            # we are removing odd edge columns in APO data
            flux = data[0][:, 12:]
            wave = data[1][:, 12:]
            cont = data[2][:, 12:]
            print(flux.shape, wave.shape, cont.shape, 'shapes of flux wave cont')
            print('min, median, max flux', np.min(flux), np.median(flux), np.max(flux))
            print('range wave', np.ptp(wave))
            print('median cont', np.median(cont))
            
            #chceck for nans in wavelength
            nan_mask = np.isnan(wave)
            print(nan_mask, 'nan mask in wave array!!!')



        # another check if the header has different keywords
        elif 'FLUX' in h and 'WAVE' in h and 'RESPONSE' in h:
            print("Detected multi-extension FITS format")
            flux = h['FLUX'].data[:, 12:]
            wave = h['WAVE'].data[:, 12:]
            cont = h['RESPONSE'].data[:, 12:]

        else:
            raise ValueError("Could not determine FITS file format — missing required HDUs or invalid structure.")

        # Post-processing
        #print if there are negative flux values
        if np.any(flux < 0):
            print(f"Warning: Negative flux values found and set to zero. Count: {(flux < 0).sum()}")
        flux[flux < 0] = 0

        # Ensure increasing wavelength order
        if wave[0, 0] > wave[0, -1]:
            print("Reversing arrays to ensure ascending wavelength")
            wave = wave[..., ::-1]
            flux = flux[..., ::-1]
            cont = cont[..., ::-1]

        # Compute BVC if missing= should be added now in reduction pipeline
        if 'BVC' not in header:
            print('BVC not in header, computing it now')
            bvc = compute_barycentric_velocity(header)
            header['BVC'] = bvc
            print(f"\033[92mComputed BVC: {bvc:.2f} km/s\033[0m")
        else:
            print(f"\033[94mBVC from header: {header['BVC']:.2f} km/s\033[0m")
            
        jd_mid = or_none(header, 'JD-MID')
        bjd_mid = or_none(header, 'BJD-MID')

        def is_valid_jd(jd_val):
            try:
                jd_val = float(jd_val)
                return np.isfinite(jd_val)
            except:
                return False

        # If neither valid JD-MID nor BJD-MID is present, fall back to DATE-OBS
        if not is_valid_jd(jd_mid) and not is_valid_jd(bjd_mid):
            print("No valid JD-MID or BJD-MID found, falling back to DATE-OBS")
            try:
                t_obs = Time(header['DATE-OBS'].strip(), format='isot', scale='utc')
                header['JD-MID'] = t_obs.jd
                print(f"\033[93mFallback JD-MID set to {t_obs.jd:.5f} from DATE-OBS\033[0m")
            except Exception as e:
                print(f"\033[91mFailed to compute fallback JD-MID: {e}\033[0m")

        return flux, wave, cont, header

    except Exception as e:
        print(f"\033[91mError loading file: {e}\033[0m")
        raise
    finally:
        h.close()



def get_star(header) -> components.Star:
    name = or_none(header, 'OBJECT')

    ra_raw = or_none(header, 'OBJ-RA')
    dec_raw = or_none(header, 'OBJ-DEC')

    coordinates = None
    if ra_raw and dec_raw:
        try:
            # Try as HH:MM:SS format
            coordinates = SkyCoord(ra_raw + ' ' + dec_raw, unit=(u.hourangle, u.deg))
        except ValueError:
            try:
                # Try as decimal degrees
                ra_deg = float(ra_raw)
                dec_deg = float(dec_raw)
                coordinates = SkyCoord(ra_deg, dec_deg, unit="deg")
            except Exception as e:
                print(f"Warning: Could not parse coordinates for {name}: {e}")
                coordinates = None

    # Get proper motion
    try:
        proper_motion = (header['S-PM-RA'], header['S-PM-DEC'])
    except Exception:
        proper_motion = (None, None)

    return components.Star(name, coordinates=coordinates, proper_motion=proper_motion)


def get_instrument(header) -> components.Instrument:
    """Determine the instrument from the header and return Instrument object

    :param header: The Fits-header.
    :type header: :class:`fits.Header`
    :return: The instrument object.
    :rtype: :class:`Instrument`
    """
    # Patch missing keywords early (for APO_SONG)
    if 'INSTRUME' not in header:
        header['INSTRUME'] = 'APO_SONG'  # Reasonable default
    if 'PROGRAM' not in header:
        header['PROGRAM'] = 'APO_SONG_REDUCED'

    if 'TELESCOP' in header:
        if 'Node 1' in header['TELESCOP'] and 'Spectrograph' in header['INSTRUM']:
            return conf.my_instruments['song_1']
        elif 'Node 1' in header['TELESCOP'] and 'Spectrograph' in header['INSTRUM']:
            return conf.my_instruments['song_2']
        elif 'Waltz' in header['TELESCOP']:
            return conf.my_instruments['waltz']
        elif 'Hamilton' in header['INSTRUME'] or 'HAMILTON' in header['PROGRAM'].upper() or \
             '3M-COUDE' in header['TELESCOP'].upper() or '3M-CAT' in header['PROGRAM'].upper():
            return conf.my_instruments['lick']
        elif 'APO SONG' in header['INSTRUME'] or 'APO_SONG_REDUCED' in header['PROGRAM']:
            return conf.my_instruments['song_3']

    else:
        if 'NEWCAM' in header['PROGRAM'] and 'hamcat' in header.get('VERSION', ''):
            return conf.my_instruments['lick']

    # If nothing matched, still fail cleanly
    raise TypeError('Could not determine instrument')


def check_iodine_cell(header):
    """Check the position and state of the I2 cell during the observation
    
    :param header: The Fits-header.
    :type header: :class:`fits.Header`
    
    :return: Whether or not the I2 cell was in the light path.
    :rtype: bool
    :return: The ID of the used I2 cell.
    :rtype: int, or None
    """
    # If the IODID keyword is set, we should be safe
    if 'IODID' in header.keys() and header['I2POS'] != 2:
        iodine_in_spectrum = True
        iodine_cell_id = header['IODID']
    # Otherwise, let's make a qualified guess based on the I2POS keyword
    else:
        # TODO: Log this event
        # Position 3 corresponds to id=1
        if header['I2POS'] == 3:
            iodine_in_spectrum = True
            iodine_cell_id = 1
        elif header['I2POS'] == 1:
            iodine_in_spectrum = True
            iodine_cell_id = 2
        else:
            # Position 2 lets the light pass through
            iodine_in_spectrum = False
            iodine_cell_id = None
    return iodine_in_spectrum, iodine_cell_id


def or_none(header, key, fallback_value=None):
    """A convenience function to prevent non-existent Fits-header cards from
    throwing up errors
    
    :param header: The Fits-header.
    :type header: :class:`fits.Header`
    :param key: The keyword of the header card of interest.
    :type key: str
    :param fallback_value: What to return if the header card does not exist
        (default: None).
    :type fallback_value: str, int, float, or None
    
    :return: The header card or the 'fallback_value'.
    :rtype: str, int, float, or None
    """
    try:
        return header[key]
    except KeyError:
        # TODO: Log this event
        return fallback_value


def get_exposuretime(header, instrument):
    """Get the exposure time from the fits header (this extra function is 
    neccessary to make old Lick spectra work smoothly)
    
    """
    if 'SONG' in instrument.name:
        return or_none(header, 'EXPTIME')
    elif 'EXPOSURE' in header and 'Lick' in instrument.name:
        # sometimes it's in millisecs - let's try and catch most of these times
        if header['EXPOSURE'] > 3600.:
            return header['EXPOSURE'] / 1000.
        else:
            return header['EXPOSURE']
    elif 'EXPTIME' in header and 'Lick' in instrument.name:
        return header['EXPTIME']
    else:
        return None

def get_barytime(header, instrument):
    """Get the barycentric timestamp (JD or BJD) from the header"""

    def is_valid_jd(value):
        try:
            value = float(value)
            return np.isfinite(value)
        except Exception:
            return False

    if 'SONG' in instrument.name or 'APO SONG' in instrument.name:
        print('Instrument is SONG/APO SONG!!!')
        # Prefer BJD-MID
        if 'BJD-MID' in header and is_valid_jd(header['BJD-MID']):
            return float(header['BJD-MID'])
        elif 'JD-MID' in header and is_valid_jd(header['JD-MID']):
            return float(header['JD-MID'])
        else:
            print("Missing or invalid BJD-MID / JD-MID in header. Falling back to DATE-OBS.")
            try:
                t_obs = Time(header['DATE-OBS'].strip(), format='isot', scale='utc')
                return t_obs.jd
            except Exception as e:
                print(f"Failed to parse DATE-OBS: {e}")
                return None

    elif 'Lick' in instrument.name:
        try:
            date = header['DATE-OBS'].strip()[:10]
            stime = header['MP-START'].strip()
            mtime = header['MP-MID'].strip()
            time_start = Time(date + 'T' + stime, format='isot', scale='utc')
            bary_date = Time(date + 'T' + mtime, format='isot', scale='utc')
            if time_start > bary_date:
                bary_date += TimeDelta(1., format='jd')
            return bary_date.jd
        except Exception as e:
            print(f"Failed to parse Lick MID-TIME: {e}")
            return None

    print("Instrument not recognized. Could not determine bary time.")
    return None

