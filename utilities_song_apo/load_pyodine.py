from os.path import splitext, abspath
import astropy.io.fits
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits as pyfits
from astropy.time import Time, TimeDelta
import h5py
from pyodine import components
from astropy.coordinates import EarthLocation, SkyCoord
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from astropy.modeling import models, fitting
from scipy.interpolate import UnivariateSpline


from utilities_song_apo import conf
# code taken from PyAstronomy to convert vacuum to air wavelengths
import numpy as np

"""def vac_to_air(wvl_angstrom):
    Used in testing APO SONG iodine data
    Convert vacuum wavelength to air wavelength using the Ciddor (1996) formula.
    
    Parameters
    ----------
    wvl_angstrom : float or array
        Wavelength in vacuum [Angstrom]
    
    Returns
    -------
    wvl_air : float or array
        Wavelength in air [Angstrom]
    
    wvl = np.array(wvl_angstrom, dtype=float)
    
    # Ciddor constants (from Applied Optics 35, 1566–1573, 1996)
    k0 = 238.0185
    k1 = 5792105.0
    k2 = 57.362
    k3 = 167917.0

    s2 = (1e4 / wvl) ** 2  # wavenumber squared in um^-2

    # Ciddor refractive index for standard air (CO2 = 450 ppm)
    n = (k1 / (k0 - s2) + k3 / (k2 - s2)) / 1e8 + 1.0

    return wvl / n"""


class IodineTemplate(components.IodineAtlas):
    """The iodine template class to be used in the modelling

    :param iodine_cell_id: The iodine cell ID to identify the I2 template
        spectrum by in the :module:`conf`, or the direct pathname to the I2
        template spectrum.
    :type iodine_cell_id: int or str
    I think apo iodine is 20:50
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
            if 'flux_norm' in h:
                flux = h['flux_norm'][()]
            elif 'flux_normalized' in h:
                flux = h['flux_normalized'][()]
            else:
                raise KeyError("Neither 'flux_norm' nor 'flux_normalized' found in HDF5 file.")
            wave = h['wavelength_air'][()]    # originally: wavelength_air
            #check if this is APO data
            if 'APO_I2_cleaned_normalized.h5' in self.orig_filename or 'aposong' in self.orig_filename.lower():
                print('This is APO SONG iodine data ')
            #if 'APO_I2_from_FITS_cleaned.h5' in self.orig_filename or 'aposong' in self.orig_filename.lower(): no longer needed fixed from frank? 2024-10-15
                #flux = flux / np.max(flux)  # normalize
                #print('This is APO SONG iodine data and its been normalized and converted to air wavelengths')
                
            # Ensure increasing wavelength
            if wave[0] > wave[-1]:
                wave = wave[::-1]
                flux = flux[::-1]
        super().__init__(flux, wave)


class ObservationWrapper(components.Observation):
    """A wrapper for the representation of APO_SONG observation spectra

    :param filename: The filename of the observation to load.
    :type filename: str
    :param instrument: The instrument used to obtain the observation. If None,
        the information is drawn from the Fits-header (default).
    :type instrument: :class:`Instrument`
    :param star: The star of the observation. If None, the information is
        drawn from the Fits-header (default).
    :type star: :class:`Star`
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

         #Weights added. Using this formula from dop code for now
            #(the value of 0.008 is the flatfield noise - should be changed)
        """if weight is None:# or len(weight) is not len(self.flux):
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
        self.iodine_in_spectrum, self.iodine_cell_id = True, 1 #check_iodine_cell(header)

        # Camera details
        self.exp_time = get_exposuretime(header, self.instrument)  # or_none(header, 'EXPOSURE')
        self.flux_level = None      # FIXME: Define a flux measure
        self.gain = None            # FIXME: Not in header
        self.readout_noise = None   # FIXME: Not in header
        self.dark_current = None    # FIXME: Not in header

        # Timing
        try:
            self.date = header['DATE-OBS'].strip()
            self.time_start = Time(self.date, format='isot', scale='utc')
            print(self.time_start, 'time_start')

        except ValueError:  # aka if date is not in isot format (eg. in older fit headers)
            yy = self.date[6:8]
            #print(yy)
            # adding a 20, or 19 to yy, depending on cetury
            if float(yy) <= 50:
                yyyy = '20' + yy
            else:
                yyyy = '19' + yy
                #print(yyyy)
            mm = self.date[3:5]
            dd = self.date[0:2]

            # reshaping old date to isot format by hand (maybe theres a better solution)
            if (self.date[2] == '/') and (self.date[5] == '/') and (float(yy) >= 13):
                self.date = yyyy + '-' + mm + '-' + dd
                self.time_start = Time(self.date, format='isot', scale='utc')


        self.time_weighted = None

        self.bary_date = or_none(header, 'LICKJD')#'MID-JD')
        self.bary_vel_corr = or_none(header, 'LICKBVC')#'BVC')#
        self.bary_date = or_none(header, 'BJD-MID')
        self.bary_vel_corr = or_none(header, 'BVC')
        if self.bary_vel_corr is None:
            try:
                print('Computing barycentric correction--THIS IS APO DATA')
                ra = float(header['RA'])
                dec = float(header['DEC'])
                coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

                if self.time_start is not None:
                    location = EarthLocation(lat=32.7804*u.deg, lon=-105.8203*u.deg, height=2788*u.m)
                    bary_corr = coord.radial_velocity_correction(obstime=self.time_start, location=location)
                    self.bary_vel_corr = bary_corr.to(u.km/u.s).value
                else:
                    print("DATE-OBS is missing or invalid. Can't compute barycentric correction.")
                    self.bary_vel_corr = None

            except Exception as e:
                print(f" Error computing barycentric correction: {e}")
                self.bary_vel_corr = None
        
        
        #self.topo_bary_factor = or_none(header, 'BVCFACT')
        #self.mjd_corr = or_none(header, 'MID-JD')#'MBJD')
        #self.moon_vel = or_none(header, 'MOONVEL') * 1000.  # convert to m/s
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

def compute_barycentric_velocity(header, apo=True):
    try:
        ra = float(header['RA'])
        dec = float(header['DEC'])
        date_obs = header['DATE-OBS']

        time_obs = Time(date_obs, format='isot', scale='utc')
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)

        # APO location from webpage converted to degree Latitude 32° 46' 49" N, Longitude 105° 49' 13" W Elevation 2788 meters
        location = EarthLocation(lat=32.7804*u.deg, lon=-105.8203*u.deg, height=2788*u.m)
        barycorr = coord.radial_velocity_correction(obstime=time_obs, location=location)

        return barycorr.to(u.km/u.s).value  # returns value in m/s

    except Exception as e:
        print(f" Could not compute barycentric velocity: {e}")
        return None
    



from scipy.interpolate import UnivariateSpline

"""def normalize_orders(flux, wave, max_iter=5, sigma_clip_thresh=3.0, spline_s=0.01):

    Robust normalization to provide usable continuum for pyodine.

    Parameters:
    - flux: 2D array (orders x pixels)
    - wave: 2D array (orders x pixels)
    - max_iter: sigma clip iterations
    - sigma_clip_thresh: threshold for sigma clipping
    - spline_s: smoothing factor for UnivariateSpline

    Returns:
    - norm_flux: normalized flux
    - continuum: continuum used for normalization
 
    n_orders, n_pixels = flux.shape
    norm_flux = np.zeros_like(flux)
    continuum = np.ones_like(flux)

    for i in range(n_orders):
        f = np.nan_to_num(flux[i], nan=0.0, posinf=0.0, neginf=0.0)
        w = wave[i]

        # Flip if needed
        if w[0] > w[-1]:
            w = w[::-1]
            f = f[::-1]

        # Masking: finite + increasing wavelength + basic signal threshold
        mask = (f > 0) & np.isfinite(f) & np.isfinite(w) & (np.diff(w, prepend=w[0] - 1e-5) > 0)

        if np.sum(mask) < 10:
            print(f"[Order {i}] Too few valid points, defaulting to flat continuum.")
            cont = np.ones_like(f)
        else:
            try:
                # Sigma clipping
                f_clipped = sigma_clip(f[mask], sigma=sigma_clip_thresh, maxiters=max_iter)
                good = ~f_clipped.mask

                w_fit = w[mask][good]
                f_fit = f[mask][good]

                # Remove duplicates
                w_fit, idx = np.unique(w_fit, return_index=True)
                f_fit = f_fit[idx]

                if len(w_fit) < 10:
                    raise ValueError("Too few points after clipping")

                spline = UnivariateSpline(w_fit, f_fit, s=spline_s * len(f_fit))
                cont = spline(w)
            except Exception as e:
                print(f"[Order {i}] Spline fit failed: {e}. Using flat continuum.")
                cont = np.ones_like(f)

        # Normalize and assign
        cont[cont <= 0] = 1.0  # Avoid divide-by-zero or negative continuum
        norm_flux[i] = f / cont
        continuum[i] = cont

    return norm_flux, continuum"""

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

        # Case 1: PyODINE output format (single HDU, 3D cube)
        if h[0].data is not None and h[0].data.ndim == 3:
            print("Detected PyODINE single-HDU format")
            data = h[0].data
            if data.shape[0] < 3:
                raise ValueError("Expected 3 layers (flux, wave, cont) in data cube.")
            flux = data[0][:, 12:]
            wave = data[1][:, 12:]
            cont = data[2][:, 12:]
            print(flux.shape, wave.shape, cont.shape, 'shapes of flux wave cont')
            print('min, median, max flux', np.min(flux), np.median(flux), np.max(flux))
            print('range wave', np.ptp(wave))
            print('median cont', np.median(cont))

        # Case 2: Original format with named HDUs remove first 12 columns in APO data for weird edge
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

        # Compute BVC if missing
        if 'BVC' not in header:
            print('BVC not in header, computing it now')
            bvc = compute_barycentric_velocity(header)
            header['BVC'] = bvc
            print(f"\033[92mComputed BVC: {bvc:.2f} km/s\033[0m")
        else:
            print(f"\033[94mBVC from header: {header['BVC']:.2f} km/s\033[0m")

        return flux, wave, cont, header

    except Exception as e:
        print(f"\033[91mError loading file: {e}\033[0m")
        raise
    finally:
        h.close()

   
    
    
    
    
def load_file3(filename) -> components.Observation:
    """A convenience function to load observation data from file

    :param filename: The filename of the observation to load.
    :type filename: str

    :return: The flux of the observation spectrum.
    :rtype: ndarray
    :return: The wavelengths of the observation spectrum.
    :rtype: ndarray
    :return: The continuum flux of the observation spectrum.
    :rtype: ndarray
    :return: The Fits-header.
    :rtype: :class:`fits.Header`
    """
    try:
        print(f"Trying to load file: {filename}")

        ext = splitext(filename)[1]
        if ext == '.fits':
            # Load the file
        
            h = astropy.io.fits.open(filename)
            print(h, 'h')
            header = h[0].header
            print(header, 'header'  )
            #print(h[1].data)
            compute_barycentric_velocity(header)

            #flux = h[0].data[0]
            flux = h['FLUX'].data[:,12:]#[0].data
            cont = h['RESPONSE'].data[:,12:]#[:,12:]#[0].data
            wave = h['WAVE'].data[:,12:]#[:,12:]#h[0].data[3]
            # make negative flux values zero
            flux[flux < 0] = 0
            
            
            # attempting to just get flat orders from 0-1 for cross correlation
            #flux, cont = normalize_orders(flux, wave)
            """for i in [10, 20, 30]:
                plt.figure()
                plt.plot(wave[i], flux[i], label='Original')
                plt.plot(wave[i], cont[i], label='Continuum', color='orange')
                plt.plot(wave[i], flux[i]/cont[i], label='Normalized', color='green')
                plt.axhline(1, color='gray', linestyle='--')
                plt.title(f"Order {i}")
                plt.legend()
                plt.show()"""
            #cont_est = median_filter(flux, size=(1, 201))  # per order, 201-pixel smoothing
            #flux = flux/cont_est # attempt to normalize here
            print(header, 'header')
            #check and see since its not in old data
            if 'BVC' not in header:
                print('BVC not in header, computing it now')
                bary_vel_corr_apo = compute_barycentric_velocity(header)
                header['BVC'] = bary_vel_corr_apo
            else:
                print('BVC found in header')
                bvc = header.get('BVC')
                #print this in green
                print(f"\033[92m{bvc}\033[0m", 'bvc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

            bary_vel_corr_apo = compute_barycentric_velocity(header)
            #print this in red
            print(f"\033[91m{bary_vel_corr_apo}\033[0m", 'bary_vel_corr_apo!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
            #header['BVC'] = bary_vel_corr_apo
            #check if the computed BVC makes sense
            #if bvc - bary_vel_corr_apo > 0.1 or bvc - bary_vel_corr_apo < -0.1:
                #print(f" Warning: Header BVC={bvc:.2f} km/s differs from computed BVC={bary_vel_corr_apo:.2f} km/s")
                # check that BVC is somewhere in the plus minus 10 km/s range
            if bary_vel_corr_apo is not None:
                if -10. < bary_vel_corr_apo < 10.:
                    print(f" Barycentric velocity correction computed: {bary_vel_corr_apo:.2f} km/s")
                else:
                    print(f" Computed BVC={bary_vel_corr_apo:.2f} km/s seems off, please check!")
            else:
                print(" Problem with barycentric velocity correction computation.")
    

        
            
            
            #cont = np.ones_like(flux)
            if wave[0, 0] > wave[0, -1]:  # Check first order
                print('Reversing wavelength and flux arrays to be ascending')
                wave = wave[..., ::-1]
                flux = flux[..., ::-1]
                cont = cont[..., ::-1]  # already replaced with ones anyway
            print(f"Loaded FITS file with shape {flux.shape}")
            #print(flux, 'flux')
            print(wave, 'wave')
            #print(cont, 'cont')
            #weight = None
            # making sure teh continuum is just 1s in shape of flux for rightnow
            #cont = np.ones(flux.shape)
            h.close()
            # TODO: Check for `songwriter` signature
            return flux, wave, cont, header
        else:
            # Unsupported file format
            raise TypeError('Unsupported file format (%s)' % ext)
    except IOError:
        print('Could not open file %s' % filename)
    except TypeError as e:
        print('TypeError')
        print(e.args[0])


def get_star(header) -> components.Star:
    """Create a star object based on header data

    :param header: The Fits-header.
    :type header: :class:`fits.Header`

    :return: The star object.
    :rtype: :class:`Star`
    """
    # TODO: Load stars from some kind of catalog based on name instead?

    name = or_none(header, 'OBJECT')
    try:
        coordinates = SkyCoord(
            header['RA'].strip() + ' ' + header['OBJ-DEC'].strip(),
            unit=(u.hourangle, u.deg)
        )
    except:# KeyError:
        # TODO: Log this event
        coordinates = None
    # Get the proper motion vector
    try:
        proper_motion = (header['RA_PM'], header['DEC_PM'])
    except:# KeyError:
        # TODO: Log this event
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
    necessary to make old Lick spectra work smoothly)
    """
    if 'SONG' in instrument.name:
        return or_none(header, 'EXPOSURE') or or_none(header, 'EXPTIME')

    elif 'EXPOSURE' not in header and 'EXPTIME' in header:
        header['EXPOSURE'] = header['EXPTIME']

    elif 'EXPOSURE' in header and 'Lick' in instrument.name:
        # sometimes it's in milliseconds - let's try and catch most of these times
        if header['EXPOSURE'] > 3600.:
            return header['EXPOSURE'] / 1000.
        else:
            return header['EXPOSURE']

    elif 'EXPTIME' in header and 'Lick' in instrument.name:
        return header['EXPTIME']

    return None


def get_barytime(header, instrument):
    """Get the date and time of the weighted midpoint from the fits header
    (this extra function is neccessary to make old Lick spectra work smoothly)\

    """
    if 'SONG' in instrument.name:
        return or_none(header, 'BJD-MID')
    elif 'Lick' in instrument.name:
        # in Lick the MID-time is only given in hrs, mins, secs
        # so we create the MID-JD manually here
        date  = header['DATE-OBS'].strip()[:10]
        stime = header['MP-START'].strip()
        mtime = header['MP-MID'].strip()
        time_start = Time(date+'T'+stime, format='isot', scale='utc')
        bary_date  = Time(date+'T'+mtime, format='isot', scale='utc')
        # check whether time_weighted is on the next full day
        if time_start > bary_date:
            bary_date += TimeDelta(1., format='jd')
        return bary_date.jd
