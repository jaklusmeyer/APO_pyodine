from os.path import splitext, abspath#, isdir

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits as pyfits
import h5py
from astropy.time import Time, TimeDelta
import glob
from CERES import components

from utilities_lick import conf


class LickScienceSpec(components.ScienceSpec):

    # Custom properties
    _spec = None    # Internal storage of spectral flux

    def __init__(self, filename, instrument=None, star=None):
        flux, header = load_raw_file(filename)

        super().__init__(flux)
        self.orig_header = header
        self.orig_filename = abspath(filename)

        self.instrument = instrument or get_instrument(header)
        self.star = star or get_star(header)
        self.iodine_in_spectrum = True # for Lick science-specs always true
        #self.iodine_in_spectrum, self.iodine_cell_id = check_iodine_cell(header)

        # Camera details
        self.exp_time = or_none(header, 'EXPOSURE')
        self.flux_level = None      # FIXME: Define a flux measure
        self.gain = None            # FIXME: Not in header
        self.readout_noise = None   # FIXME: Not in header
        self.dark_current = None    # FIXME: Not in header

        # Timing
        self.date, self.time_start, self.time_weighted = times_from_header(header)
        
        # TODO: bary_date_weighted and deal with times above
        # TODO: Implement flux check
        # TODO: Re-calculate BVC


class CalibrationSpec(components.ScienceSpec):
    
    def __init__(self, filename=None, instrument=None, flux=None, header=None):
        if type(flux) is not np.ndarray:
            flux, header = load_raw_file(filename)
        super().__init__(flux)
        self.orig_header = header
        self.orig_filename = None if filename is None else abspath(filename)

        self.instrument = instrument or get_instrument(header)

        # Camera details
        self.flux_level = None      # FIXME: Define a flux measure
        self.gain = None            # FIXME: Not in header
        self.readout_noise = None   # FIXME: Not in header
        self.dark_current = None    # FIXME: Not in header        


class NightLog(components.Log_archive):
    """
        A nightly log of all spectra to be processed
    """
    def __init__(self, directoryname=None):
        self.directoryname = directoryname
        self.images = {
                'biases': [],
                'darks': [],
                'flats': [],
                'thars': [],
                'i2s': [],
                'sciens': [],
                }
        self.infos = {
                }
        for key in self.images.keys():
            self.infos[key+'_mjd'] = []
        #self.biases = {}
        #self.darks  = {}
        #self.flats  = {}
        #self.thars  = {}
        #self.i2s    = {}
        #self.sciens = {}
        #self.imagetypes = ['biases', 'darks', 'flats', 'thars', \
        #                   'i2s', 'sciens']
    
    def classify_files(self, constraint='*fits'):
        """
            Classifies all files in a directory
            and writes a night log of science images
        """
        all_files = sorted(glob.glob(self.directoryname + '/*fits'))
        for filename in all_files:
            with pyfits.open(filename) as h:
                header = h[0].header
                date, time_start, time_weighted = times_from_header(header)
                if header['OBJECT'].lower().strip() == 'bias':
                    self.images['biases'].append(filename)
                    self.infos['biases_mjd'].append(time_weighted.mjd)
                if header['OBJECT'].lower().strip() == 'dark':
                    self.images['darks'].append(filename)
                    self.infos['darks_mjd'].append(time_weighted.mjd)
                elif header['OBJECT'].lower().strip() == 'wideflat':
                    self.images['flats'].append(filename)
                    self.infos['flats_mjd'].append(time_weighted.mjd)
                elif header['OBJECT'].lower().strip() == 'i2':
                    self.images['i2s'].append(filename)
                    self.infos['i2s_mjd'].append(time_weighted.mjd)
                elif header['OBJECT'].lower().strip() == 'thar':
                    self.images['thars'].append(filename)
                    self.infos['thars_mjd'].append(time_weighted.mjd)
                else:
                    self.images['sciens'].append(filename)
                    self.infos['sciens_mjd'].append(time_weighted.mjd)
                    #self.sciens[filename]['object'] = header['OBJECT'].upper()
                    #self.sciens[filename]['exp_time'] = float(header['EXPTIME'])
        """ OLD:
        for filename in all_files:
            with pyfits.open(filename) as h:
                header = h[0].header
                date, time_start, time_weighted = times_from_header(header)
                if header['OBJECT'].lower().strip() == 'bias':
                    self.biases[filename] = {}
                    self.biases[filename]['mjd'] = time_weighted.mjd
                elif header['OBJECT'].lower().strip() == 'wideflat':
                    self.flats[filename] = {}
                    self.flats[filename]['mjd'] = time_weighted.mjd
                elif header['OBJECT'].lower().strip() == 'thar':
                    self.thars[filename] = {}
                    self.thars[filename]['mjd'] = time_weighted.mjd
                elif header['OBJECT'].lower().strip() == 'i2':
                    self.i2s[filename] = {}
                    self.i2s[filename]['mjd'] = time_weighted.mjd
                else:
                    self.sciens[filename] = {}
                    self.sciens[filename]['object'] = header['OBJECT'].upper()
                    self.sciens[filename]['mjd_start'] = time_start.mjd
                    self.sciens[filename]['mjd_weighted'] = time_weighted.mjd
                    self.sciens[filename]['exp_time'] = float(header['EXPOSURE'])
        """
    
    def __getitem__(self, *args):
        typename, = args
        if isinstance(typename, str):
            if typename in self.images.keys():
                return self.images[typename]
            elif typename == 'all':
                allfiles = []
                for key in self.images.keys():
                    allfiles += self.images[key]
                return allfiles
            else:
                raise ValueError(
                        'Imagetype not existing. Choose one of {}'.format(
                        [key for key in self.images.keys()])
                        )
        else:
            raise ValueError(
                    'Please enter a string. Choose one of {}'.format(
                    [key for key in self.images.keys()])
                    )        
    
    def __len__(self):
        length = 0
        for key in self.images.keys():
            length += len(self.images[key])
        return length
    
    def __str__(self):
        outstring = '<NightLog (sciens:'
        for spec in self.images['sciens']:
            outstring += '\n\t{},'.format(spec)
        outstring += '\n  '
        for key in self.images.keys():
            if key != 'sciens':
                outstring += '{}  '.format(key)
        outstring += '\n  '
        
        outstring += '  {}       {}      {}      {}     {}   )>'.format(
                len(self.images['biases']), len(self.images['darks']),
                len(self.images['flats']), len(self.images['thars']),
                len(self.images['i2s']))
        return outstring


def load_raw_file(filename) -> components.ScienceSpec:
    try:
        ext = splitext(filename)[1]
        if ext == '.fits':
            # Load the fits file
            with pyfits.open(filename) as h:
                flux = h[0].data
                header = {}
                for key, value in h[0].header.items():
                    if key != '':
                        header[key] = value
            return flux, header
        elif ext == '.h5':
            # Load the hdf5 file
            with h5py.File(filename, 'r') as f:
                flux = f['data'][:]
                header = {}
                for key, value in f['data'].attrs.items():
                    header[key] = value
            return flux, header
        else:
            # Unsupported file format
            raise TypeError('Unsupported file format (%s)' % ext)
    except IOError:
        print('Could not open file %s' % filename)
    except TypeError as e:
        print(e.args[0])


def get_star(header) -> components.Star:
    """
        Create a star object based on header data
    """
    # TODO: Load stars from some kind of catalog based on name instead?
    name = or_none(header, 'OBJECT')
    try:
        coordinates = SkyCoord(
            header['RA'] + ' ' + header['DEC'],
            unit=(u.hourangle, u.deg)
        )
    except KeyError:
        # TODO: Log this event
        coordinates = None
    # Get the proper motion vector
    try:
        proper_motion = (header['RA_PM'], header['DEC_PM'])
    except KeyError:
        # TODO: Log this event
        proper_motion = (None, None)
    return components.Star(name, coordinates=coordinates, proper_motion=proper_motion)


def get_instrument(header) -> components.Instrument:
    """
        Determine the instrument from the header and return Instrument object
    """
    if 'INSTRUME' in header.keys():
        if 'Hamilton Spec.' in header['INSTRUME'] or \
        'Hamilton Spectrograph' in header['INSTRUME']:
            return conf.my_instruments['lick']
    if 'VERSION' in header.keys():
        if 'hamcat' in header['VERSION']:
            return conf.my_instruments['lick']
    if 'TELESCOP' in header.keys():
        if 'Waltz Telescope' in header['TELESCOP']:
            return conf.my_instruments['waltz']
    else:
        # TODO: Log this event
        raise TypeError('Could not determine instrument')


def or_none(header, key, fallback_value=None):
    try:
        return header[key]
    except KeyError:
        # TODO: Log this event
        return fallback_value


def times_from_header(header):
    """
        Return modified Julian date from header
    """
    date    = header['DATE-OBS'].strip()
    mtime   = header['MP-MID'].strip()
    time_start = Time(date, format='isot', scale='utc')
    if mtime == '00:00:00.000000':
        # Most probably the time point was not taken
        try:
            exp_time = header['EXPOSURE']
            if 'millisec' in header.comments['EXPOSURE'].lower():
                exp_time = exp_time / 1000
        except:
            exp_time = float(header['EXPTIME'])
            if 'millisec' in header.comments['EXPTIME'].lower():
                exp_time = exp_time / 1000
        time_weighted = time_start + TimeDelta(exp_time/2., format='sec')
    else:
        time_weighted = Time(date[:10]+'T'+mtime, format='isot', scale='utc')
        if time_start > time_weighted:
            time_weighted += TimeDelta(1., format='jd')
    return date, time_start, time_weighted


def rotate_and_save(NightLogObject, n_rotations=1):
    """
        Rotate the images of all files in the NightLog-Object
        and save them in the same directory (but under new names).
    """
    allfiles = NightLogObject['all']
    for filename in allfiles:
        flux, header = load_raw_file(filename)
        spectrum = components.ScienceSpec(flux)
        spectrum.rotate(n_rotations=n_rotations)
        new_filename = splitext(filename)[0] + '_rot' + \
                        splitext(filename)[1]
        spectrum.save_spectrum(filename=new_filename, header=header,
                               fmt=splitext(filename)[1][1:])


def MedianCombine(imagelist, zero_level='none'):
    """
        Median combine a list of images
    """
    n = len(imagelist)
    if n==0:
        raise ValueError('Empty list provided!')
    
    im1, header1 = load_raw_file(imagelist[0])
    """ Where do we get the RON from?! """
    #ron1,gain1 = float(h1['CCDRON']),float(h1['BSCALE'])
    ron1, gain1 = 8.0, 4.0
    im1 = OverscanTrim(im1).astype('float')
    
    if zero_level != 'none':
        zero_data, header = load_raw_file(zero_level)
        im1 -= zero_data.astype('float')
    
    factor = 1.25
    if n < 3:
        factor = 1
    ron1 = factor * ron1 / np.sqrt(n)
    
    if n > 1:
        for i in range(n-1):
            im, header = load_raw_file(imagelist[i+1])
            if zero_level == 'none':
                im1 = np.dstack((im1, OverscanTrim(im)))
            else:
                im1 = np.dstack((im1, OverscanTrim(im) - zero_data))
        im1 = np.median(im1, axis=2)
    Spec = CalibrationSpec(flux=im1, header=header1)
    Spec.readout_noise, Spec.gain = ron1, gain1
    return Spec, ron1, gain1


def OverscanTrim(d):
    """
        Overscan correct and Trim a refurbished CORALIE image
    """
    """
    # bias has no significant structure, so a single median suffices, I think
    # overscan = [0:49] [2097:2145]
    data = d[:,51:2099]
    ov1  = d[:,:50]
    ov2  = d[:,2098:]
    ov = np.vstack((ov1.transpose(),ov2.transpose()))
    overscan = median(ov,axis=0)
    overscan = np.array([overscan]*2048)
    overscan = overscan.transpose()
    newdata = data - overscan
    """
    newdata = d.copy()
    return newdata
