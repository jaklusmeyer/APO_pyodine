from os.path import splitext, abspath#, isdir

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits as pyfits
import h5py
from astropy.time import Time, TimeDelta
import glob

from CERES import components

from utilities_waltz import conf


class WaltzSpectrum(components.ScienceSpec):
    
    def __init__(self, filename=None, instrument=None, obj=None, flux=None, 
                 header=None, rotate=False, flip=False):
        # With given filename, load an existing spectrum
        if filename is not None:
            flux, header = load_raw_file(filename)
            super().__init__(flux)
            self.orig_filename = abspath(filename)
        elif flux is not None:
            super().__init__(flux)
            self.orig_filename = None
        else:
            raise AttributeError('No filename and flux are given.')
        
        if rotate is True:
            self.rotate()
        if flip is True:
            self.flip()
        
        self.orig_header = header
        if self.orig_header is None:
            self.orig_header = {}
            self.star = obj
            self.instrument = instrument
            self.iodine_in_spectrum = None
            
            # Camera details
            self.exp_time = None
            self.flux_level = None      # FIXME: Define a flux measure
            self.gain, self.readout_noise = None, None
            
            # Timing
            self.date, self.time_start, self.time_weighted = None, None, None
        else:            
            self.instrument = get_instrument(header)
            self.star = get_star(header)
            #self.iodine_in_spectrum = None # for Lick science-specs always true
            self.iodine_in_spectrum = check_iodine_cell(header)
            
            # Camera details
            self.exp_time = or_none(header, 'EXPTIME')
            self.flux_level = None      # FIXME: Define a flux measure
            self.gain, self.readout_noise = get_gain_ron(header)
            
            # Timing
            self.date, self.time_start, self.time_weighted = times_from_header(header)
            self.date_bjd = or_none(header, 'DATE-BJD')
        
    def save_spectrum(self, filename=None, header=None, fmt=None):
        """
           Save spectrum in fits or hdf5 format. Extending same method in
           components.ScienceSpec with the header.
        """
        if not header:
            header = self.orig_header
        super(WaltzSpectrum, self).save_spectrum(filename, header, fmt)
        

"""
class WaltzScienceSpec(components.ScienceSpec):

    # Custom properties
    _spec = None    # Internal storage of spectral flux

    def __init__(self, filename, instrument=None, star=None, rotate=False, flip=False):
        flux, header = load_raw_file(filename)

        super().__init__(flux)
        # If rotate is true, rotate the spectrum by 90deg
        if rotate:
            self.rotate()
        # If flip is true, flip the spectrum in axis 0 
        # (should be cross-disp. direction)
        if flip:
            self.flip()
        
        self.header = header
        self.orig_filename = abspath(filename)

        self.instrument = instrument or get_instrument(header)
        self.star = star or get_star(header)
        self.iodine_in_spectrum = True # for Lick science-specs always true
        #self.iodine_in_spectrum, self.iodine_cell_id = check_iodine_cell(header)

        # Camera details
        self.exp_time = or_none(header, 'EXPTIME')
        self.flux_level = None      # FIXME: Define a flux measure
        self.gain, self.readout_noise = get_gain_ron(header)
        self.dark_current = None    # FIXME: Not in header

        # Timing
        self.date, self.time_start, self.time_weighted = times_from_header(header)
        
        # TODO: bary_date_weighted and deal with times above
        # TODO: Implement flux check
        # TODO: Re-calculate BVC
    
    def save_spectrum(self, filename=None, header=None, fmt=None):
        """"""
           Save spectrum in fits or hdf5 format. Extending same method in
           components.ScienceSpec with the header.
        """"""
        if not header:
            header = self.header
        super(WaltzScienceSpec, self).save_spectrum(filename, header, fmt)


class CalibrationSpec(components.ScienceSpec):
    
    def __init__(self, filename=None, instrument=None, flux=None, header=None, rotate=False, flip=False):
        if type(flux) is not np.ndarray:
            flux, header = load_raw_file(filename)
        super().__init__(flux)
        self.header = header
        self.orig_filename = None if filename is None else abspath(filename)
        
        if instrument:
            self.instrument = instrument 
        elif header:
            self.instrument = get_instrument(header)
        else:
            self.instrument = None
        
        # Camera details
        self.flux_level = None      # FIXME: Define a flux measure
        if header:
            self.gain, self.readout_noise = get_gain_ron(header)
        else:
            self.gain, self.readout_noise = None, None
        self.dark_current = None    # FIXME: Not in header
        
    def save_spectrum(self, filename=None, header=None, fmt=None):
        """"""
           Save spectrum in fits or hdf5 format. Extending same method in
           components.ScienceSpec with the header.
        """"""
        if not header:
            header = self.header
        super(CalibrationSpec, self).save_spectrum(filename, header, fmt)
"""

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
                'science': [],
                }
        self.infos = {}
        for key in self.images.keys():
            self.infos[key+'_mjd'] = []
        self.infos['science_obj'] = []
        #self.biases = {}
        #self.darks  = {}
        #self.flats  = {}
        #self.thars  = {}
        #self.i2s    = {}
        #self.science = {}
        #self.imagetypes = ['biases', 'darks', 'flats', 'thars', \
        #                   'i2s', 'science']
    
    def classify_files(self, constraint='*fits'):
        """
            Classifies all files in a directory
            and writes a night log of science images
        """
        all_files = sorted(glob.glob(self.directoryname + '/' + constraint))
        for filename in all_files:
            if filename not in self['all']:
                with pyfits.open(filename) as h:
                    header = h[0].header
                    date, time_start, time_weighted = times_from_header(header)
                    # Biases of the Andor CCD have exposure times of 1e-5 in header
                    if header['SHUTTER'] == 2 and \
                    header['EXPTIME'] <= 1.2E-5:
                        self.images['biases'].append(filename)
                        self.infos['biases_mjd'].append(time_weighted.mjd)
                    elif header['SHUTTER'] == 2 and \
                    header['EXPTIME'] != 0.0:
                        self.images['darks'].append(filename)
                        self.infos['darks_mjd'].append(time_weighted.mjd)
                    elif header['SHUTTER'] == 0 and \
                    header['PCUPSTAT'].lower().strip() == 'out' and \
                    header['LEDPOWER'].lower().strip() == 'on' and \
                    header['I2CESTAT'].lower().strip() == 'out':
                        self.images['flats'].append(filename)
                        self.infos['flats_mjd'].append(time_weighted.mjd)
                    elif header['SHUTTER'] == 0 and \
                    header['PCUPSTAT'].lower().strip() == 'out' and \
                    header['LEDPOWER'].lower().strip() == 'on' and \
                    header['I2CESTAT'].lower().strip() == 'in':
                        self.images['i2s'].append(filename)
                        self.infos['i2s_mjd'].append(time_weighted.mjd)
                    elif header['SHUTTER'] == 0 and \
                    header['PCUPSTAT'].lower().strip() == 'out' and \
                    header['THAPOWER'].lower().strip() == 'on':
                        self.images['thars'].append(filename)
                        self.infos['thars_mjd'].append(time_weighted.mjd)
                    elif header['SHUTTER'] == 0 and \
                    header['PCUPSTAT'].lower().strip() == 'in':
                        self.images['science'].append(filename)
                        self.infos['science_mjd'].append(time_weighted.mjd)
                        self.infos['science_obj'].append(get_star(header).name)
                        #self.sciens[filename]['object'] = header['OBJECT'].upper()
                        #self.sciens[filename]['exp_time'] = float(header['EXPTIME'])
                    """
                    elif header['SHUTTER'] == 0 and \
                    header['PCUPSTAT'].lower().strip() == 'out' and \
                    header['CLMRSTAT'].lower().strip() == 'in' and \
                    header['I2CESTAT'].lower().strip() == 'out':
                        self.images['flats'].append(filename)
                        self.infos['flats_mjd'].append(time_weighted.mjd)
                    elif header['SHUTTER'] == 0 and \
                    header['PCUPSTAT'].lower().strip() == 'out' and \
                    header['CLMRSTAT'].lower().strip() == 'in' and \
                    header['I2CESTAT'].lower().strip() == 'in':
                        self.images['i2s'].append(filename)
                        self.infos['i2s_mjd'].append(time_weighted.mjd)
                    elif header['SHUTTER'] == 0 and \
                    header['PCUPSTAT'].lower().strip() == 'out' and \
                    header['CLMRSTAT'].lower().strip() == 'out':
                        self.images['thars'].append(filename)
                        self.infos['thars_mjd'].append(time_weighted.mjd)
                    """
                    
    
    def __getitem__(self, *args):
        typename, = args
        if isinstance(typename, str):
            if typename in self.images.keys():
                return self.images[typename]
            elif typename in self.infos.keys():
                return self.infos[typename]
            elif typename == 'all':
                allfiles = []
                for key in self.images.keys():
                    allfiles += self.images[key]
                return allfiles
            else:
                raise ValueError(
                        'Image- or infotype not existing. Choose one of {}, {}, or "all".'.format(
                        [key for key in self.images.keys()], [key for key in self.infos.keys()])
                        )
        else:
            raise ValueError(
                    'Please enter a string. Choose one of {}, {}, or "all".'.format(
                    [key for key in self.images.keys()], [key for key in self.infos.keys()])
                    )        
    
    def __len__(self):
        length = 0
        for key in self.images.keys():
            length += len(self.images[key])
        return length
    
    def __str__(self):
        outstring = '<NightLog (science:'
        for spec in self.images['science']:
            outstring += '\n\t{},'.format(spec)
        outstring += '\n  '
        for key in self.images.keys():
            if key != 'science':
                outstring += '{}  '.format(key)
        outstring += '\n  '
        
        outstring += '  {}       {}      {}      {}     {}   )>'.format(
                len(self.images['biases']), len(self.images['darks']),
                len(self.images['flats']), len(self.images['thars']),
                len(self.images['i2s']))
        return outstring
    
    def sort_by_time(self, typename='all'):
        if isinstance(typename, str):
            if typename in self.images.keys():
                if len(self.infos[typename+'_mjd']) > 1:
                    sorted_ind = np.argsort(self.infos[typename+'_mjd'])
                    self.images[typename] = [self[typename][arg] for arg in sorted_ind]
                    self.infos[typename+'_mjd'] = [self.infos[typename+'_mjd'][arg] for arg in sorted_ind]
            elif typename == 'all':
                for key in self.images.keys():
                    self.sort_by_time(key)
            else:
                raise ValueError(
                        'Imagetype not existing. Choose one of {}'.format(
                        [key for key in self.images.keys()])
                        )
    
    def save(self, filename='Log.txt'):
        outstring = ''
        with open(filename, 'w') as f:
            for key in self.images.keys():
                if len(self.images[key]) > 0:
                    outstring = key + ':\n'
                    for i, specname in enumerate(self.images[key]):
                        outstring += '\t{}\t{}\n'.format(specname, self.infos[key+'_mjd'][i])
                    f.write(outstring)
    
    def remove_by_name(self, name):
        for key in self.images.keys():
            if name in self[key]:
                ind = self[key].index(name)
                del self[key][ind]
                for infokey in self.info.keys():
                    if key in infokey:
                        del self.info[infokey][ind]
                return
        print('Element {} was not found in this LobObject!'.format(name))


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
    try:
        name = header['OBJECT']
        name = name.strip()
    except Exception as e:
        print(e)
        name = None
    try:
        coordinates = SkyCoord(
            header['RA'] + ' ' + header['DEC'],
            unit=(u.hourangle, u.deg)
        )
    except Exception as e:
        # TODO: Log this event
        print(e)
        coordinates = None
    # Get the proper motion vector
    try:
        proper_motion = (header['RA_PM'], header['DEC_PM'])
    except Exception as e:
        # TODO: Log this event
        print(e)
        proper_motion = (None, None)
    return components.Star(name, coordinates=coordinates, proper_motion=proper_motion)


def get_instrument(header) -> components.Instrument:
    """
        Determine the instrument from the header and return Instrument object
    """
    try:
        if 'INSTRUME' in header and \
        'Hamilton Spec.' in header['INSTRUME']:
            return conf.my_instruments['lick']
        elif 'TELESCOP' in header and \
        'Waltz Telescope' in header['TELESCOP']:
            return conf.my_instruments['waltz']
    except Exception as e:
        # TODO: Log this event
        #raise TypeError('Could not determine instrument')
        print(e)
        return None


def or_none(header, key, fallback_value=None):
    try:
        return header[key]
    except KeyError as e:
        print(e)
        # TODO: Log this event
        return fallback_value


def times_from_header(header):
    """
        Return modified Julian date from header
    """
    try:
        date    = header['DATE-OBS'][:10]
        stime   = header['DATE-OBS']
        mtime   = header['OBS-MID']
        time_start = Time(stime, format='isot', scale='utc')
        if '---' in mtime:
            exp_time = header['EXPTIME']
            time_weighted = time_start + TimeDelta(exp_time/2., format='sec')
        else:
            time_weighted = Time(mtime, format='isot', scale='utc')
    except Exception as e:
        print(e)
        date, time_start, time_weighted = None, None, None
    return date, time_start, time_weighted


def get_gain_ron(header):
    """
       Return gain and readout noise from header
    """
    try:
        gain_name = header['PREGNAME'].strip()
        # Readout noise depends on the HSS used
        # (from the Andor iKon L936 manual)
        #readout_noises = {
        #        5.0: 31.5,
        #        3.0: 11.7,
        #        1.0: 7.0,
        #        0.05: 2.9
        #        }
        hss = float(header['HSSPEED'])
        # Both gain and readout noise depend
        # on the gain and HSS settings
        gain = CCD_gain[round(hss, 2)][gain_name]
        readout_noise = CCD_ron[round(hss, 2)][gain_name]
        #readout_noise = readout_noises[round(hss, 2)]
        print(gain, readout_noise)
        return gain, readout_noise
    except Exception as e:
        print(e)
        return None, None


def check_iodine_cell(header):
    """
       Check whether iodine cell was in or out of spectrum.
    """
    try:
        i2state = header['I2CESTAT'].strip()
        if 'in' in i2state:
            return True
        else:
            return False
    except Exception as e:
        print(e)
        return None


def rotate_and_save(NightLogObject, n_rotations=1, name_addon=''):
    """
        Rotate the images of all files in the NightLog-Object
        and save them in the same directory (maybe under new names).
    """
    allfiles = NightLogObject['all']
    for filename in allfiles:
        flux, header = load_raw_file(filename)
        spectrum = components.ScienceSpec(flux)
        spectrum.rotate(n_rotations=n_rotations)
        new_filename = splitext(filename)[0] + name_addon + \
                        splitext(filename)[1]
        spectrum.save_spectrum(filename=new_filename, header=header,
                               fmt=splitext(filename)[1][1:])


def flip_and_save(NightLogObject, name_addon=''):
    """
        Flip the images of all files in the NightLog-Object
        (in cross-disp. direction)
        and save them in the same directory (but under new names).
    """
    allfiles = NightLogObject['all']
    for filename in allfiles:
        flux, header = load_raw_file(filename)
        spectrum = components.ScienceSpec(flux)
        spectrum.flip()
        new_filename = splitext(filename)[0] + name_addon + \
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
    gain1, ron1 = get_gain_ron(header1)
    #ron1, gain1 = 8.0, 4.0
    im1 = OverscanTrim(im1).astype('float')
    
    if zero_level != 'none':
        zero_data, header = load_raw_file(zero_level)
        im1 -= zero_data.astype('float')
    
    # Does this make sense?
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
    Spec = WaltzSpectrum(flux=im1, header=header1)
    Spec.readout_noise, Spec.gain = ron1, gain1
    return Spec


CCD_ron = {
        5.0: {'1x': 61.0,
              '2x': 43.0,
              '4x': 32.6
                },
        3.0: {'1x': 36.7,
              '2x': 14.1,
              '4x': 11.4
                },
        1.0: {'1x': 10.0,
              '2x': 7.7,
              '4x': 6.5
                },
        0.05: {'1x': 4.0,
               '2x': 3.2,
               '4x': 2.9
                }
        }

CCD_gain = {
        5.0: {'1x': 9.5,
              '2x': 4.9,
              '4x': 2.6
                },
        3.0: {'1x': 4.7,
              '2x': 2.4,
              '4x': 1.3
                },
        1.0: {'1x': 4.1,
              '2x': 2.2,
              '4x': 1.2
                },
        0.05: {'1x': 4.0,
               '2x': 2.2,
               '4x': 1.2
                }
        }


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
