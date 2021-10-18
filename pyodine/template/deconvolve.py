import sys
import numpy as np
from ..components import Spectrum, TemplateChunk
from ..template.base import StellarTemplate, StellarTemplate_Chunked
from ..lib import misc

from progressbar import ProgressBar
#from ..models.lsf import LSFModel
#from scipy import signal


#class Deconvolver:
#    """Abstract base class"""
#    def deconvolve_obs(self, normalized_observation, velocity_offset):
#        """Deconvolve an Observation (multiorder)"""
#        raise NotImplementedError
#
#    def deconvolve_single(self, normalized_spectrum, order):
#        """Deconvolve a single spectrum"""
#        raise NotImplementedError


class SimpleDeconvolver():#Deconvolver):
    """Deconvolver to create a :class:'StellarTemplate' with stitched chunks
    
    A simple deconvolver, that deconvolves all chunks of a stellar spectrum
    with the help of a LSF (either directly from the O-star model, or handed in 
    manually with lsf_fixed), stitches the deconvolved chunks together and thus
    produces a 2D template (orders and pixels).
    
    Args:
        ostar_chunks (:class:'ChunkArray'): Chunks of the modelled O-star
            observation.
        ostar_model (:class:'SimpleModel'): The employed model for the fitting
            of the O-star observation.
        ostar_params (list): The :class:'ParameterSet' objects for each O-star
            chunk containing the best-fit results.
    """
    def __init__(self, ostar_chunks, ostar_model, ostar_params):
        self.ostar_chunks = ostar_chunks
        self.ostar_model = ostar_model
        self.ostar_params = ostar_params

    def deconvolve_obs(self, normalized_observation, velocity_offset, bary_v,
                       lsf_fixed=None, deconv_pars=None):
        """Deconvolve all orders in 'normalized_observation' and return a
        :class:'StellarTemplate' object
        
        Args:
            normalized_observation (:class:'NormalizedObservation'): The 
                normalized stellar observation spectrum (without I2).
            velocity_offset (float): The velocity-offset of the stellar 
                template to the reference spectrum.
            bary_v (float): The barycentric velocity of the stellar template.
            lsf_fixed (Optional[ndarray[nr_chunks,nr_pix_lsf]]): If an array 
                with a pre-defined lsf (e.g. smoothed lsf) is given, this one
                is used in the deconvolution. Defaults to None.
            deconv_pars (Optional[dict]): A set of deconvolution parameters. If
                None is given, a hardcoded set is used (Default).
        
        Return:
            :class:'StellarTemplate': The deconvolved stellar template with
                chunks within an order stitched together.
        """
        
        if deconv_pars == None:
            deconv_pars = {'osample_temp': 10.0,
                           'jansson_niter': 1200,
                           'jansson_zerolevel': 0.00,
                           'jansson_contlevel': 1.02,
                           'jansson_conver': 0.1,
                           'jansson_chi_change': 1e-6,
                           'lsf_conv_width': 6.
                           }
        
        template = StellarTemplate(normalized_observation, velocity_offset=velocity_offset, 
                                   bary_vel_corr=bary_v, osample=deconv_pars['osample_temp'])
        
        for i in normalized_observation.orders:
            print('Deconvolve order ', i)
            sys.stdout.flush()  # TODO: Logging
            template[i] = self.deconvolve_single(normalized_observation[i], i, lsf_fixed, deconv_pars)
        return template

    def deconvolve_single(self, normalized_spectrum, order, lsf_fixed, deconv_pars):
        """Deconvolve a single order using LSFs from fitted O-star chunks and
        return a new :class:'Spectrum' object
        
        Args:
            normalized_spectrum (:class:'Spectrum'): An order of the
                normalized stellar observation.
            order (int): The order number to work on.
            lsf_fixed (ndarray[nr_chunks,nr_pix_lsf]): If an array 
                with a pre-defined lsf (e.g. smoothed lsf) is given, this one
                is used in the deconvolution. If None, the lsfs are constructed
                from the best-fit parameters from the O-star modelling.
            deconv_pars (dict): A set of deconvolution parameters.
        
        Return:
            :class:'Spectrum': The deconvolved order.
        """

        osample = deconv_pars['osample_temp'] #10.0    # Oversampling factor in Jansson deconvolution
        niter = deconv_pars['jansson_niter'] #1200      # This many iterations
        zerolevel = deconv_pars['jansson_zerolevel'] #0.00  # The zeropoint -- should always be zero, unless you know better
        contlevel = deconv_pars['jansson_contlevel'] #1.02  # Continuum level
        conver = deconv_pars['jansson_conver'] #0.2      # Convergence parameter in Jansson deconvolution #0.02

        # Pixel vector for the LSF
        nlsf = deconv_pars['lsf_conv_width']  # TODO: Warn if more than chunk padding
        nlsf_fine = int(nlsf * osample)
        pix_lsf = np.linspace(-int(nlsf), int(nlsf), 2 * nlsf_fine + 1)  # TODO: Use lsf_model.generate_x(osample_factor)?

        # Generate oversampled pixel vector before the loop in order
        # to ensure identical sampling in overlapping regions.
        npix = len(normalized_spectrum)
        full_pix_fine = np.linspace(0, npix, int(npix* osample + 1))

        # Results go here
        deconv_chunks = list()

        for i in self.ostar_chunks.get_order_indices(order):
            ochunk = self.ostar_chunks[i]

            # Oversampled pixel grid
            first = int(np.searchsorted(full_pix_fine, ochunk.padded.abspix[0], 'left'))
            last = int(np.searchsorted(full_pix_fine, ochunk.padded.abspix[-1], 'left'))
            indices = np.arange(first, last + 1, dtype='int')

            abspix_fine = full_pix_fine[indices]
            pix_fine = abspix_fine - ochunk.abspix[0] + ochunk.pix[0]

            # Get the corresponding pixels from the normalized spectrum
            normspec = normalized_spectrum[ochunk.padded.abspix]
            flux_fine = misc.rebin(ochunk.padded.abspix, normspec.flux, abspix_fine)

            # Get and evaluate the wavelength model
            wave_params = self.ostar_params[i].filter(prefix='wave')
            wave_fine = self.ostar_model.wave_model.eval(pix_fine, wave_params)

            # Get and evaluate the LSF model
            lsf_params = self.ostar_params[i].filter(prefix='lsf')
            if lsf_fixed is None:
                lsf = self.ostar_model.lsf_model.eval(pix_lsf, lsf_params)  # TODO: Maybe better to use the same x-vector as in the fit?
            else:
                lsf = lsf_fixed[i]
                #lsf = self.ostar_model.lsf_model.eval(lsf_fixed[i], lsf_params)
            
            flux_deconv = jansson(flux_fine, lsf, niter, a=zerolevel, b=contlevel, 
                                  delta=conver, chi_change=deconv_pars['jansson_chi_change'])  # b = max(flux)?
            jj = slice(nlsf_fine, -nlsf_fine)
            deconv_chunks.append({
                'flux': flux_deconv[jj],
                'wave': wave_fine[jj],
                'indices': indices[jj],
            })

        full_flux = -np.ones(len(full_pix_fine))
        full_wave = -np.ones(len(full_pix_fine))

        # Keep track of the first non-overlapping pixel
        alone_start_right = 0

        # Assume that the chunk list is ordered by wavelength
        for k in range(len(deconv_chunks) - 1):
            
            left = deconv_chunks[k]
            right = deconv_chunks[k + 1]

            # Offset between local and global indices
            offset_left = left['indices'][0]
            offset_right = right['indices'][0]

            # First index of the ranges without overlap
            alone_start_left = alone_start_right  # Previous dc2 value
            alone_start_right = int(np.searchsorted(right['indices'], left['indices'][-1], 'right'))

            # First index of the range with overlap in the left chunk
            overlap_start = len(left['indices']) - alone_start_right

            # Build the full flux vector and wavelength vector
            # Before overlap
            ii = slice(offset_left + alone_start_left, offset_left + overlap_start)
            full_flux[ii] = left['flux'][alone_start_left:overlap_start]
            full_wave[ii] = left['wave'][alone_start_left:overlap_start]

            # During overlap, weighted mean shifting linearly from left to right
            wleft = np.linspace(1, 0, alone_start_right)
            wright = np.linspace(0, 1, alone_start_right)
            ii = slice(offset_left + overlap_start, offset_right + alone_start_right)
            full_flux[ii] = wleft * left['flux'][overlap_start:] + wright * right['flux'][:alone_start_right]
            full_wave[ii] = wleft * left['wave'][overlap_start:] + wright * right['wave'][:alone_start_right]

        # Fill in data after the last overlap
        last = deconv_chunks[-1]
        ii = last['indices'][alone_start_right:]
        full_flux[ii] = last['flux'][alone_start_right:]
        full_wave[ii] = last['wave'][alone_start_right:]

        # Return the valid data points
        ii = np.where(full_flux != -1)
        return Spectrum(full_flux[ii], full_wave[ii])


class ChunkedDeconvolver():#Deconvolver):
    """Deconvolver to create a :class:'StellarTemplate_Chunked'
    
    A deconvolver, that deconvolves all chunks of a stellar spectrum
    with the help of a LSF (either directly from the O-star model, or handed in 
    manually with lsf_fixed), but does not stitch them together. Rather it 
    returns a :class:'StellarTemplate_Chunked' object.
    
    Args:
        ostar_chunks (:class:'ChunkArray'): Chunks of the modelled O-star
            observation.
        ostar_model (:class:'SimpleModel'): The employed model for the fitting
            of the O-star observation.
        ostar_params (list): The :class:'ParameterSet' objects for each O-star
            chunk containing the best-fit results.
    """
    def __init__(self, ostar_chunks, ostar_model, ostar_params):
        self.ostar_chunks = ostar_chunks
        self.ostar_model = ostar_model
        self.ostar_params = ostar_params

    def deconvolve_obs(self, normalized_observation, velocity_offset, bary_v, 
                       weights=None, lsf_fixed=None, deconv_pars=None):
        """Deconvolve all orders in 'normalized_observation' and return a
        :class:'StellarTemplate_Chunked' object
        
        Args:
            normalized_observation (:class:'NormalizedObservation'): The 
                normalized stellar observation spectrum (without I2).
            velocity_offset (float): The velocity-offset of the stellar 
                template to the reference spectrum.
            bary_v (float): The barycentric velocity of the stellar template.
            weights (Optional[ndarray[nr_chunks]]): If an array with chunk
                weights is supplied, these are passed into the
                :class:'StellarTemplate_Chunked' object. If None, weights are
                calculated analytically (default).
            lsf_fixed (Optional[ndarray[nr_chunks,nr_pix_lsf]]): If an array 
                with a pre-defined lsf (e.g. smoothed lsf) is given, this one
                is used in the deconvolution. Defaults to None.
            deconv_pars (Optional[dict]): A set of deconvolution parameters. If
                None is given, a hardcoded set is used (Default).
        
        Return:
            :class:'StellarTemplate_Chunked': The deconvolved stellar template.
        """
        
        if deconv_pars == None:
            deconv_pars = {'osample_temp': 10.0,
                           'jansson_niter': 1200,
                           'jansson_zerolevel': 0.00,
                           'jansson_contlevel': 1.02,
                           'jansson_conver': 0.2,
                           'jansson_chi_change': 1e-6,
                           'lsf_conv_width': 10.
                           }
        
        template = StellarTemplate_Chunked(normalized_observation, velocity_offset=velocity_offset, 
                                           bary_vel_corr=bary_v, osample=deconv_pars['osample_temp'])
        
        with ProgressBar(max_value=len(self.ostar_chunks), redirect_stdout=True) as bar:
            bar.update(0)
            print('Deconvolving chunks...')

            for i in range(len(self.ostar_chunks)):
                #print('Deconvolve chunk ', i)
                sys.stdout.flush()  # TODO: Logging
                temp_chunk = self.deconvolve_single_chunk(normalized_observation, i, deconv_pars, 
                                                          weights, lsf_fixed)
                template.append(temp_chunk)
                bar.update(i+1)
        
        return template

    def deconvolve_single_chunk(self, normalized_observation, i, deconv_pars, weights, lsf_fixed):
        """Deconvolve a single chunk using the lsf from the fitted O-star chunks 
        (or the fixed lsf if given) and return it as a :class:'TemplateChunk'
        object
        
        Args:
            normalized_observation (:class:'NormalizedObservation'): The 
                normalized stellar observation spectrum (without I2).
            i (int): The chunk number to work on.
            deconv_pars (dict): A set of deconvolution parameters.
            weights (ndarray[nr_chunks]): An array with chunk weights. If None, 
                weights are calculated analytically.
            lsf_fixed (ndarray[nr_chunks,nr_pix_lsf]): If an array 
                with a pre-defined lsf (e.g. smoothed lsf) is given, this one
                is used in the deconvolution. If None, the lsfs are constructed
                from the best-fit parameters from the O-star modelling.
        
        Return:
            :class:'TemplateChunk': The deconvolved chunk spectrum.            
        """

        osample = deconv_pars['osample_temp'] #10.0    # Oversampling factor in Jansson deconvolution
        niter = deconv_pars['jansson_niter'] #1200      # This many iterations
        zerolevel = deconv_pars['jansson_zerolevel'] #0.00  # The zeropoint -- should always be zero, unless you know better
        contlevel = deconv_pars['jansson_contlevel'] #1.02  # Continuum level
        conver = deconv_pars['jansson_conver'] #0.2      # Convergence parameter in Jansson deconvolution #0.02

        # Pixel vector for the LSF
        nlsf = deconv_pars['lsf_conv_width']  # TODO: Warn if more than chunk padding
        nlsf_fine = int(nlsf * osample)
        pix_lsf = np.linspace(-int(nlsf), int(nlsf), 2 * nlsf_fine + 1)  # TODO: Use lsf_model.generate_x(osample_factor)?
        
        ochunk = self.ostar_chunks[i]
        normalized_spectrum = normalized_observation[ochunk.order]
        
        # Generate oversampled pixel vector before the loop in order
        # to ensure identical sampling in overlapping regions.
        npix = len(normalized_spectrum)
        full_pix_fine = np.linspace(0, npix, int(npix* osample + 1))
        
        # Oversampled pixel grid
        first = int(np.searchsorted(full_pix_fine, ochunk.padded.abspix[0], 'left'))
        last = int(np.searchsorted(full_pix_fine, ochunk.padded.abspix[-1], 'left'))
        indices = np.arange(first, last + 1, dtype='int')

        abspix_fine = full_pix_fine[indices]
        pix_fine = abspix_fine - ochunk.abspix[0] + ochunk.pix[0]

        # Get the corresponding pixels from the normalized spectrum
        normspec = normalized_spectrum[ochunk.padded.abspix]
        flux_fine = misc.rebin(ochunk.padded.abspix, normspec.flux, abspix_fine)

        # Get and evaluate the wavelength model
        wave_params = self.ostar_params[i].filter(prefix='wave')
        wave_fine = self.ostar_model.wave_model.eval(pix_fine, wave_params)

        # Get and evaluate the LSF model
        lsf_params = self.ostar_params[i].filter(prefix='lsf')
        if lsf_fixed is None:
            lsf = self.ostar_model.lsf_model.eval(pix_lsf, lsf_params)  # TODO: Maybe better to use the same x-vector as in the fit?
        else:
            lsf = lsf_fixed[i]
        
        flux_deconv = jansson(flux_fine, lsf, niter, a=zerolevel, b=contlevel, 
                              delta=conver, chi_change=deconv_pars['jansson_chi_change'])  # b = max(flux)?
        
        w0 = wave_params['intercept']
        w1 = wave_params['slope']
        order = ochunk.order
        pix0 = ochunk.abspix[0]
        
        if weights is None:
            # Analytic chunk weights - need to reduce dispersion by template oversampling factor
            chunk_weight = misc.analytic_chunk_weights(flux_deconv, w0, w1/osample)
        else:
            chunk_weight = weights[i]
        
        # deconvolved flux, fine wavelength grid, fine pixel grid, wave intercept, dispersion (normal sampling),
        # order, pixel 0 of chunk within order, and chunk weight
        temp_chunk = TemplateChunk(flux_deconv, wave_fine, pix_fine, w0, w1, order, pix0, chunk_weight)
        return temp_chunk


def jansson(observed, lsf, niter, a=0.0, b=1.0, delta=0.1, chi_change=1e-8):
    """Jansson deconvolution algorithm
    
    P. Heeren: Added chi_change (from Lick dop code)
    
    Args:
        observed (ndarray[nr_pix]): The flux values of the observed spectrum.
        lsf (ndarray[nr_pix_lsf]): The Line Spread Function for the 
            deconvolution (assumed normalized to sum(lsf)=1.0).
        niter (int): Maximum number of iterations for the algorithm.
        a (Optional[float]): Zero-level of the deconvolved spectrum. Defaults
            to 0.0 (to enforce positivity).
        b (Optional[float]): Continuum level of the deconvolved spectrum.
            Defaults to 1.0 (to enforce the continuum).
        delta (Optional[float]): Convergence parameter. Defaults to 0.1.
        chi_change (Optional[float]): Minimum change of red. Chi**2 in Jansson 
            deconvolution before iterations are stopped.
    
    Return:
        ndarray[nr_pix]: The deconvolved spectrum.
    """
    
    old = observed
    kernel = lsf
    guess = np.convolve(old, kernel, 'same')

    old = guess + \
        delta * (1.0 - np.abs(guess - (a + b) / 2.) * 2. / (b - a)) * \
        (observed - np.convolve(guess, kernel, 'same'))   # Relaxation funtion
    
    chi = 1.
    k = -1
    while k <= niter:
        k += 1
        old_chi = chi
        # Relaxation function
        relax = delta * (1.0 - np.abs(old - (a + b) / 2.0) * 2.0 / (b - a))
        guess = np.convolve(old, kernel, 'same')
        convol1 = np.convolve(observed, kernel, 'same')
        convol2 = np.convolve(guess, kernel, 'same')
        # Update `old` with the new estimate
        old += relax * (convol1 - convol2)
        chi = np.std((observed-old)**2)
        if abs((old_chi/chi)-1.) < chi_change:
            #print(k, (old_chi/chi))
            k = niter + 1
    """
    # This is how it is done in Lick dop code
    spec = observed
    dsp = spec
    inpr = lsf
    fake_sp = np.convolve(dsp, inpr, 'same')
    diff = dsp - fake_sp
    
    chi = 1.
    k = -1
    while k <= niter:
        k += 1
        old_chi = chi
        # Smoothing
        diff = signal.wiener(diff, mysize=15, noise=None)
        # Relaxation function
        rcor = delta * (1.0 - np.abs(dsp - (a + b) / 2.0) * 2.0 / (b - a))
        dsp = dsp + rcor * diff
        fake_sp = np.convolve(dsp, inpr, 'same')
        old = dsp
        
        diff = spec - fake_sp
        
        chi = np.std(diff**2)
        if ((old_chi/chi)-1.) < chi_change:
            print(k)
            k = niter + 1
    """
    return old  # TODO: Return only the valid range?
