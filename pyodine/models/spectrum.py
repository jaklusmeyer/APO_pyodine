import numpy as np
from scipy.interpolate import splrep, splev
import logging
import sys

from ..lib.misc import rebin, osample
from .base import DynamicModel, ParameterSet


class SimpleModel(DynamicModel):
    """A working implementation of a :class:'DynamicModel'
    
    This child class incorporates methods to guess fitting parameters for each
    chunk and build a model spectrum.
    """
    
    param_names = ['velocity', 'tem_depth', 'iod_depth']
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')

    def guess_params(self, chunk):
        """Guess all model parameters
        
        Basically this just calls the respective methods of the underlying
        submodels.
        
        :param chunk: The chunk for which to guess the parameters.
        :type chunk: :class:`Chunk`
        
        :return: The parameter guesses.
        :rtype: :class:`ParameterSet`
        """
        
        params = ParameterSet(velocity=0., tem_depth=1., iod_depth=1.)  # TODO: Improve these guesses!
        params.add(self.lsf_model.guess_params(chunk), prefix='lsf')
        params.add(self.wave_model.guess_params(chunk), prefix='wave')
        params.add(self.cont_model.guess_params(chunk), prefix='cont')
        return params

    def eval(self, chunk, params, require=None, chunk_ind=None):#, lsf_fixed=None):
        """Evaluate model for one :class:`Chunk` using a given 
        :class:`ParameterSet`
        
        :param chunk: Evaluate the model for this chunk.
        :type chunk: :class:`Chunk`
        :param params: The model parameters to use.
        :type params: :class:`ParameterSet`
        :param require: If require='full', the iodine atlas and Doppler-shifted 
            template is required to cover the full range of the chunk. Default 
            is None.
        :type require: str, or None
        :param chunk_ind: The index of the chunk within the observation.
        :type chunk_ind: int, or None
        
        :return: The model spectrum for this chunk.
        :rtype: ndarray[nr_pix]
        """

        lsf_params = params.filter(prefix='lsf')
        wave_params = params.filter(prefix='wave')
        cont_params = params.filter(prefix='cont')

        # Generate the "observed" wavelength grid, to be returned in the end
        wave_obs = self.wave_model.eval(chunk.pix, wave_params)

        # Generate the "fine" wavelength grid, used for convolution.
        # Extends beyond the chunk limits as defined by chunk.padding.
        pix_fine = osample(chunk.padded.pix, self.osample_factor)
        wave_fine = self.wave_model.eval(pix_fine, wave_params)

        #######################################################################
        ## IODINE:
        #######################################################################

        # Load iodine atlas
        iod = self.iodine_atlas.get_wavelength_range(
            wave_fine[0],
            wave_fine[-1],
            require=require
        )
        # Ensure "normalization" to mean value 1.0
        # FIXME: Normalization might depend on selected wavelength range?
        flux_iod = iod.flux / np.mean(iod.flux)

        # Scale depth of iodine atlas
        flux_iod = params['iod_depth'] * (flux_iod - 1.0) + 1.0

        # Interpolate iodine atlas to the fine grid
        # (Extrapolation may happen, but if keyword `require` is set to 'full',
        #  the iodine atlas will only load if it covers the full wavelength range)
        tck = splrep(iod.wave, flux_iod, s=0)
        iod_fine = splev(wave_fine, tck, der=0)
        if any(np.isnan(iod_fine)):
                logging.error('NaN value detected in iodine function.')

        #######################################################################
        ## STELLAR TEMPLATE:
        #######################################################################

        # Load stellar template
        if self.stellar_template is None:
            tem_fine = np.ones(len(pix_fine))  # For O-star fitting
        else:
            # Calculate relativistic doppler factor
            beta = params['velocity'] / 299792458.
            doppler = np.sqrt((1. + beta) / (1. - beta))
            
            if chunk_ind is None:
                # Fetch the relevant part of the template
                # (padded wavelength range with velocity shift)
                tem = self.stellar_template.get_wavelength_range(
                    wave_fine[0] / doppler,
                    wave_fine[-1] / doppler,
                    require=require
                )
            else:
                tem = self.stellar_template[chunk_ind]

            # Ensure "normalization" to mean value 1.0
            # FIXME: Do something more sophisticated here
            flux_tem = tem.flux / np.mean(tem.flux)

            # Scale depth of stellar template
            flux_tem = params['tem_depth'] * (flux_tem - 1.0) + 1.0

            # Interpolate shifted stellar template to the fine grid
            # (Extrapolation may happen, but if keyword `require` is set to
            #  'full', the iodine atlas will only load if it covers the full
            #  wavelength range)
            tck = splrep(tem.wave * doppler, flux_tem, s=0)
            tem_fine = splev(wave_fine, tck, der=0)
            if any(np.isnan(tem_fine)):
                logging.error('NaN value detected in template function.')
        
        #######################################################################
        ## LSF
        #######################################################################
        
        # If fixed LSF given, use that one:
        """
        if self.lsf_array is not None:
            lsf = self.lsf_model.eval(self.lsf_array, lsf_params)
        else:
            x_lsf = self.lsf_model.generate_x(
                    self.osample_factor, self.conv_width)
            lsf = self.lsf_model.eval(x_lsf, lsf_params)
        """
        
        x_lsf, lsf = self.eval_lsf(params)
        
        #######################################################################
        ## CONVOLUTION:
        #######################################################################

        # Convolve the model with the LSF
        spec_fine = np.convolve(iod_fine * tem_fine, lsf, 'same')
        
        if any(np.isnan(spec_fine)):
            logging.error('NaN value detected in oversampled model function.')
        
        # Resample back to original grid, given by wavelength
        spec_obs = rebin(wave_fine, spec_fine, wave_obs)
        
        # Apply continuum
        spec_obs *= self.cont_model.eval(chunk.pix, cont_params)

        return spec_obs
