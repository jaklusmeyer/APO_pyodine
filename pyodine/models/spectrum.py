import numpy as np
from scipy.interpolate import splrep, splev

from ..lib.misc import rebin, osample
from .base import DynamicModel, ParameterSet
import matplotlib.pyplot as plt


class SimpleModel(DynamicModel):
    """A working implementation of a :class:'DynamicModel'
    
    This child class incorporates methods to guess fitting parameters for each
    chunk and build a model spectrum.
    """

    param_names = ['velocity', 'tem_depth', 'iod_depth']
    
    # I moved the osample_factor to the base class DynamicModel
    #options = {
    #    'osample_factor': 4,  # Default oversampling factor # used to be 4.0
    #}

    def guess_params(self, chunk):
        """Guess all model parameters
        
        Basically this just calls the respective methods of the underlying
        submodels.
        
        Args:
            chunk (:class:'Chunk'): The chunk for which to guess the parameters.
        
        Return:
            :class:'ParameterSet': The parameter guesses.
        """
        
        params = ParameterSet(velocity=0., tem_depth=1., iod_depth=1.)  # TODO: Improve these guesses!
        params.add(self.lsf_model.guess_params(chunk), prefix='lsf')
        params.add(self.wave_model.guess_params(chunk), prefix='wave')
        params.add(self.cont_model.guess_params(chunk), prefix='cont')
        return params

    def eval(self, chunk, params, require=None, chunk_ind=None):#, lsf_fixed=None):
        """Evaluate model for one :class:'Chunk' using a given 
        :class:'ParameterSet'
        
        Args:
            chunk (:class:'Chunk'): Evaluate the model for this chunk.
            params (:class:'ParameterSet'): The model parameters to use.
            require (Optional[str]): If require='full', the iodine atlas and 
                Doppler-shifted template is required to cover the full range 
                of the chunk. Default is None.
            chunk_ind (Optional[int]): The index of the chunk within the
                observation.
        
        Return:
            ndarray[nr_pix]: The model spectrum for this chunk.
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

        # TODO: Iodine corresponds to osample=10.0. Use that as the fine grid?

        #
        # IODINE:
        #

        # Load iodine atlas
        # TODO: Optimize by loading this only once?
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
                print('NaN value detected in iodine function.')
                print(iod_fine)

        #
        # STELLAR TEMPLATE:
        #

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
                print('NaN value detected in template function.')
                """
                print(np.where(np.isnan(flux_tem)))
                print(len(tck[0]), len(tck[1]))
                print(len(tem.wave), len(flux_tem))
                
                plt.plot(wave_fine, iod_fine, alpha=0.5)
                plt.plot(chunk.wave, chunk.flux/np.max(chunk.flux), alpha=0.5)
                plt.plot(tem.wave*doppler, flux_tem, alpha=0.5)
                plt.show()"""
                

        #
        #  CONVOLUTION:
        #

        # Convolve with the LSF
        # If fixed LSF given, use that one:
        if self.lsf_array is not None:
            lsf = self.lsf_model.eval(self.lsf_array, lsf_params)
        else:
            x_lsf = self.lsf_model.generate_x(self.osample_factor, self.conv_width)
            lsf = self.lsf_model.eval(x_lsf, lsf_params)
        spec_fine = np.convolve(iod_fine * tem_fine, lsf, 'same')
        if any(np.isnan(spec_fine)):
            print('NaN value detected in model function (1).')
            print(spec_fine)
        # Resample back to original grid, given by wavelength
        spec_obs = rebin(wave_fine, spec_fine, wave_obs)
        if any(np.isnan(spec_obs)):
            print('NaN value detected in model function (2).')
            print(spec_obs)
        # Apply continuum
        spec_obs *= self.cont_model.eval(chunk.pix, cont_params)
        if any(np.isnan(spec_obs)):
            print('NaN value detected in model function (3).')
            print(spec_obs)

        return spec_obs
