import numpy as np
from .base import ParameterSet
from .shapes import LinearStaticModel


class WavelengthModel:
    pass


class LinearWaveModel(WavelengthModel, LinearStaticModel):
    """A linear wavelength model"""
    @staticmethod
    def guess_params(chunk):
        """Make an educated guess of the wavelength parameters for a given chunk
        
        Args:
            chunk (:class:'Chunk'): The chunk for which to guess the parameters.
        
        Returns:
            :class:'ParameterSet': The guessed parameters (wavelength zero
                point and dispersion).
        """
        p = np.polyfit(chunk.pix, chunk.wave, 1)
        return ParameterSet(intercept=p[1], slope=p[0])