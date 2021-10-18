import numpy as np
from .base import ParameterSet
from .shapes import LinearStaticModel


class ContinuumModel:
    pass


class LinearContinuumModel(ContinuumModel, LinearStaticModel):
    """A linear continuum model"""
    @staticmethod
    def guess_params(chunk):
        """Make an educated guess of the continuum parameters for a given chunk
        
        Args:
            chunk (:class:'Chunk'): The chunk for which to guess the parameters.
        
        Returns:
            :class:'ParameterSet': The guessed parameters (continuum zero
                point and slope).
        """
        if chunk.cont is not None:
            # Fit a straight line to the continuum and scale to fit the flux
            p = np.polyfit(chunk.pix, chunk.cont, 1)
            intercept = np.median(chunk.flux)
            slope = p[0] / p[1] * intercept
            return ParameterSet(intercept=intercept, slope=slope)
        else:
            # Fit a straight line to the spectral flux
            p = np.polyfit(chunk.pix, chunk.flux, 1)
            return ParameterSet(intercept=p[1], slope=p[0])
