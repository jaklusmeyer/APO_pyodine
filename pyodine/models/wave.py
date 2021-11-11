import numpy as np
from .base import ParameterSet
from .shapes import LinearStaticModel, ParabolicStaticModel


class LinearWaveModel(LinearStaticModel):
    """A linear wavelength model"""
    @staticmethod
    def guess_params(chunk):
        """Make an educated guess of the wavelength parameters for a given chunk
        
        :param chunk: The chunk for which to guess the parameters.
        :type chunk: :class:`Chunk`
        
        :return: The guessed parameters (wavelength zero point and slope).
        :rtype: :class:`ParameterSet`
        """
        p = np.polyfit(chunk.pix, chunk.wave, 1)
        return ParameterSet(intercept=p[1], slope=p[0])


class ParabolicWaveModel(ParabolicStaticModel):
    """A 2nd degree polynomial wavelength model"""
    @staticmethod
    def guess_params(chunk):
        """Make an educated guess of the wavelength parameters for a given chunk
        
        :param chunk: The chunk for which to guess the parameters.
        :type chunk: :class:`Chunk`
        
        :return: The guessed parameters (wavelength zero point and slope).
        :rtype: :class:`ParameterSet`
        """
        p = np.polyfit(chunk.pix, chunk.wave, 2)
        return ParameterSet(p0=p[0], p1=p[1], p2=p[2])