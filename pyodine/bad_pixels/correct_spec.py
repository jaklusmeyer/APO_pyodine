#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:12:20 2020

@author: paul
"""

import numpy as np
from ..components import Observation, NormalizedObservation

def correct_spectrum(spec, weights, orders):
    """Correct bad pixel regions in spectra
    Used to correct regions of spectra with weights=0 due to bad pixels etc.
    These regions are simply linearly interpolated from their boundaries.
    
    Args:
        spec (:class:'Observation'): The spectrum to correct.
        weights (ndarray[nr_ord,nr_pix]): The mask which marks the bad pixels
            (being 0 there).
        orders (ndarray[nr_ord_correct]): This array indicates which orders
            to correct.
    Return:
        :class:'Observation': The corrected spectrum.
    """
    for i, o in enumerate(orders):
        spec_order = spec[o]
        ind = np.where(weights[o]==0.)[0]
        
        if len(ind) > 0:
            start_pix = [ind[0]-1]
            end_pix = []
            for j in range(len(ind)-1):
                if ind[j+1] - ind[j] > 1:
                    start_pix.append(ind[j+1]-1)
                    end_pix.append(ind[j]+1)
            end_pix.append(ind[-1]+1)
            
            for j in range(len(start_pix)):
                spec_order.flux[start_pix[j]:end_pix[j]] = np.linspace(
                    spec_order.flux[start_pix[j]], spec_order.flux[end_pix[j]], end_pix[j]-start_pix[j])
            
            if isinstance(spec, NormalizedObservation):
                spec[o] = spec_order.flux
            elif isinstance(spec, Observation):
                spec[o].flux = spec_order.flux
        else:
            print('No zero weights: Order {}!'.format(o))
    
    return spec