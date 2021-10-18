import numpy as np
import lmfit

from .base import Fitter, FitResult
from ..models.base import ParameterSet
from ..components import NoDataError


class LmfitWrapper(Fitter):
    """A wrapper for the LM-fitting of the chunks
    
    Args:
        model (:class:'SimpleModel'): A :class:'SimpleModel',
            which contains all the required submodels that are used in the
            fitting procedure.
    """
    
    def __init__(self, model):
        self.model = model

    @staticmethod
    def convert_params(params, to_lmfit=False, from_lmfit=False):
        """
        Convert between :class:'ParameterSet' and :class:'lmfit.Parameters'
        
        Args:
            params (:class:'ParameterSet' or :class:'lmfit.Parameters'): 
                An object containing fitting parameters.
            to_lmfit (Optional[bool]): Convert input to :class:lmfit.Parameters'.
            from_lmfit (Optional[bool]): Convert input to :class:'ParameterSet'.
        
        Returns:
            object: The converted parameter object.
        """
        if to_lmfit:
            lmfit_params = lmfit.Parameters()
            lmfit_params.add_many(*params.items())
            return lmfit_params
        elif from_lmfit:
            pyodine_params = ParameterSet()
            for k in params:
                pyodine_params[k] = params[k].value
            return pyodine_params
        else:
            raise ValueError('Either `to_lmfit` or `from_lmfit` must be true!')
    
    
    def fit_ostar(self, chunk, weight=None):
        """Convenience function for fitting with fixed velocity and no template
        
        """
        params = self.model.guess_params(chunk)
        lmfit_params = self.convert_params(params, to_lmfit=True)
        # Fix parameters
        # TODO: Let the model supply a list of parameters to fix for o-stars
        lmfit_params['velocity'].vary = False
        lmfit_params['tem_depth'].vary = False
        lmfit_params['iod_depth'].set(min=0.1)
        return self.fit(chunk, lmfit_params, weight=weight)
    

    def fit(self, chunk, lmfit_params, weight=None, chunk_ind=None, **kwargs):
        """Fit the chunk and return the best-fit result
        
        Args:
            chunk (:class:'Chunk'): The chunk to be modelled.
            lmfit_params (:class:'lmfit.Parameters'): The parameter object,
                defining starting values, limits, etc.
            weight (Optional[ndarray[npix]]): Pixel weights to use in the
                model evaluation.
            chunk_ind (Optional[int]): Index of the chunk, to grab the
                respective chunk from the template.
        
        Returns:
            :class:'LmfitResult': The best-fit result.
        """
        
        # Add-on: pixel weights in the fitting function, as used in dop code
        def func(lmfit_params, x, weight, chunk_ind):
            params = self.convert_params(lmfit_params, from_lmfit=True)
            if isinstance(weight, (list, np.ndarray)):
                return (self.model.eval(chunk, params, require=None, chunk_ind=chunk_ind) - chunk.flux) * np.sqrt(np.abs(weight))
            else:
                return self.model.eval(chunk, params, require=None, chunk_ind=chunk_ind) - chunk.flux

        # Carry out the fit
        try:
            # Make sure that the initial parameter guesses are consistent with
            # template and iodine atlas coverage
            params = self.convert_params(lmfit_params, from_lmfit=True)
            self.model.eval(chunk, params, require='full', chunk_ind=chunk_ind)
            # Carry out the fit
            # Add-on: pixel weights in the fitting function, as used in dop code
            lmfit_result = lmfit.minimize(func, lmfit_params, args=[chunk.pix, weight, chunk_ind], **kwargs)
            # Make sure that the fitted parameters are consistent with
            # template and iodine atlas coverage
            new_params = self.convert_params(lmfit_result.params, from_lmfit=True)
            self.model.eval(chunk, new_params, require='full', chunk_ind=chunk_ind)
            # Return output as LmfitResult object
            return self.LmfitResult(chunk, self.model, lmfit_result, chunk_ind=chunk_ind)
        except NoDataError:
            print('No Data!')
            return self.LmfitResult(chunk, self.model, None, chunk_ind=chunk_ind)


    class LmfitResult(FitResult):
        """Results from LmfitWrapper
        
        Args:
            chunk (:class:'Chunk'): The chunk which was modelled.
            model (:class:'SimpleModel'): The :class:'SimpleModel' object used 
                in the fitting procedure.
            lmfit_result (:class:'lmfit.MinimizerResult'): The best-fit results
                from the modelling. None if it failed.
            chunk_ind (Optional[int]): The index of the modelled chunk.
        """
        def __init__(self, chunk, model, lmfit_result, chunk_ind=None):
            self.chunk = chunk
            self.model = model
            self.lmfit_result = lmfit_result
            self.chunk_ind = chunk_ind

        @property
        def params(self):
            """Return a :class:'ParameterSet' with the fitted parameters"""
            if self.lmfit_result is not None:
                return LmfitWrapper.convert_params(self.lmfit_result.params, from_lmfit=True)
            else:
                return ParameterSet({p: np.NaN for p in self.model.all_param_names})

        @property
        def errors(self):
            """Return a dictionary of standard errors for the fitted parameters"""
            if self.lmfit_result is not None:
                lp = self.lmfit_result.params
                return {p: lp[p].stderr for p in lp}
            else:
                return {p: np.NaN for p in self.model.all_param_names}

        @property
        def init_params(self):
            """Return a dictionary of initial values"""
            if self.lmfit_result is not None:
                lp = self.lmfit_result.params
                params = ParameterSet()
                for n in list(self.lmfit_result.params.keys()):
                    if lp[n].vary is False:
                        params[n] = lp[n].value
                    else:
                        ii = self.lmfit_result.var_names.index(n)
                        params[n] = self.lmfit_result.init_vals[ii]
                return params
            else:
                return ParameterSet(
                    {p: np.NaN for p in self.model.all_param_names}
                )

        @property
        def report(self):
            """Return a fit report"""
            if self.lmfit_result is not None:
                return lmfit.fit_report(self.lmfit_result)
            else:
                return 'Chunk failed...'

        @property
        def redchi(self):
            """Return the red. Chi**2 of the fit"""
            if self.lmfit_result is not None:
                return self.lmfit_result.redchi
            else:
                return np.NaN

        @property
        def neval(self):
            """Return the number of evaluations of the fit"""
            if self.lmfit_result is not None:
                return self.lmfit_result.nfev
            else:
                return 0
        
        @property
        def medcounts(self):
            """Return the median counts of the chunk"""
            return np.median(self.chunk.flux)
    
    
    def fit_lsfs(self, lsf_model, params):
        """Fit the lsf model of this initialized fitter object to another
        lsf, defined by the input arguments.
        
        Args:
            lsf_model (:class:'LSFModel'): The LSF model to fit to.
            params (:class:'ParameterSet'): The LSF parameters to evaluate the
                supplied LSF model.
        
        Returns:
            :class:'ParameterSet': The best-fit parameters of the fit.
        """
        
        def fit_func(lmpars, x):
            pars = self.convert_params(lmpars, from_lmfit=True)
            return lsf_y - self.model.lsf_model.eval(x, pars)
        
        x = self.model.lsf_model.generate_x(self.model.osample_factor, self.model.conv_width)
        
        lsf_y = lsf_model.eval(x, params)        
        
        lmpars = self.convert_params(self.model.lsf_model.guess_params(0), to_lmfit=True)
    
        lmfit_result = lmfit.minimize(fit_func, lmpars, args=[x])
        
        return self.convert_params(lmfit_result.params, from_lmfit=True)
        
        