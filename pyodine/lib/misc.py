import numpy as np
from scipy.interpolate import interp1d
from scipy.special import erfinv
import os
import json
import logging
import logging.config

_c = 299792458  # m/s


def printLog(filename=None, *args):
    """Function to print to the terminal and at the same time to a specified file
    
    Args:
        filename (str): The name of the file to print to.
        *args (str): String argument(s) to print.
    """
    print(*args)
    if isinstance(filename, str) and filename != '':
        try:
            with open(filename, 'a') as f:
                print(*args, file=f)
        except Exception as e:
            print(e)


def setup_logging(config_file=None, level=logging.INFO, error_log=None,
                  info_log=None, quiet=False):
    """Setup logging configuration
    https://fangpenlin.com/posts/2012/08/26/good-logging-practice-in-python/
    """
    
    log_handlers = []
    
    # Is the configuration file (json format) there?
    if isinstance(config_file, str) and os.path.exists(config_file):
        try:
            # Try and load the configuration dictionary
            with open(config_file, 'rt') as f:
                config = json.load(f)
            
            # If you want errors logged
            if isinstance(error_log, str):
                
                # Create directory structure if non-existent yet
                error_log_dir = os.path.dirname(error_log)
                if error_log_dir != '' and not os.path.exists(error_log_dir):
                    os.makedirs(error_log_dir)
                
                config['handlers']['error_file_handler']['filename'] = error_log
                log_handlers.append('error_file_handler')
            else:
                del config['handlers']['error_file_handler']
            
            # If you want info logged
            if isinstance(info_log, str):
                
                # Create directory structure if non-existent yet
                info_log_dir = os.path.dirname(info_log)
                if info_log_dir != '' and not os.path.exists(info_log_dir):
                    os.makedirs(info_log_dir)
                
                config['handlers']['info_file_handler']['filename'] = info_log
                config['handlers']['info_file_handler']['level'] = level
                log_handlers.append('info_file_handler')
            else:
                del config['handlers']['info_file_handler']
            
            # If you want info printed or not 
            # (if not, errors will still be printed!)
            if quiet:
                #del config['handlers']['console']
                config['handlers']['console']['level'] = logging.ERROR
            else:
                config['handlers']['console']['level'] = level
            log_handlers.append('console')
            
            config['root']['handlers'] = log_handlers
            config['root']['level'] = level
            
            logging.config.dictConfig(config)
            
        except Exception as e:
            logging.basicConfig(level=level)
            logging.error('Logger could not be configured from config file', 
                          exc_info=True)
    else:
        if not quiet:
            logging.basicConfig(level=level)
        else:
            logging.basicConfig(level=logging.ERROR)
        logging.warning('No config file supplied or config file not found.')


def return_existing_files(filenames):
    """Check the input filenames and only return those that actually exist
    
    :param filenames: The full pathname(s) of the input file(s).
    :type filenames: str, list, ndarray, tuple
    
    :return: A list of existing files.
    :rtype: list
    :return: A list of non-existing files.
    :rtype: list
    """
    
    if isinstance(filenames, str):
        filenames = [filenames]
    
    if isinstance(filenames, (list,tuple,np.ndarray)):
        existing_files = [f for f in filenames if os.path.isfile(f)]
        bad_files = [f for f in filenames if not os.path.isfile(f)]
    else:
        raise ValueError('No files supplied! (Either of str, list, ndarray, tuple)')
    
    return existing_files, bad_files


def findwave(multiorder, wavelength):
    for i, spec in enumerate(multiorder):
        if spec.wave[0] < wavelength < spec.wave[-1]:
            return i
    raise ValueError('Could not find a corresponding order for wavelength %s Å' % wavelength)


# def sp_resample(deltav, wave):
#     """
#         Return a vector of wavelengths between `wl_start` and `wl_end` with
#         constant velocity spacing `deltav` (m/s)
#     """
#     npt = int(np.log(wave[-1] / wave[0]) / np.log(1.0 + deltav / _c))
#     return wave[0] * (1.0 + deltav / _c)**np.arange(npt, dtype='float')


def osample(x, factor):
    """Linear oversampling of a vector
    
    Args:
        x (ndarray[nr_pix_in]): Input vector.
        factor (int): Oversampling factor.
    
    Return:
        ndarray[nr_pix_out]: Resulting oversampled vector.
    """
    n = int(np.round(factor * (len(x) - 1) + 1))
    return np.linspace(x[0], x[-1], n)


def normgauss(x, fwhm):
    """A normalized Gaussian, defined by its FWHM
    
    Args:
        x (ndarray[nr_pix]): Input vector to sample over.
        fwhm (float): Full-width half-maximum of the Gaussian.
    
    Return:
        ndarray[nr_pix]: The normalized Gaussian.
    """
    y = np.exp(-2.77258872223978123768 * x**2 / fwhm**2)
    # Make sure that the sum equals one
    return y / np.sum(y)


def rebin(wold, sold, wnew):
    """Interpolates OR integrates a spectrum onto a new wavelength scale, depending
    on whether number of pixels per angstrom increases or decreases. Integration
    is effectively done analytically under a cubic spline fit to old spectrum.

    Ported to from rebin.pro (IDL) to Python by Frank Grundahl (FG).
    Original program written by Jeff Valenti.

    IDL Edit History:
    ; 10-Oct-90 JAV Create.
    ; 22-Sep-91 JAV Translated from IDL to ANA.
    ; 27-Aug-93 JAV Fixed bug in endpoint check: the "or" was essentially an "and".
    ; 26-Aug-94 JAV Made endpoint check less restrictive so that identical old and
    ;       new endpoints are now allowed. Switched to new Solaris library
    ;       in call_external.
    ; Nov01 DAF eliminated call_external code; now use internal idl fspline
    ; 2008: FG replaced fspline with spline
    
    Args:
        wold (ndarray[nr_pix_in]): Input wavelength vector.
        sold (ndarray[nr_pix_in]): Input spectrum to be binned.
        wnew (ndarray[nr_pix_out]): New wavelength vector to bin to.
    
    Return:
        ndarray[nr_pix_out]: Newly binned spectrum.
    """

    def idl_rebin(a, shape):
        sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] // shape[1]
        return a.reshape(sh).mean(-1).mean(1)

    # Determine spectrum attributes.
    nold   = np.int32(len(wold))            # Number of old points
    nnew   = np.int32(len(wnew))            # Number of new points
    psold  = (wold[nold - 1] - wold[0]) / (nold - 1)  # Old pixel scale
    psnew  = (wnew[nnew - 1] - wnew[0]) / (nnew - 1)  # New pixel scale

    # Verify that new wavelength scale is a subset of old wavelength scale.
    if (wnew[0] < wold[0]) or (wnew[nnew - 1] > wold[nold - 1]):
        logging.warning('New wavelength scale not subset of old.')

    # Select integration or interpolation depending on change in dispersion.

    if psnew < psold:

        # Pixel scale decreased ie, finer pixels
        # Interpolating onto new wavelength scale.
        dum  = interp1d(wold, sold)  # dum  = interp1d( wold, sold, kind='cubic' ) # Very slow it seems.
        snew = dum(wnew)

    else:

        # Pixel scale increased ie more coarse
        # Integration under cubic spline - changed to interpolation.

        xfac = np.int32(psnew / psold + 0.5) 	# pixel scale expansion factor

        # Construct another wavelength scale (W) with a pixel scale close to that of
        # the old wavelength scale (Wold), but with the additional constraint that
        # every XFac pixels in W will exactly fill a pixel in the new wavelength
        # scale (Wnew). Optimized for XFac < Nnew.

        dw   = 0.5 * (wnew[2:] - wnew[:-2])        # Local pixel scale

        pre  = np.float(2.0 * dw[0] - dw[1])
        post = np.float(2.0 * dw[nnew - 3] - dw[nnew - 4])

        dw   = np.append(dw[::-1], pre)[::-1]
        dw   = np.append(dw, post)
        w    = np.zeros((nnew, xfac), dtype='float')

        # Loop thru subpixels
        for i in range(0, xfac):
            w[:, i] = wnew + dw * (np.float(2 * i + 1) / (2 * xfac) - 0.5)   # pixel centers in W

        nig  = nnew * xfac			# Elements in interpolation grid
        w    = np.reshape(w, nig)   # Make into 1-D

        # Interpolate old spectrum (Sold) onto wavelength scale W to make S. Then
        # sum every XFac pixels in S to make a single pixel in the new spectrum
        # (Snew). Equivalent to integrating under cubic spline through Sold.

        # dum    = interp1d( wold, sold, kind='cubic' ) # Very slow!
        # fill_value in interp1d added to deal with w-values just outside the interpolation range
        dum    = interp1d(wold, sold, fill_value="extrapolate")
        s      = dum(w)
        s      = s / xfac				# take average in each pixel
        sdummy = s.reshape(nnew, xfac)
        snew   = xfac * idl_rebin(sdummy, [nnew, 1])
        snew   = np.reshape(snew, nnew)

    return snew


# Chauvenet criterion, as implemented in dop code
def chauvenet_criterion(residuals, iterate=True):
    """Chauvenet's criterion, as implemented in dop code (from the Kgiant server
    of the Landessternwarte Heidelberg):
    Find elements that lie too far away from the others.
    Updated: nan-values in the residuals are immediately marked as bad ones.
    
    Args:
        residuals (ndarray[nr_pix]): Input vector with residuals to analyze.
        iterate (Optional[bool]): Iteratively continue to throw out elements?
            Defaults to True.
    
    Return:
        ndarray[nr_pix]: A mask of same length as input residuals array, with 
            ones where the criterion was passed, and zeros where it failed.
        ndarray[nr_good]: An array with the indices where mask is True.
        ndarray[nr_bad]: An array with the indices where mask is False.
    """
    # WE ONLY CARE ABOUT FINITE VALUES...
    fin = np.where(np.isfinite(residuals))
    N = len(fin[0])
    if N == 0:
        raise ValueError('No finite values in the residuals!')
    # CHECK FOR CASE WHERE THERE'S ONLY ONE DATUM...
    # MAKE USER THINK ABOUT DREADFUL MISTAKE...
    if N == 1:
        raise ValueError('Chauvenets criterion cant be applied to only one datum!')
    
    #mean = np.mean(residuals[fin])
    #rms  = np.std(residuals[fin])
    mean = np.nanmean(residuals)
    rms  = np.nanstd(residuals)
    
    # FIND INVERSE ERROR FUNCTION OF (1 - 0.5/N)...
    # MAKE A MASK OF POINTS THAT PASSED CHAUVENET'S CRITERION...
    #mask = np.abs(residuals[fin]-mean) < 1.05*rms*erfinv(1. - 0.5/N)
    #mask = np.abs(residuals[fin]-mean) < 1.4142135623730950488*rms*erfinv(1. - 0.5/N)
    mask = np.abs(residuals-mean) < 1.4142135623730950488*rms*erfinv(1. - 0.5/N)

    # DO YOU REALLY WANT TO ITERATE...
    # BEVINGTON AND TAYLOR SAY YOU SHOULDN'T...
    if iterate == True:
        indx = np.where(mask == True)
        if len(indx[0]) == 0:
            raise ValueError('All data have failed Chauvenets criterion.',
                             'Check that your model is appropriate for these data.')

        if len(indx[0]) < N:
            #iter_mask, iter_mask_true, iter_mask_false = chauvenet_criterion(residuals[fin][indx])
            iter_mask, iter_mask_true, iter_mask_false = chauvenet_criterion(residuals[indx])
            mask[indx] = mask[indx] & iter_mask

    # RETURN THE INDICES WHERE YOU'VE PASSED CHAUVENET'S CRITERION... 
    return mask, np.where(mask == True), np.where(mask == False)


# Smooth over given chunk LSFs
def smooth_lsf(chunk_arr, pixel_avg, order_avg, order_dist, fit_results, redchi2=None, 
               osample=None, lsf_conv_width=None):
    """Smooth over all chunk LSFs, with given 'radii' in orders & order pixels. 
    LSFs are weighted by red. Chi2 values of modeled chunks.
    Implemented with great parallels to the dop IDL code dop_psf_smooth.pro.
    
    Args:
        chunk_arr (:class:'ChunkArray'): An array of chunks.
        pixel_avg (int): Pixels to smooth over in dispersion direction.
        order_avg (int): Orders to smooth over in cross-dispersion direction.
        order_dist (int): Approximate pixel distance between orders.
        fit_results (:class:'LmfitResult'): The fit result which holds the LSFs 
            to smooth over and red. Chi2s.
        redchi2 (Optional[ndarray[nr_chunks]]): An array of red. Chi2 values
            to use instead of the ones from the fit result.
        osample (Optional[int]): Oversampling to use; if not given, use the one 
            from the fit results model.
        lsf_conv_width (Optional[int]): Number of pixels to evaluate the LSF on
            (towards either side).
    
    Return:
        ndarray[nr_chunks, nr_pix]: An array with the smoothed LSFs for all
            chunks.
    """
    # First we set up an easy chunk numpy-array with orders and center pixels
    # to find relevant smoothing indices for each chunk later
    chunk_arr_simple = []
    chunk_lsfs = []
    for i, chunk in enumerate(chunk_arr):
        chunk_arr_simple += [[chunk.order, chunk.abspix[int(len(chunk)/2.)]]]
        if lsf_conv_width is None:
            lsf_conv_width = fit_results[0].model.conv_width
        # Also we already evaluate all individual LSFs here, then we can later
        # simply pick them
        if osample is None:
            x_lsf, lsf = fit_results[i].fitted_lsf
        else:
            x_lsf = fit_results[i].model.lsf_model.generate_x(osample, conv_width=lsf_conv_width)
            lsf = fit_results[i].model.lsf_model.eval(x_lsf, fit_results[i].params.filter('lsf'))
        chunk_lsfs.append(lsf)
        
    chunk_arr_simple = np.array(chunk_arr_simple)
    chunk_lsfs = np.array(chunk_lsfs)
    
    # Now check if redchi2 is given. If not, use the redchi2 from fit_results
    if redchi2 is None or len(redchi2) is not len(chunk_arr_simple):
        redchi2 = np.zeros(len(chunk_arr_simple))
        for i, result in enumerate(fit_results):
            redchi2[i] = fit_results[i].redchi
    else:
        redchi2 = np.array(redchi2)
    
    # Now loop over chunks
    lsfs_smoothed = []
    for i, chunk in enumerate(chunk_arr):
        # These are the chunk indices to smooth over
        xsm = np.where( (chunk_arr_simple[:,0] <= chunk.order + order_avg) &
                        (chunk_arr_simple[:,0] >= chunk.order - order_avg) &
                        (chunk_arr_simple[:,1] <= chunk.abspix[int(len(chunk)/2.)] + pixel_avg) &
                        (chunk_arr_simple[:,1] >= chunk.abspix[int(len(chunk)/2.)] - pixel_avg))
        #print(xsm)
        
        if len(xsm[0]) == 0:
            logging.warning('Chunk {}: No chunks to smooth over. Using input chunk.'.format(i))
            # Use lsf from fit result for that chunk
            lsfs_smoothed.append(chunk_lsfs[i])
        
        else:
            # Compute chunk weights for order averaging
            ord_sep = np.abs(chunk_arr_simple[xsm[0],0] - chunk.order) * order_dist
            ord_wt = 1. / np.sqrt(1. * ord_sep/len(chunk) + 1.)
            pix_sep = np.abs(chunk_arr_simple[xsm[0],1] - chunk.abspix[int(len(chunk)/2.)])
            pix_wt = 1. / np.sqrt(pix_sep/len(chunk) + 1.)
            #print(ord_sep, ord_wt, pix_sep, pix_wt)
            wt = 1. / redchi2[xsm[0]]
            wt = ord_wt * pix_wt * wt
            #print(redchi2[xsm[0]])
            #print(wt)
            # Check for nans
            nan_ind = np.argwhere(np.isnan(wt))
            if len(nan_ind) > 0:
                #print('Chunk {}, subchunk {}: nan wt'.format(i, xsm[0][nan_ind]))
                wt[nan_ind] = 0.
            
            # Now get all the lsfs
            lsf_array = []
            for j, xs in enumerate(xsm[0]):
                lsf = chunk_lsfs[xs]
                xneg = np.where(lsf < 0.0)
                if len(xneg[0]) > 0:
                    lsf[xneg[0]] = 0.0
                xnan = np.argwhere(np.isnan(lsf))
                if len(xnan) > 0:
                    print('Chunk {}, lsf {}: nan'.format(i, xs))
                else:
                    # No extra shifting to center here
                    lsf_array.append(lsf * wt[j])
            
            lsf_array = np.array(lsf_array)        
            lsf_av = np.median(lsf_array, axis=0)
            """
            if osample is None:
                lsf_av = fit_results[i].model.osample_factor * lsf_av / np.sum(lsf_av)
            else:
                lsf_av = osample * lsf_av / np.sum(lsf_av)"""
            lsf_av = lsf_av / np.sum(lsf_av) # we use a different normalization -> change that?!
            lsfs_smoothed.append(lsf_av)\
    
    lsfs_smoothed = np.array(lsfs_smoothed)
    
    return lsfs_smoothed


def smooth_parameters_over_orders(parameters, par_name, chunks, 
                                   nr_orders, nr_chunks_order, deg=2):
    """For a given parameter (par_name) in a list of :class:'ParameterSet'
    objects, fit a polynomial of given degree over the central chunk pixels 
    within each order and return the evaluated results.
    
    Args:
        parameters (list): A list of :class:'ParameterSet' objects.
        par_name (str): Parameter key to smooth.
        chunks (:class:'ChunkArray'): The chunks of the observation.
        nr_orders (int): The number of orders.
        nr_chunks_order (int): The number of chunks in each order.
        deg (Optional[int]): Degree of the polynomial (default: 2).
    
    Return:
        ndarray[nr_chunks]: A flattened array with the fit results of all chunks.
    """
    
    pfits = []
    for o in range(nr_orders):
        par_data = np.array([parameters[i][par_name] \
                             for i in range(o*nr_chunks_order,(o+1)*nr_chunks_order)])
        ch_pix = np.array([chunks[i].abspix[0] + (len(chunks[i]) // 2) \
                           for i in range(o*nr_chunks_order,(o+1)*nr_chunks_order)])
        
        # fit_polynomial returns fitted y-values and coefficients - only use the first
        pfits.append(fit_polynomial(ch_pix, par_data, deg=deg)[0])
        
    return np.array(pfits).flatten()


def smooth_fitresult_over_orders(fit_result, par_name, deg=2):
    """For a given parameter (par_name) in a fit_result, fit a polynomial of 
    given degree over the central chunk pixels within each order and 
    return the evaluated results.
    
    Args:
        fit_result (list): A list of :class:'LmfitResult' objects from the fit.
        par_name (str): Parameter key to smooth.
        deg (Optional[int]): Degree of the polynomial (default: 2).
    
    Return:
        ndarray[nr_chunks]: A flattened array with the fit results of all chunks.
    """
    pfits = []
    for o in range(fit_result[0].chunk.order, fit_result[-1].chunk.order+1):
        par_result_o = []
        ch_pix_o = []
        for i, res in enumerate(fit_result):
            if res.chunk.order == o:
                par_result_o.append(res.params[par_name])
                ch_pix_o.append(res.chunk.abspix[0] + (len(res.chunk) // 2))
        
        par_result_o = np.array(par_result_o)
        ch_pix_o = np.array(ch_pix_o)
        
        # fit_polynomial returns fitted y-values and coefficients - only use the first
        pfits.append(fit_polynomial(ch_pix_o, par_result_o, deg=deg)[0])
    
    return np.array(pfits).flatten()


def fit_polynomial(x_data, y_data, deg=2):
    """Fit a polynomial to data
    
    This routine masks out NaNs from the fit, and returns the input data
    if something goes wrong.
    
    Args:
        x_data (np.ndarray, list): The input x-data.
        y_data (np.ndarray, list): The input y-data to be modelled.
        deg (Optional[int]): Degree of the polynomial (default: 2).
    
    Return:
        ndarray: The fitted polynomial data, or the input y-data if something
            went wrong.
    """
    try:
        if isinstance(x_data, list):
            x_data = np.array(x_data)
        elif isinstance(y_data, list):
            y_data = np.array(y_data)
        
        idx = np.where(np.isfinite(y_data))
        n = len(x_data[idx])
        
        pfit, stats = np.polynomial.Polynomial.fit(
                x_data[idx], y_data[idx], deg, full=True, window=(0, n), domain=(0, n))
        return pfit(x_data), pfit.coef
    except Exception as e:
        logging.error('Polynomial fitting failed. Falling back to input values', 
                      exc_info=True)
        return y_data, None


# Weights/error calculation
def analytic_chunk_weights(chunk_flux, wave_intercept, wave_slope):
    """Calculate analytic weight of a chunk, based on its spectral content
    
    Implemented with great parallels to the dop IDL code dsst_m.pro, of the part
    based on equations in Bulter&Marcy (1996): sigma(mean) = 1./sqrt(sum(1/sig^2)).
    
    Args:
        chunk_flux (ndarray[nr_pix]): Flux values of the pixels in a chunk.
        wave_intercept (float): The wavelength intercept w0 of that chunk.
        wave_slope (float): The wavelength dispersion w1 of that chunk.
    
    Return:
        float: The weight of the chunk, representing its RV content.
    """
    eps = np.sqrt(chunk_flux[:-1])
    # slope: dI/d(pix)
    slope_int_per_pix = chunk_flux[1:] - chunk_flux[0:-1]
    # slope in real intensity per m/s
    slope_int_per_ms  = slope_int_per_pix * (wave_intercept / (_c * wave_slope))
    # Eqn. 5 in error write up
    weight = np.sum((slope_int_per_ms / eps)**2)
    
    return weight

