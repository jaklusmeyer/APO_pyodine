import numpy as np
from .base import StaticModel, ParameterSet


class LSFModel:
    """Abstract base class for the LSF models
    """
    @staticmethod
    def generate_x(osample_factor, conv_width=10.):
        """Generate a pixel vector to sample the LSF over
        
        Args:
            osample_factor (float): The oversample factor of the model.
            conv_width (Optional[float]): The number of pixels on either side.
                Defaults to 10.
        
        Return:
            ndarray[(2*int(conv_width)*int(osample_factor)+1)]: The evaluated
                pixel vector.
        """
        #return np.linspace(-4.0, 4.0, 8 * osample_factor + 1)
        return np.linspace(-conv_width, conv_width, 
                           int(2*conv_width) * int(osample_factor) + 1)


class SingleGaussian(LSFModel, StaticModel):
    """The LSF model of a Single Gaussian
    
    1 free parameter: The FWHM of the Gaussian.
    """
    param_names = ['fwhm']

    @staticmethod
    def eval(x, params):
        """Evaluate the LSF
        
        Args:
            x (ndarray): The pixel vector over which to evaluate the LSF.
            params (:class:'ParameterSet'): The LSF parameters.
        
        Return:
            ndarray: The normalized LSF.
        """
        # A single gaussian defined by its FWHM
        y = np.exp(-2.77258872223978123768 * x**2. / params['fwhm']**2)
        # Make sure that the sum equals one
        return y / np.sum(y)  # FIXME: Normalize to unit area?

    @staticmethod
    def guess_params(chunk):
        """Guess the LSF parameters for a given chunk
        
        At the moment, this returns just fixed values, independent of the chunk.
        
        Args:
            chunk (:class:'Chunk'): The chunk for which to make the guess.
        
        Return:
            :class:'ParameterSet': The guessed LSF parameters.
        """
        return ParameterSet(fwhm=2.0)  # FIXME: Make a better guess
    
    @staticmethod
    def name():
        """The name of the LSF as a string
        
        Return:
            str: The LSF name.
        """
        return 'SingleGaussian'


class SuperGaussian(LSFModel, StaticModel):
    """The LSF model of a Super Gaussian
    
    4 free parameters: Central sigma, central exponent, left and right Gaussian
        amplitudes.
    """
    param_names = ['sigma', 'exponent', 'left', 'right']

    @staticmethod
    def eval(x, params):
        # Include fixed satellite parameters
        a = np.array([params['left'], 1.0, params['right']])
        b = np.array([-1.0, 0.0, 1.0])
        c = np.array([2.5, params['sigma'], 2.5])
        n = np.array([2.0, params['exponent'], 2.0])

        # Supergauss function
        def func(x):
            xarr = np.transpose([x] * 3)
            f = np.sum(a * np.exp(-0.5 * (np.abs(xarr - b) / c) ** n), axis=1)
            # Replace negative values with zero
            f[np.where(f < 0.0)] = 0.0
            return f

        # Evaluate function
        y = func(x)
        
        # added for NaN-debugging
        y_sum = np.sum(y)
        if np.isnan(y_sum):
            print(y)
        #if y_sum==0:
        #    print('Sum of lsf is 0. Setting to 1e-4.')
        #    print(y)
        #    y_sum = 1e-4
        
        # Calculate centroid and re-center the LSF
        offset = np.sum(x * y) / y_sum #np.sum(y)
        y = func(x + offset)
        
        # added for NaN-debugging
        if any(np.isnan(y)):
            print('NaN value detected in un-normalized lsf function.')
            print('Sum of y: ', y_sum)
            print(params)
        y_sum = np.sum(y)
        #if y_sum==0:
        #    print('Sum of lsf is 0. Setting to 1e-4.')
        #    y_sum = 1e-4
        y = y / y_sum
        if any(np.isnan(y)):
            print('NaN value detected in lsf function.')
            print('Sum of y: ', y_sum)
            print(params)
        
        # Make sure that the sum equals one
        return y# / np.sum(y)

    @staticmethod
    def guess_params(chunk):
        return ParameterSet(
            sigma=1.0, exponent=1.9, left=0.2, right=0.2) # left, right 0.0
        # FIXME: Make a better qualified guess
    
    @staticmethod
    def name():
        return 'SuperGaussian'


class MultiGaussian_SONG(LSFModel, StaticModel):
    """The LSF model of a Multi Gaussian as used in SONG
    
    10 free parameters: Amplitudes of the Satellite Gaussians (5 left, 5 right).
    """
    param_names = [
        'left_5', 'left_4', 'left_3', 'left_2', 'left_1',
        'right1', 'right2', 'right3', 'right4', 'right5',
    ]

    @staticmethod
    def eval(x, params):
        # Temporary hack. Static methods can't see the rest of the class
        param_names = [
            'left_5', 'left_4', 'left_3', 'left_2', 'left_1',
            'right1', 'right2', 'right3', 'right4', 'right5',
        ]
        # convert input dict to list (nope, this is not pretty..)
        params = np.array([params[k] for k in param_names])

        # Set up parameter vectors, including central gaussian
        a = np.array([
            params[0], params[1], params[2], params[3], params[4],
            1.0,
            params[5], params[6], params[7], params[8], params[9],
        ])
        b = np.array([
            -2.6, -2.0, -1.6, -1.2, -0.7,
            0.0,
            0.7, 1.2, 1.6, 2.0, 2.6,
        ])
        c = np.array([
            0.7, 0.7, 0.7, 0.7, 0.7,
            0.4,
            0.7, 0.7, 0.7, 0.7, 0.7,
        ])
        n = 11

        # Multigauss function
        def func(x):
            xarr = np.repeat([x], n, axis=0)
            f = np.sum(a * np.exp(-0.5 * ((np.transpose(xarr) - b) / c)**2.), axis=1)
            #f[np.where(f < 0.0)] = 0.0
            return f

        # Evaluate function and find centroid
        y = func(x)
        
        # added for NaN-debugging
        y_sum = np.sum(y)
        #if y_sum==0:
        #    print('Sum of lsf is 0. Setting to 1e-4.')
        #    print(y)
        #    y_sum = 1e-4
        
        # Calculate centroid and re-center the LSF
        offset = np.sum(x * y) / y_sum #np.sum(y)  # TODO: Is this the correct way of weighting?
        y = func(x + offset)
        
        # added for NaN-debugging
        if any(np.isnan(y)):
            print('NaN value detected in un-normalized lsf function.')
        y_sum = np.sum(y)
        #if y_sum==0:
        #    print('Sum of lsf is 0. Setting to 1e-4.')
        #    y_sum = 1e-4
        
        y = y / y_sum
        if any(np.isnan(y)):
            print('NaN value detected in lsf function.')
            print('Sum of y: ', y_sum)
            print(params)
        
        return y

    @staticmethod
    def guess_params(chunk):
        # These are the median parameters from the cf of rs15.31
        # 0.4    0.0820768    0.0557928     0.167839     0.417223     0.453222
        #        0.400129     0.390609     0.138321    0.0599234    0.0588782
        return ParameterSet(
            left_5=0.1, left_4=0.2, left_3=0.3, left_2=0.5, left_1=0.7,
            right1=0.7, right2=0.5, right3=0.3, right4=0.2, right5=0.1
        )  # FIXME: Make a better guess
    
    @staticmethod
    def name():
        return 'MultiGaussian_SONG'


class MultiGaussian_SONGnew(LSFModel, StaticModel):
    """The LSF model of a Multi Gaussian as used in SONG (new)
    
    10 free parameters: Amplitudes of the Satellite Gaussians (5 left, 5 right).
    """
    param_names = [
        'left_5', 'left_4', 'left_3', 'left_2', 'left_1',
        'right1', 'right2', 'right3', 'right4', 'right5',
    ]

    @staticmethod
    def eval(x, params):
        # Temporary hack. Static methods can't see the rest of the class
        param_names = [
            'left_5', 'left_4', 'left_3', 'left_2', 'left_1',
            'right1', 'right2', 'right3', 'right4', 'right5',
        ]
        # convert input dict to list (nope, this is not pretty..)
        params = np.array([params[k] for k in param_names])

        # Set up parameter vectors, including central gaussian
        a = np.array([
            params[0], params[1], params[2], params[3], params[4],
            1.0,
            params[5], params[6], params[7], params[8], params[9],
        ])
        b = np.array([
            -2.9, -2.5, -1.9, -1.4, -1.0,
            0.0,
            1.0, 1.4, 1.9, 2.5, 2.9,
        ])
        c = np.array([
            0.9, 0.9, 0.9, 0.9, 0.9,
            0.6,
            0.9, 0.9, 0.9, 0.9, 0.9,
        ])
        n = 11

        # Multigauss function
        def func(x):
            xarr = np.repeat([x], n, axis=0)
            f = np.sum(a * np.exp(-0.5 * ((np.transpose(xarr) - b) / c)**2.), axis=1)
            f[np.where(f < 0.0)] = 0.0
            return f

        # Evaluate function and find centroid
        y = func(x)
        
        # added for NaN-debugging
        y_sum = np.sum(y)
        #if y_sum==0:
        #    print('Sum of lsf is 0. Setting to 1e-4.')
        #    print(y)
        #    y_sum = 1e-4
        
        # Calculate centroid and re-center the LSF
        offset = np.sum(x * y) / y_sum #np.sum(y)  # TODO: Is this the correct way of weighting?
        y = func(x + offset)
        
        # added for NaN-debugging
        if any(np.isnan(y)):
            print('NaN value detected in un-normalized lsf function.')
        y_sum = np.sum(y)
        #if y_sum==0:
        #    print('Sum of lsf is 0. Setting to 1e-4.')
        #    y_sum = 1e-4
        
        y = y / y_sum
        if any(np.isnan(y)):
            print('NaN value detected in lsf function.')
            print('Sum of y: ', y_sum)
            print(params)
        
        return y

    @staticmethod
    def guess_params(chunk):
        # These are the median parameters from the cf of rs15.31
        # 0.4    0.0820768    0.0557928     0.167839     0.417223     0.453222
        #        0.400129     0.390609     0.138321    0.0599234    0.0588782
        return ParameterSet(
            left_5=0.1, left_4=0.2, left_3=0.3, left_2=0.5, left_1=0.7,
            right1=0.7, right2=0.5, right3=0.3, right4=0.2, right5=0.1
        )  # FIXME: Make a better guess
    
    @staticmethod
    def name():
        return 'MultiGaussian_SONGnew'


class MultiGaussian(LSFModel, StaticModel):
    """The LSF model of a Multi Gaussian
    
    10 free parameters: Amplitudes of the Satellite Gaussians (5 left, 5 right).
    """
    param_names = [
        'left_5', 'left_4', 'left_3', 'left_2', 'left_1',
        'right1', 'right2', 'right3', 'right4', 'right5',
    ]

    @staticmethod
    def eval(x, params):
        # Temporary hack. Static methods can't see the rest of the class
        param_names = [
            'left_5', 'left_4', 'left_3', 'left_2', 'left_1',
            'right1', 'right2', 'right3', 'right4', 'right5',
        ]
        # convert input dict to list (nope, this is not pretty..)
        params = np.array([params[k] for k in param_names])

        # Set up parameter vectors, including central gaussian
        a = np.array([
            params[0], params[1], params[2], params[3], params[4],
            1.0,
            params[5], params[6], params[7], params[8], params[9],
        ])
        # In Butler 1996: Gaussians placed at 0.5 pixels apart
        # This is from the cf's
        b = np.array([
            -2.4, -2.1, -1.6, -1.1, -0.6,
            0.0,
            0.6, 1.1, 1.6, 2.1, 2.4
            ])
        c = np.array([
            0.3, 0.3, 0.3, 0.3, 0.3,
            0.4,
            0.3, 0.3, 0.3, 0.3, 0.3
            ])
        n = 11

        # Multigauss function
        def func(x):
            xarr = np.repeat([x], n, axis=0)
            f = np.sum(a * np.exp(-0.5 * ((np.transpose(xarr) - b) / c)**2.), axis=1)
            #f[np.where(f < 0.0)] = 0.0
            return f

        # Evaluate function and find centroid
        y = func(x)
        
        # added for NaN-debugging
        y_sum = np.sum(y)
        #if y_sum==0:
        #    print('Sum of lsf is 0. Setting to 1e-4.')
        #    print(y)
        #    y_sum = 1e-4
        
        # Calculate centroid and re-center the LSF
        offset = np.sum(x * y) / y_sum #np.sum(y)  # TODO: Is this the correct way of weighting?
        y = func(x + offset)
        
        # added for NaN-debugging
        if any(np.isnan(y)):
            print('NaN value detected in un-normalized lsf function.')
        y_sum = np.sum(y)
        #if y_sum==0:
        #    print('Sum of lsf is 0. Setting to 1e-4.')
        #    y_sum = 1e-4
        
        y = y / y_sum
        if any(np.isnan(y)):
            print('NaN value detected in lsf function.')
            print('Sum of y: ', y_sum)
            print(params)
        
        return y# / np.sum(y)  # FIXME: Normalize to unit area?

    @staticmethod
    def guess_params(chunk):
        # These are the median parameters from the cf of rs15.31
        # 0.4    0.0820768    0.0557928     0.167839     0.417223     0.453222
        #        0.400129     0.390609     0.138321    0.0599234    0.0588782
        return ParameterSet(
            left_5=0.1, left_4=0.2, left_3=0.3, left_2=0.5, left_1=0.7,
            right1=0.7, right2=0.5, right3=0.3, right4=0.2, right5=0.1
        )  # FIXME: Make a better guess
    
    @staticmethod
    def name():
        return 'MultiGaussian'


class MultiGaussian_Lick(LSFModel, StaticModel):
    """The LSF model of a Multi Gaussian as used in Lick (algorithm employed
    as in dop code)
    
    10 free parameters: Amplitudes of the Satellite Gaussians (5 left, 5 right).
    """
    param_names = [
        'left_5', 'left_4', 'left_3', 'left_2', 'left_1', 
        'right1', 'right2', 'right3', 'right4', 'right5'
        ]
    
    @staticmethod
    def eval(x, params):
        # Temporary hack. Static methods can't see the rest of the class
        param_names = [
            'left_5', 'left_4', 'left_3', 'left_2', 'left_1', 
            'right1', 'right2', 'right3', 'right4', 'right5'
            ]
        # convert input dict to list (nope, this is not pretty..)
        params = np.array([params[k] for k in param_names])

        # Set up parameter vectors, including central gaussian
        a = np.array([
            params[0], params[1], params[2], params[3], params[4],
            1.0,
            params[5], params[6], params[7], params[8], params[9]
            ])
        # In Butler 1996: Gaussians placed at 0.5 pixels apart
        # This is from the cf's
        b = np.array([
            -2.4, -2.1, -1.6, -1.1, -0.6,
            0.0,
            0.6, 1.1, 1.6, 2.1, 2.4
            ])
        c = np.array([
            0.3, 0.3, 0.3, 0.3, 0.3,
            0.4,
            0.3, 0.3, 0.3, 0.3, 0.3
            ])
        #n = 11

        # Gaussian function, computed as in Lick dop code
        def func(x, cntr=None):
            y = np.zeros(len(x))
            if cntr is None:
                cntr = 0.0
            cen = 0. - cntr
            # First central Gaussian
            cent_wid = c[5] * 5.
            if cent_wid < 2.:
                cent_wid = 2.  # define the central gaussian over this restricted pixel interval
            xx = np.where((x >= cen-cent_wid) & (x <= cen+cent_wid))
            y[xx] = np.exp(-0.5 * ((x[xx] - cen) / c[5])**2.)
            
            # Now surrounding Gaussians
            for i in range(len(a)):
                if i != 5:
                    cen = b[i] - cntr
                    gd_range = 5. * c[i]
                    xx = np.where((x >= cen-gd_range) & (x <= cen+gd_range))
                    y[xx] += a[i] * np.exp(-0.5 * ((x[xx] - cen) / c[i])**2.)
            # Normalize
            #dx = (x[-1] - x[0]) / len(x)
            #y[np.where(y < 0.0)] = 0.0
            
            return y / np.sum(y) #(dx * np.sum(y))

        # Evaluate function and find centroid
        y = func(x)
        # SHIFT TO CENTER - CHECK THIS - MAY WANT TO ELIMINATE!
        """
        fwhm = 0.5 * np.max(y)
        x2 = np.where((y >= fwhm) & (np.abs(x) < 6.)) # peak points +/- from center
        
        if len(x2[0]) >= 3:
            #print('>=3')
            dd = np.where(y[x2] == np.max(y[x2]))
            ndd = len(dd[0])
            if ndd <= 2:
                #print('ndd<=2')
                dd = dd[0][0]
            if ndd > 2:
                #print('ndd>2')
                dd = dd[0][int(ndd/2.)]
            cntr = x[x2[0][dd]]
            #print(cntr)
            if abs(cntr) >= 0.1 and abs(cntr) < 1.2:
                #print('Shift LSF.')
                y = func(x, cntr=cntr)
                #y = y / np.sum(np.abs(y))#(dx * np.sum(y))
                #print('New: ', x[np.argmax(y)])
        """
        return y

    @staticmethod
    def guess_params(chunk):
        # These are the median parameters from the cf of rs15.31
        # 0.4    0.0820768    0.0557928     0.167839     0.417223     0.453222
        #        0.400129     0.390609     0.138321    0.0599234    0.0588782
        return ParameterSet(
            left_5=0.1, left_4=0.2, left_3=0.3, left_2=0.4, left_1=0.5,
            right1=0.5, right2=0.4, right3=0.3, right4=0.2, right5=0.1
        )  # FIXME: Make a better guess
    
    @staticmethod
    def name():
        return 'MultiGaussian_Lick'


class FixedLSF(LSFModel, StaticModel):
    """The LSF model for a fixed LSF
    
    3 free parameters: Amplitude of the LSF, order and pixel0 of the respective
        LSF chunk.
    """
    param_names = [
        'amplitude', 'order', 'pixel0'
        ]
    
    @staticmethod
    def eval(x, params):
        """The parameter x serves here not as the x-vector over which to
        compute the LSF, but instead is the LSF itself!
        Small hack for the time being.
        """
        lsf_fixed = x[params['order'], params['pixel0']]
        return lsf_fixed * params['amplitude']
        

    @staticmethod
    def guess_params(chunk):
        return ParameterSet(
                amplitude=1., order=chunk.order, pixel0=chunk.abspix[0]
                )  # FIXME: Make a better guess
    
    @staticmethod
    def name():
        return 'FixedLSF'


model_index = {
        'SingleGaussian': SingleGaussian,
        'SuperGaussian': SuperGaussian,
        'MultiGaussian_SONG': MultiGaussian_SONG,
        'MultiGaussian_SONGnew': MultiGaussian_SONGnew,
        'MultiGaussian': MultiGaussian,
        'MultiGaussian_Lick': MultiGaussian_Lick,
        'FixedLSF': FixedLSF}


class LSF_Array:
    """A convenience class to enable modeling of a fixed LSF
    
    Needs to be supplied with a full lsf_array, and respective orders
    and pixels arrays describing the position of each LSF chunk.
    
    Args:
        lsf_array (ndarray[nr_chunks_total,nr_lsf_pix]): An array of 
            evaluated LSFs for all chunks in all orders.
        orders (ndarray[nr_chunks_total]): An array of order numbers for all
            chunks.
        pixels (ndarray[nr_chunks_total]): An array of pixel numbers for all
            chunks.
    """
    def __init__(self, lsf_array, orders, pixels):
        self.lsf_array = lsf_array
        self.orders = orders
        self.pixels = pixels
    
    def __getitem__(self, args):
        """The dedicated get-method
        
        Return the LSF of a desired chunk.
        
        Args:
            args (tuple): Should be a 2-entry tuple with desired order and
                pixel number.
        
        Return:
            ndarray[nr_lsf_pix]: An array containing the LSF of the desired
                chunk.
        """
        if len(args) != 2:
            raise IndexError('Two indices expected.')
        order, pixel = args
        try:
            return self.lsf_array[np.where((self.orders==order) & (self.pixels==pixel))[0][0]]
        except Exception as e:
            raise e