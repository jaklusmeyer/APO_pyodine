import numpy as np
import logging
import sys

from pyodine.components import Chunk, ChunkArray


def simple(obs, width=91, padding=0, orders=None, chunks_per_order=None, pix_offset=0):
    """
    A simple chunking algorithm that splits the given orders of the
    observation in a number of fixed-size chunks (pixel space),
    leaving an equal amount of unused pixels in each end of the order.
    The chunks are not allowed to overlap, but if padding>0, the padded chunks
    will extend into their neighbour chunks (necessary for convolution).
    With chunks_per_order the number of chunks within an order can be 
    constrained; if it is not given, the maximum possible number of chunks 
    using the other parameters is generated.
    """

    # In case no orders were submitted, chunk all orders
    if orders is None:
        orders = slice(None)
    # In case only one order is submitted
    if type(orders) is int:
        orders = [orders]

    # Number of chunks per order
    max_chunks_per_order = int((obs.npix - 2 * padding) / width)
    if chunks_per_order is None:
        chunks_per_order = max_chunks_per_order
    elif chunks_per_order > max_chunks_per_order:
        raise ValueError('Cannot construct more than {} chunks per order with the given parameters!'.format(
                max_chunks_per_order))

    # Pixel offset of first chunk
    offset = int((obs.npix - 2 * padding - chunks_per_order * width) / 2) + int(padding)

    chunks = ChunkArray()
    for i in orders:
        for j in range(chunks_per_order):
            # Create a new chunk and calculate pixels
            pixels = offset + j * width + np.arange(width, dtype='int')
            chunk = Chunk(obs, i, pixels, padding)
            chunks.append(chunk)

    return chunks


def user_defined(obs, width=91, padding=0, orders=None, chunks_per_order=None, pix_offset0=None):
    """The standard chunking algorithm for template creation
    
    A simple chunking algorithm that splits the given orders of the
    observation in a number of fixed-size chunks (pixel space),
    leaving an equal amount of unused pixels in each end of the order.
    The chunks are not allowed to overlap, but if padding>0, the padded chunks
    will extend into their neighbour chunks (necessary for convolution).
    With chunks_per_order the number of chunks within an order can be 
    constrained; if it is not given, the maximum possible number of chunks 
    using the other parameters is generated.
    
    Args:
        obs (:class:'Observation'): The observation which will be chunked.
        width (Optional[int]): The chunk width in pixels. Defaults to 91.
        padding (Optional[int]): The padding width on either chunk side in
            pixels. Defaults to 0 (but you should make it bigger!).
        orders (Optional[int,list]): The order(s) which should be used. If None,
            all orders of the observation are used (default). (?!)
        chunks_per_order (Optional[int]): The number of chunks per order. If
            None, the maximum number of chunks of given width fitting into the
            order is used (default).
        pix_offset0 (Optional[int]): From which pixel to start the first chunk.
            If None, the chunking region is centered within the order (default).
    
    Return:
        :class:'ChunkArray': The created chunks.
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')

    # In case no orders were submitted, chunk all orders
    if orders is None:
        orders = obs.orders #slice(None)
    # In case only one order is submitted
    if type(orders) is int:
        orders = [orders]
    
    # Number of chunks per order
    max_chunks_per_order = int((obs.npix - 2 * padding) / width)
    if chunks_per_order is None:
        chunks_per_order = max_chunks_per_order
    elif chunks_per_order > max_chunks_per_order:
        logging.info('')
        logging.info('Warning! Max. nr. of chunks without cutting down is {}!'.format(
                max_chunks_per_order))

    # Pixel offset of first chunk
    if pix_offset0 is None:
        offset = int((obs.npix - 2 * padding - chunks_per_order * width) / 2) + int(padding)
    else:
        offset = pix_offset0

    chunks = ChunkArray()
    for i in orders:
        endpix2 = offset
        for j in range(chunks_per_order):
            # Create a new chunk and calculate pixels
            pixels = endpix2 + np.arange(width, dtype='int')#offset + j * width + np.arange(width, dtype='int')
            
            startpix = pixels[0]
            endpix = startpix + width
            # Make sure that chunks do not extend past the order edges
            if startpix < 0:
                logging.info('Startpixel (order {}, chunk {}): {}'.format(
                        i, j, startpix))
                logging.info('-> Correcting to 0.')
                startpix2 = 0
                width2 = width + startpix
            elif endpix > len(obs[i]):
                logging.info('Endpixel (order {}, chunk {}): {}'.format(
                        i, j, endpix))
                logging.info('-> Correcting to maximum pixel in order.')
                startpix2 = startpix
                width2 = width - (endpix - len(obs[i]))
            else:
                startpix2 = startpix
                width2 = width
            
            # Adapt the padding to not extend past the order edges
            if startpix2 - padding < 0:
                padding2 = startpix2
            else:
                padding2 = padding
            if startpix2 + width2 + padding2 > len(obs[i]):
                padding2 = len(obs[i]) - (startpix2 + width2)
            
            pixels2 = startpix2 + np.arange(width2)
            
            chunk = Chunk(obs, i, pixels2, padding2)
            chunks.append(chunk)
            endpix2 = pixels2[-1] + 1

    return chunks


def edge_to_edge(obs, width=91, padding=0, orders=None, chunks_per_order=None):
    """
        A variation of the "simple" algorithm. Instead of leaving unused pixels
        at the order edges, it allows chunks to overlap. The number of pixels
        defined by the "padding" keyword will still be left at the ends.

        Useful for generating templates.
    """

    # In case no orders were submitted, chunk all orders
    if orders is None:
        orders = slice(None)
    # In case only one order is submitted
    if type(orders) is int:
        orders = [orders]

    # Number of chunks per order
    if chunks_per_order is None:
        chunks_per_order = int(np.ceil((obs.npix - 2 * padding) / width))

    # Starting pixels of the chunks
    startpix = np.linspace(
        padding,
        obs.npix - padding - width,
        chunks_per_order,
        dtype='int'
    )

    chunks = ChunkArray()
    for i in orders:
        for j in range(chunks_per_order):
            # Create a new chunk and calculate pixels
            pixels = startpix[j] + np.arange(width, dtype='int')
            chunk = Chunk(obs, i, pixels, padding)
            chunks.append(chunk)

    return chunks


def wave_defined(obs, temp, width=91, padding=0, orders=None, order_correction=0,
                 delta_v=None):
    """The standard chunking algorithm for the observation modelling
    
    Chunk the observation to defined wavelength sections, corresponding
    to template chunks, using the relative barycentric velocity of 
    observation to template. This way chunks of a series of observations will
    always hold the same wavelength information (apart from the RV shift).
    
    Built similarly as the chunking algorithm in the dop code package.
    
    Args:
        obs (:class:'Observation'): The observation which will be chunked.
        temp (:class:'StellarTemplate_Chunked'): The corresponding deconvolved
            stellar template.
        width (Optional[int]): Desired width of the chunks, which needs to 
            correspond to template chunk width. If None, the template chunk
            width is chosen automatically (default).
        padding (Optional[int]): The padding width on either chunk side in
            pixels. Defaults to 0 (but you should make it bigger!).
        orders (Optional[int, list]): The order(s) which should be used. If None,
            the same orders as in the template are used (default).
        order_correction (Optional[int]): In case of a change of the zeroth 
            extracted order over time, this mismatch between template and 
            observation can be corrected (obs_orders + order_correction = 
            template_orders). Defaults to 0.
        delta_v (Optional[float]): If supplied, this will be used as relative 
            velocity between template and observation instead of the relative
            barycentric velocity. Defaults to None.
    
    Return:
        :class:'ChunkArray': The created chunks.
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    # Check if desired chunk width corresponds to template chunk width
    temp_width = temp[1].pix0 - temp[0].pix0
    if width == None:
        width = temp_width
    if width != temp_width:
        raise ValueError(
                'Desired chunk width does not correspond to template chunk width: {}'.format(
                        temp_width))
    c = 299792458.0
    
    # Template orders
    temp_orders = temp.orders_unique
    
    # In case no orders were submitted, chunk only the same orders as in template
    if orders is None:
        orders = temp_orders
    # In case only one order is submitted
    elif isinstance(orders, int):
        orders = [orders]
    
    chunks = ChunkArray()
    
    # Calculate barycentric velocity shift between template and observation
    # (or use the shift supplied as argument to this function)
    if isinstance(delta_v, float):
        init_dv = delta_v
    else:
        init_dv = (temp.bary_vel_corr - obs.bary_vel_corr)
    init_z  = init_dv / c
    logging.info('')
    logging.info('Barycentric redshift between template and observation: ')
    logging.info('v = {}, z = {}\n'.format(init_dv, init_z))
    
    # One order at a time
    for o in orders:
        order_ind = temp.get_order_indices(o)
        for i in order_ind:
            
            shft_wav = temp[i].w0 + init_z * temp[i].w0
            diff = np.abs(shft_wav - obs[o+order_correction].wave)
            pix_ind = np.argmin(np.abs(diff)) # pixel with closest wavelength
            
            startpix = round(pix_ind - width // 2)
            endpix = startpix + width
            # Make sure that chunks do not extend past the order edges
            if startpix < 0:
                logging.info('Startpixel (order {}, chunk {}): {}'.format(
                        o, i, startpix))
                logging.info('-> Correcting to 0.')
                startpix2 = 0
                width2 = width + startpix
            elif endpix > len(obs[o+order_correction]):
                logging.info('Endpixel (order {}, chunk {}): {}'.format(
                        o, i, endpix))
                logging.info('-> Correcting to maximum pixel in order.')
                startpix2 = startpix
                width2 = width - (endpix - len(obs[o+order_correction]))
            else:
                startpix2 = startpix
                width2 = width
            
            # Adapt the padding to not extend past the order edges
            if startpix2 - padding < 0:
                padding2 = startpix2
            else:
                padding2 = padding
            if startpix2 + width2 + padding2 > len(obs[o+order_correction]):
                padding2 = len(obs[o+order_correction]) - (startpix2 + width2)
            
            # Create a new chunk and calculate pixels
            pixels = startpix2 + np.arange(width2, dtype='int')
            chunk = Chunk(obs, o+order_correction, pixels, padding2)
            chunks.append(chunk)
    
    return chunks


def user_defined2(obs, wave_dict, padding=0):
    """An algorithm to create completely user-defined chunks, i.e. start and
    end wavelengths for all chunks are given by the user
    
    :param obs: The observation which will be chunked.
    :type obs: :class:`Observation`
    :param wave_dict: A dictionary with start and end wavelengths for all 
        chunks (arrays under 'start_wave' and 'end_wave').
    :type wave_dict: dict
    :param padding: The padding width on either chunk side in pixels. Defaults 
        to 0 (but you should make it bigger!).
    :type padding: int
    
    :return: The created chunks.
    :rtype: :class:`ChunkArray`
    """
    
    # Setup the logging if not existent yet
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(stream=sys.stdout, level=logging.INFO, 
                            format='%(message)s')
    
    # Unpack dictionary wavelengths
    start_wavelengths = np.array(wave_dict['start_wave'])
    end_wavelengths   = np.array(wave_dict['end_wave'])
    
    # Make sure that both arrays are of same size
    if len(start_wavelengths) != len(end_wavelengths):
        raise ValueError(
                'Different number of start and end wavelengths supplied: {} and {}!'.format(
                        len(start_wavelengths), len(end_wavelengths)))
    
    chunks = ChunkArray()
    for start_wave, end_wave in zip(start_wavelengths, end_wavelengths):
        # Return the order with best wavelength coverage from the observation
        order, coverage = obs.check_wavelength_range(start_wave, end_wave)
        
        # If no valid order could be returned, or wavelength range is not
        # fully covered, raise an error
        if order is None or coverage != 1.0:
            raise ValueError(
                    'No (complete) coverage for desired wavelength range: {} - {}'.format(
                            start_wave, end_wave))
        
        # Find the corresponding start and end pixels within the order
        start_pix = np.searchsorted(obs[order].wave, start_wave, side='right') - 1
        end_pix   = np.searchsorted(obs[order].wave, end_wave, side='left') + 1
        
        pixels = np.arange(start_pix, end_pix, dtype='int')
        
        # Adapt the padding to not extend past the order edges
        if start_pix - padding < 0:
            padding2 = start_pix
        else:
            padding2 = padding
        if end_pix + padding2 > len(obs[order]):
            padding2 = len(obs[order]) - end_pix
        
        chunk = Chunk(obs, order, pixels, padding2)
        chunks.append(chunk)

    return chunks