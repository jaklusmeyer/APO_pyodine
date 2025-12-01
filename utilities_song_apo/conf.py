from pyodine.components import Instrument

import os

i2_dir_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../iodine_atlas')

# Used by get_instrument below...  # TODO: Make a settings file
my_instruments = {
    'song_1': Instrument(
            'SONG Hertzsprung spectrograph (Tenerife)',
            latitude=28.2983,
            longitude=-16.5094,  # East longitude
            altitude=2400.0
    ),
    'song_2': Instrument(
            'SONG China spectrograph (Delingha)',
            latitude=37.378001,
            longitude=97.73167,  # East longitude
            altitude=3200.    # From https://arxiv.org/pdf/1602.00838
    ),
    'lick': Instrument(
            'Hamilton spectrograph (Lick observatory)',
            latitude=37.34139,
            longitude=238.35722,
            altitude=1283.
    ),
    'waltz': Instrument(
            'Waltz spectrograph (LSW Heidelberg)',
            latitude=49.398611,
            longitude=8.720833,
            altitude=560.
    ),
    'song_3': Instrument(
    'SONG USA spectrograph (Apache Point Observatory)',
    latitude=32.7803,
    longitude=-105.8203,  # West longitude as negative
    altitude=2788  # Approximate elevation in meters
)
}
#2_dir_path_apo='/Users/jaklusmeyer/code/NMSU/SONG/'
# List of iodine atlas locations ##jessica fix 15OCT2024 - I2_r.h5 file has been masked so negative values are zero
my_iodine_atlases = {
            1: os.path.join(i2_dir_path, 'song_iodine_cell_01_65C.h5'),
                2: os.path.join(i2_dir_path, 'Fischer_Cell_May2022_downsampled3.h5'),
                3: os.path.join(i2_dir_path, 'APO_I2_cleaned_normalized.h5'),
                4: os.path.join(i2_dir_path, 'APO_I2_cleaned_normalized_shifted.h5')
}
                #4: os.path.join(i2_dir_path_apo, 'Paul_iodine_spectrum_APOSONGTEST_normalized.h5'),
                #5: os.path.join(i2_dir_path_apo, 'I2072109_003_r_testing.h5')