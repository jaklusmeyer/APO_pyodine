.. _new_utilities:

Setup your instrument's ``utilities`` module
============================================

In the development of **pyodine** we really wanted to make the software highly flexible - we have therefore tried to strictly separate all universally valid routines from code that might be different between instruments. Each instrument therefore comes with its own ``utilities`` module (e.g. ``utilities_song`` and ``utilities_lick``), which contains interchangeable code and instrument-specific parameters (see :ref:`overview_utilities`); so, when you adapt **pyodine** to your own project (let's call it *TANGO* for *Tricky Acronym of a New Great Observatory*), you should add your own module ``utilities_tango``. This is probably easiest done by

1. copy-pasting another ``utilities`` module (e.g. ``utilities_song`` to ``utilities_tango``), and then

2. adapting everything in the newly created module to your needs.

In the following, we assume you have done that and explain what might have to be changed in your new module.

.. _new_utilities_conf:

Adapt ``conf.py``
-----------------

The ``utilities_tango`` module contains a file ``conf.py``, which holds geographical information on the instrument and the pathname to the I2 atlas of the project (see also :ref:`overview_utilities_conf`). So, you should add *TANGO* to the supported instruments in the first dictionary, and the I2 atlas path for *TANGO* to the second dictionary (we assume here that you've placed the I2 atlas into the directory 'iodine_atlas' within the **pyodine** repository):
::

    my_instruments = {
                #other instruments,
                'tango': Instrument(
                                'TANGO project',
                                 latitude=0.0,
                                 longitude=0.0,
                                 altitude=0.0
                )
    }

    my_iodine_atlases = {
                #other I2 atlases,
                3: os.path.join(i2_dir_path, 'tango_i2.h5')
    }

.. _new_utilities_load_pyodine:

Adapt ``load_pyodine.py``
-------------------------

The file ``load_pyodine.py`` contains all the necessary code to correctly read data products of the instrument in question (see :ref:`overview_utilities_load_pyodine`). The first class :class:`IodineTemplate` specifies how to load the I2 atlas data, which depends on the file format that the data is stored in (more on that in :ref:`new_i2atlas_format`). But basically, if you have saved your *TANGO* I2 atlas in the same way as done for the other instruments, you do not have to change anything here.

The second class :class:`ObservationWrapper` is used to store the observations acquired with your instrument in memory and make them accessible for **pyodine** - i.e. both the actual spectral data of the (already extracted and wavelength calibrated) Echelle orders as well as supplementary information such as date and object name of an observation. The class contains a number of properties that should definitely be filled out:

* ``instrument`` should contain the initiliazed :class:`pyodine.components.Instrument` object of your project (*TANGO* in this case, as defined in ``conf.py``, see :ref:`new_utilities_conf`);

* ``star`` should contain an initialized :class:`pyodine.components.Star` class of the observed object (star name and coordinates);

* ``bary_date`` holds the barycentric date of the exposure, usually the photon-weighted midpoint;

* ``bary_vel_corr`` should be filled with an estimate of the barycentric velocity correction at the time of ``bary_date`` (this barycentric velocity is used to correctly place the chunks in the modelling, and a precision of :math:`\mathcal{O}(10 \,\mathrm{m/s})` suffices).

Usually, observations are saved in FITS format, in which case the spectral data is loaded from the FITS data unit, and the supplementary information from the FITS header. Both for Lick and SONG, we use several helper functions (also defined in ``load_pyodine.py``) to correctly extract the required data, so best is you have a look and adapt it as needed for your own instrument!

Adapt ``pyodine_parameters.py``
-------------------------------



Adapt ``timeseries_parameters.py``
----------------------------------


