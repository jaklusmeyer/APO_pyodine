.. pyodine documentation master file, created by
   sphinx-quickstart on Mon Oct 18 21:22:18 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyodine's documentation!
===================================

The **pyodine** package helps you to determine precise Radial Velocities (RVs) from stellar spectra, using the iodine (I2) cell method developed by `Butler et al. (1996) <https://ui.adsabs.harvard.edu/abs/1996PASP..108..500B/abstract/>`_. We aimed to build it in a flexible, modular and well-structured way, so as to make it (comparably) comprehensible and allow an easy and quick adaptation to different instruments.

Thus far, we have successfully tested **pyodine** on spectra from two different instruments, the Hamilton spectrograph at Lick observatory and the SONG spectrograph site on Tenerife.

Check out the :doc:`overview` section for further information, including how to :ref:`install <installation>` the project.

.. note::
   This project is under active development.

.. toctree::
   :maxdepth: 2
   :caption: Quick start

   installation
   overview
   code_structure

.. toctree::
   :maxdepth: 1
   :caption: Tutorial

   tutorial/preparation
   tutorial/template
   tutorial/observation
   tutorial/velocities


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
