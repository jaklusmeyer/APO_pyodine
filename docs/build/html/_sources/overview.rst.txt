.. _overview:

Overview
========

Create a deconvolved stellar template
-------------------------------------

To create a deconvolved stellar template, you can use the ``pyodine_create_templates.create_template()`` function:

.. autofunction:: pyodine_create_templates.create_template


Model a single observation
--------------------------

After you have created a deconvolved stellar template, you can use it to model an observation of the same star, obtained **WITH** the I2 cell in the light path. This is done with the function ``pyodine_model_observations.model_single_observation()``:

.. autofunction:: pyodine_model_observations.model_single_observation


Model multiple observations
---------------------------

In many cases you will have more than one single observation of the same star that you want to model. Then you can take advantage of Python's parallelizing capabilities and use the function ``pyodine_model_observations.model_multi_observations()`` to model the observations on several cores at the same time:

.. autofunction:: pyodine_model_observations.model_multi_observations


Combine chunk velocities to RV timeseries
-----------------------------------------

When you've modelled a number of observations and saved the fit results, you can use the function ``pyodine_combine_vels.combine_velocity_results()`` to combine the chunk velocities from the observations to a RV timeseries of the star:

.. autofunction:: pyodine_combine_vels.combine_velocity_results
