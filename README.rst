MACA-ETr
========

Download `MACA <https://climate.northwestknowledge.net/MACA/index.php>`__ (Multivariate Adaptive Constructed Analogs) downscaled climate data and estimate ASCE standardized reference evapotranspiration.

Given a point location ``MACA-ETr`` simplifies the process of downloading downscaled climate variables from the MACA dataset via a simple Python API. Data that can be retrieved (and used to calculate reference ET) includes global climate models from the Coupled Model Inter-Comparison Project 5 (20 models) whose output is downscaled using the MACA methodologies. The MACA downscaling model is trained from both the Livneh et al. 2014 dataset and the gridMET dataset, output from both can be downloaded using the ``MACA-ETr`` package. The historical (1950-2005) and future predictions (rcp4.5 and rcp8.5 emissions trajectories from 2006-2099) datasets are included in MACA. The ``MACA-ETr`` package can also be used to calculate ASCE short and tall reference ET time series using data downloaded from MACA at a given location within the conterminous United States.  

Documentation 
-------------

Under development

Installation
------------

You may install the dependencies using the conda virtual environment (recommended), the environment file can be downloaded `here <https://raw.githubusercontent.com/WSWUP/MACA-ETr/master/environment.yml>`__ and installed and activated by

.. code-block:: bash

   conda env create -f environment.yml
   conda activate macaetr

Once activated install with PIP:

.. code-block:: bash

   pip install macaetr


