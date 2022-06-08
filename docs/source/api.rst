.. NAP documentation master file, created by
   sphinx-quickstart on Thu Jun  2 17:30:09 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

API reference
=============

Here will be described the functions of NAP, grouped by modules. 


Study class
-----------

The Study class contains information about a series of samples that are relevant to be studied together (replicates, change of only one experimental condition, etc).
This class is meant to be instanciated into objects.

.. note::

   Study class is used as an input for many NAP functions.

.. automodule:: dreem_nap.study
   :members:
   :undoc-members:
   :show-inheritance:


Plots
-----

Description of the functions shown in the :doc: gallery.

.. note:: 

   Most plots requires a Pandas dataframe, a specific sample or a study, and sometimes a (list of) construct(s).

.. automodule:: dreem_nap.plot
   :members:
   :undoc-members:
   :show-inheritance:


Data manipulation
-----------------

Just a few fonctions to navigate through your dataframes with ease, and to save data to csv files. 

.. automodule:: dreem_nap.data_manip
   :members:
   :undoc-members:
   :show-inheritance:


Data wrangler
-------------

Process DREEM's MutationHistogram pickle files, interface with the database and with local file, and get your dataframe clean.

.. automodule:: dreem_nap.data_wrangler
   :members:
   :undoc-members:
   :show-inheritance:


Database
--------

Connect, push and pull to the database. 

.. automodule:: dreem_nap.database
   :members:
   :undoc-members:
   :show-inheritance:


Utils
-----

A few low-level functions for the other modules.

.. automodule:: dreem_nap.utils
   :members:
   :undoc-members:
   :show-inheritance:


.. note::

   These modules are under active development.


NAP has its documentation hosted on Read the Docs.