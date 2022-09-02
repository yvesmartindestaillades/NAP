.. NAP documentation master file, created by
   sphinx-quickstart on Thu Jun  2 17:30:09 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

API reference
=============

Here will be described the functions of NAP, grouped by modules. 


Study class module
------------------

The Study class contains information about a series of samples that are relevant to be studied together (replicates, change of only one experimental condition, etc).
This class is meant to be instanciated into objects.

.. note::

   Study class is used as an input for many NAP functions.

.. automodule:: dreem_nap.study
   :members:
   :undoc-members:
   :show-inheritance:


Plots module
------------

Description of the functions shown in the :doc: gallery.

.. note:: 

   Most plots requires a Pandas dataframe, a specific sample or a study, and sometimes a (list of) construct(s).

.. automodule:: dreem_nap.plotter
   :members:
   :undoc-members:
   :show-inheritance:


Manipulator module
------------------------

Just a few fonctions to navigate through your dataframes with ease, and to save data to csv files. 

.. autofunction:: dreem_nap.manipulator.Manipulator.get_SCC

   .. code-block:: python

      df = salt.mani.get_SCC(samp='C6',
                              construct='9572', 
                              cols=['mut_rates','sequence','structure','cov_bases'],
                              base_type=['A','C'], 
                              index=list(range(40,50))) 
      df.to_csv('example.csv')

   ===== ======================= ======= ============ ========= 
   .     mut_rates               base    cov_bases    paired   
   ===== ======================= ======= ============ ========= 
   41    0.008445106805762544    C       1991.0       False    
   43    0.06855439642324888     C       1988.0       False    
   45    0.007948335817188276    C       1955.0       True     
   47    0.007451564828614009    A       1897.0       True     
   ===== ======================= ======= ============ ========= 



.. _loader_module:

Loader module
--------------------

Load and filter DREEM's MutationHistogram pickle files.

.. automodule:: dreem_nap.loader
   :members:
   :undoc-members:
   :show-inheritance:



Util module
------------

A few low-level functions for the other modules.

.. automodule:: dreem_nap.util
   :members:
   :undoc-members:
   :show-inheritance:


.. note::

   These modules are under active development.


NAP has its documentation hosted on Read the Docs.