============
Introduction
============

Overview
********

NAP works through two steps:

#. :ref:`Data processing`. Merge the output of DREEM with the output of RNAstructure (or any other RNA structure prediction tool), process it and push it to the database.
#. :ref:`Data analysis`. Define a study, pull the corresponding data from the database, and generate plots and csv files.

Both steps will be described here.


Data processing
***************

This section explains how NAP builds .json dataframes from sources such as DREEM and a RNA structure prediction software.
NAP's module ``data_wrangler`` merges the sources, filters out the unvalid data, and pushes the data to database using NAP's module ``database``.

.. note::
    The functions are built for specific formats, indicated below. 
    You must respect this formatting to have NAP working.

.. blockdiag::
    :align: center    
    :caption: Diagram 1: Data processing
    :width: 1200

    blockdiag {
        group {
            orientation = portrait;
            color = lightgray;
            DREEM -> 'NAP data wrangler' -> database ;
            RNAstructure -> 'NAP data wrangler';
            }
    }





DREEM
-----

Output of `Prof. Yesselman's DREEM <https://github.com/jyesselm/dreem>`_, under the  `pickle format <https://docs.python.org/3/library/pickle.html>`_.
One DREEM pickle file corresponds to one sample.

DREEM files must be stored using the following tree structure. 
Data wrangler will use ``path_to_dreem`` to read the files.

::

    model: path_to_dreem                   ex:  data/DREEM
           |-- sample_name                     |-- A1
               |-- mutation_histos.p               |-- mutation_histos.p
                                               |-- A2
                                                   |-- mutation_histos.p
                                               |-- ...

        
.. note::

    A **sample** corresponds to the content of a physical tube, read by a sequencer and processed by DREEM.
    
.. note::
    
    A **sample** comes with specific experimental conditions and a sub-library of RNA constructs.  


RNAstructure 
------------

Output of RNAstructure, or any RNA structure prediction software, under a csv format. 
The csv file has specific column names. 
Each row corresponds to a construct.

**Columns names**
    * ``construct``: name of the constructs of the sample's sub-library.
    * ``full_sequence``: sequence of the entire RNA molecule.
    * ``roi_sequence``: sequence of the ROI only.
    * ``full_deltaG``: predicted deltaG for the entire RNA molecule.
    * ``roi_deltaG``: predicted deltaG for the ROI only.
    * ``full_structure``: predicted structure for the entire RNA molecule.
    * ``roi_structure_comparison``: Comparison between the pairing-prediction of the entire RNA molecule and the pairing-prediction of the ROI only, for the ROI bases. String of '0' and '1', of same length as ROI sequence. '0' means that both predicted structures have the same pairing state for the corresponding base. '1' means that the predicted structures have diverging pairing states for this base.
    * ``roi_start_index``: index of the first base of the ROI. Index starts with a 0.
    * ``roi_end_index``: index of the last base of the ROI. Index starts with a 0.
    * ``flank``: flank.
    * ``sub-library``: name of the sub-library.

.. note::
    
    ROI corresponds to Region of Interest.


Data wrangler
-------------

Turns DREEM and RNAstructure into a .json format sample by sample, filters out invalid constructs, and pushes the sample to the database.
Every function of NAP's module data wrangler is described on page :ref:`data_wrangler <data wrangler>`.

A construct in a sample is considered valid only if every base of the ROI has a base coverage above ``min_bases_cov``.

The sample's json format structure is the following:
::

    sample
    |-- construct 1
        |-- full_sequence
        |-- roi_sequence
        |-- ...
    |-- ...
    |-- construct N
        |-- ...     


Database
--------

The database is hosted on Google Firebase.


Sample code
-----------

    *"Un bon croquis vaut mieux qu'un long discours."* (*A good sketch is worth more than a long speech.*) - Napoléon Bonaparte

Let's show a code example.





.. _diag2:

Data analysis
*************



.. blockdiag::
   :align: center
   :caption: Diagram 2: Data Analysis
   :width: 1200

   blockdiag {
       group {
        orientation = portrait;
        color = lightgray;
        database -> 'NAP plot \n NAP data manip' -> plots ;
        studies -> 'NAP plot \n NAP data manip' -> csv;
        }
    }




