============
Introduction
============

Overview
========

NAP works through two steps:

#. :ref:`Data processing`. Merge the output of DREEM with the output of RNAstructure (or any other RNA structure prediction tool), process it and push it to the database.
#. :ref:`Data analysis`. Define a study, pull the corresponding data from the database, and generate plots and csv files.

Both steps will be described here.


Terminology
===========

Let's define a few terms for a good understanding of this page.

:Sample:

    A **sample** corresponds to the content of an experimental tube, read by a sequencer and processed by DREEM.
    This tube contained the RNA sequences defined in a  **sub-library**, and was placed in certain **experimental conditions**.  


:Sub-library:

    A **sub-library** corresponds to a set of RNA sequences, each of them named with a **construct**.


:Construct:

    A **construct** is the name of a sequence in a **sub-library**, and corresponds by extension to a set of **attributes** defined in the sub-library.
    For example, the **ROI** of the sequence or the RNAstructure-predicted deltaG.

    Several samples can have the same sub-library - and therefore the same constructs. 


:Sample-construct:

    A **sample-construct** corresponds to a specific construct in a specific sample.
    Each sample-construct has a set of **attributes**.
    
    A sample-construct can be nicknamed a construct if the corresponding sample is mentioned beforehand.


:Attributes:

    **Attributes** are a set of experimental results and properties related to a **sample-construct**.
    DREEM, RNAstructure and NAP provide attributes.
    
    Examples: ``ROI_sequence``, ``mut_bases``, ``cov_bases``. 

:Study:

    A **study** is a set of **samples** that are relevant to be studied together, typically because of their experimental conditions.
    For example, a study can be a set of replicates, or the same library with temperature variating, etc.

:ROI:

    **ROI** is the acronym for Region of Interest.
    It corresponds to a sub-sequence in the RNA sequence that we want to study.



Data processing
===============

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
*****

Output of `Yesselman and Rouskin's DREEM <https://github.com/jyesselm/dreem>`_, under the  `pickle format <https://docs.python.org/3/library/pickle.html>`_.
One DREEM pickle file corresponds to one sample.

DREEM files must be stored using the following tree structure. 
Data wrangler will use ``path_to_dreem`` to read the files.

::

    model:
    path_to_dreem                   
        |-- sample_name                     
            |-- mutation_histos.p              
                                               
    ex: 

    data/DREEM
        |-- A1
            |-- mutation_histos.p
        |-- A2
            |-- mutation_histos.p
        |-- ...


.. _intro_RNAstructure:

Example code:
::
    import pandas, pickle, 
    pickles = {pickle: f"../data/DREEM/{pickle}/mutation_histos.p" for pickle in pickles_list}
    for pick in pickles:
        mhs = pickle.load(open(pickles[pick], "rb"))
        df_sample = pd.DataFrame.from_dict(mhs2dict(mhs, drop_attribute = ['structure','_MutationHistogram__bases','sequence']),
                orient='index').rename(columns={'name':'construct'})
        print(df_sample.head())


RNAstructure 
************

Output of RNAstructure, or any RNA structure prediction software, under a csv format. 
Your csv file must respect the names detailled below.
Each row corresponds to a sample-construct.
Each column corresponds to an attribute. 

**Columns names**
    * ``construct``: (str) name of this construct.
    * ``full_sequence``: (str) sequence of the entire RNA molecule.
    * ``roi_sequence``: (str) sequence of the ROI only.
    * ``full_deltaG``: (float) predicted deltaG for the entire RNA molecule.
    * ``roi_deltaG``: (float) predicted deltaG for the ROI only.
    * ``full_structure``: (str) predicted structure for the entire RNA molecule.
    * ``roi_structure_comparison``: (str) comparison between the pairing-prediction of the entire RNA molecule and the pairing-prediction of the ROI only, for the ROI bases. String of '0' and '1', of same length as ROI sequence. '0' means that both predicted structures have the same pairing state for the corresponding base. '1' means that the predicted structures have diverging pairing states for this base.
    * ``roi_start_index``: (int) index of the first base of the ROI. Index starts with a 0.
    * ``roi_end_index``: (int) index of the last base of the ROI. Index starts with a 0.
    * ``flank``: (str) flank.
    * ``sub-library``: (str) name of the sub-library.



Data wrangler
*************

NAP's module data wrangler turns DREEM and RNAstructure into a .json format sample by sample, filters out invalid sample-constructs, and pushes the sample to the database.

Every function of data wrangler is described on page :ref:`data wrangler module <data_wrangler_module>`.

Merging DREEM and RNAstructure file
...................................

For each sample, the merge between DREEM and RNAstructure file is done w.r.t their respective ``construct`` column.
The fit is inner-typed, which means that each construct must be on both files. 


The data structure of a sample is the following:

::

    |-- a_sample
        |-- a_construct
            |-- full_sequence: "ACCGACTACTATC"  # Attribute from RNAstructure.
            |-- roi_sequence: "ACTACT"
            |-- ...
            |-- cov_bases: [0, 1769, 1795, ... ,1814, 1815, 1821] # Attribute from DREEM.
            |--
            |-- min_bases_cov: 1000 # Attribute from NAP
            |--


A more complete visualisation of the data structure can be found on :ref:`database section <intro_database_structure>`.

The columns of the merged dataset corresponds to the sample-constructs attributes. They are the following:

**Columns of the dataset**
    * Every column of :ref:`RNA structure file <intro_RNAstructure>`.
    * ``num_reads``: number of reads for this construct.
    * ``num_aligned``: (int) number of reads correctly aligned, that we will use for the analysis.
    * ``start`` : (int) beginning of the index for all list[int] type attributes. Default is 1, in which case you should start reading list[int]-typed attributes such as ``info_bases`` starting from the 2nd element.
    * ``end`` : (int) beginning of the index for all list[int] type attributes. 
    * ``num_of_mutations``: (list[int]) count of how many bases mutated n times. [4, 5, 1, 0] means that 4 bases didn't mutate, 5 bases mutated once, 1 base mutated twice, and no base mutated 3 times.
    * ``mut_bases`` : (list[int]) for each base, count of mutations.
    * ``info_bases`` : (list[int]) for each base, number of valid reads. 
    * ``del_bases`` : (list[int]) for each base, count of deletions.
    * ``ins_bases`` :(list[int])  for each base, count of inserts. 
    * ``cov_bases`` : (list[int]) for each base, the base-coverage.
    * ``mod_bases_A`` : (list[int]) for each base, the number of times that it mutated to a A base.
    * ``mod_bases_C`` : (list[int]) for each base, the number of times that it mutated to a C base.
    * ``mod_bases_G`` : (list[int]) for each base, the number of times that it mutated to a G base.
    * ``mod_bases_T`` : (list[int]) for each base, the number of times that it mutated to a T base.
    * ``skips_low_mapq`` : (int) number of reads that that we don't use because the map score is too low (default is below 15)
    * ``skips_short_read`` : (int) number of reads that we don't use because they are too short.
    * ``skips_too_many_muts`` : (int) number of reads that that we don't use because they have so many mutations, and therefore we have low confidence.
    * ``cov_bases_roi`` : (int) worst base coverage among the bases of the ROI.
    * ``cov_bases_sec_half`` : (int) worst base coverage among the bases of the second half of the sequence.


.. note::

    If every sample has the same constructs, RNAstructure information will be redundant between the sample-constructs.


Filtering out invalid constructs
................................

Valid construct:
    A sample-construct is considered valid only if every base of its ROI has a base coverage above ``min_bases_cov``.

Unvalid sample-constructs are filtered out, such that each sample loaded into the database contain only constructs that passed the filter.


Pushing samples to the database
...............................

Data wrangler connects to the database, and pushes the data sample by sample onto the database. 
Data is organised by folders and subfolders.

If when pushing a sample, a file of the same name exists in the same folder, it will be overwritten.

Most of the information is on the :ref:`section database <intro_database>`.


Sample code
...........

    *"Un bon croquis vaut mieux qu'un long discours."* (*A good sketch is worth more than a long speech.*) - Napoléon Bonaparte

Let's show a code example.





.. _intro_database:

Database
********

.. note::

    NAP's database is a module used by NAP's data_wrangler, but rarely used by the user itself.
    You only need to know how the credentials works and how the database is structured.   


.. _intro_database_structure:

Structure
.........

The database is hosted on Google Firebase. It uses the .json format.

A database root folder is called a `folder`, and corresponds to a project, a user, a version, etc.
In a folder is stored the data of a project, using the following structure:

::

    my_project_1
    |-- sample_1
        |-- construct 1
            |-- full_sequence
            |-- roi_sequence
            |-- ...
        |-- ...
        |-- construct N
            |-- ...     
    |-- sample_2
        |-- ...
    |-- ...

It is possible to create different folders and subfolders using ``/``, such as: ``my_project_2/user_1/version v2.0``:

::

    my_project_1
    |-- version_v1.0
        |-- ...    
    |-- version_v2.0    
        |-- ...    
    my_project_2
    |-- user_1
        |-- version_v1.0
            |-- ...    
        |-- version_v2.0    
            |-- ...    
    |-- user_2    
        |-- ...      
    ...


Credentials
...........

The :ref:`database.connect() <database_module>` function uses credentials to access the database, under the form of a dictionary.
Please email `yves@martin.yt <mailto:yves@martin.yt>`_ to get this your credentials.
You can also create your own database for free on `Google Firebase <https://firebase.google.com/>`_.


Example:
::

    >>> from dreem_nap import database
    >>> import json
    >>> # Firebase credentials file
    >>> firebase_credentials_file = 'data/credentials_firebase.json'
    >>> with open(firebase_credentials_file) as file:
    >>>     firebase_credentials = json.load(file)
    >>> # Give credentials to connect to firebase
    >>> database.connect(firebase_credentials)
    Initiated connection to Firebase!
    >>> database.connect(firebase_credentials)
    Re-used the previous Firebase connection










.. _diag2:

Data analysis
=============



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




