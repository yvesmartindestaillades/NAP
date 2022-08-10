========
Overview
========

Introduction
============

NAP works through two steps:

#. :ref:`Data processing`. Merge the output of DREEM with the output of RNAstructure (or any other RNA structure prediction tool), process it and push it to the database.
#. :ref:`Data analysis`. Define a study, pull the corresponding data from the database, and generate plots and csv files.



.. blockdiag::
    :align: center    
    :width: 400

    blockdiag {
        group {
            orientation = landscape;
            color = lightgray;
            'NAP data wrangler' -> database -> 'NAP plot';
            }
    }


Both steps will be described here.




Terminology
===========

    *Semantically speaking, what I love is agreeing on semantics, not debating them* - Matthew Allan

Let's define a few terms for a good understanding of this page.

:DREEM:
    **DREEM** is an experimental method used to identify unpaired bases of RNA strands.


:Sample:

    A **sample**  is a specific set of experimental procedures applied to a **library**.
    

:Library:

    A **library** corresponds to a set of **constructs**.

    A **sub-library** is a subset of a library.
    

:Construct:

    A **construct** is a DNA/RNA/protein molecule with a specific (usually artificial) sequence, which can have attributes like a name, structures, functions, etc.    
    Several samples can have the same sub-library - and therefore the same constructs. 


:Sample-construct:

    A **sample-construct** corresponds to one construct of a sample.
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

:DREEM pickle file:

    **DREEM pickle file** is the output of DREEM pipeline.
    It consists in a class call ``MutationHistogram`` that contains one **sample**.
    
    ``MutationHistogram`` is compressed using the Python module ``pickle``.





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
One DREEM pickle file contains the data of one sample.

This data needs to be read by ``data_wrangler``.

.. note::

    All you need to know is how to give ``data_wrangler`` a dictionary ``pickles`` structured such that `pickles[sample] = [path_to_the_corresponding_pickle_file]`

To make DREEM's files easy, we suggest the following tree structure: 

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


And to generate the dictionary, use the following code:
::

    >>> path_to_dreem = 'data/DREEM'
    >>> samples_list = ['A1','A2','B3']
    >>> pickles = {sample: f"{path_to_dreem}/{sample}/mutation_histos.p" for sample in samples_list}
    >>> print(pickles)
    {'A1': 'data/DREEM/A1/mutation_histos.p', 'A2': 'data/DREEM/A2/mutation_histos.p', 'B3': 'data/DREEM/B3/mutation_histos.p'}


Just a bit of code to illustrate how ``data_wrangler`` will use the ``pickles`` dictionary.
::
    
    >>> import pandas as pd
    >>> import pickle
    >>> from dreem_nap.data_wrangler import mhs2dict
    >>> 
    >>> samples_list = ['A1']
    >>> pickles = {sample: f"data/DREEM/{sample}/mutation_histos.p" for sample in samples_list}
    >>> for pick in pickles:
    ...     mhs = pickle.load(open(pickles[pick], "rb"))
    ...     df_sample = pd.DataFrame.from_dict(mhs2dict(mhs, drop_attribute = ['structure','_MutationHistogram__bases','sequence']),
    ...             orient='index').rename(columns={'name':'construct'})
    ...     print(df_sample.head())
    
    [5 rows x 19 columns]
      construct data_type  num_reads  num_aligned  ...                                        mod_bases_T skips_low_mapq skips_short_read skips_too_many_muts
    1         1       DMS          7            0  ...  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...              7                0                   0
    2         2       DMS         89            6  ...  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...             83                0                   0
    3         3       DMS         11            0  ...  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...             11                0                   0
    4         4       DMS        138            1  ...  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...            137                0                   0
    5         5       DMS          5            1  ...  (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...              4                0                   0


.. _intro_RNAstructure:

RNAstructure 
************

Output of RNAstructure, or any RNA structure prediction software, under a csv format. 
Your csv file must respect the names detailled below.
Each row corresponds to a sample-construct.
Each column corresponds to an attribute. 

**Columns names**
    * ``construct``: (str) name of this construct.
    * ``sequence``: (str) sequence of the entire RNA molecule.
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
            |-- sequence: "ACCGACTACTATC"  # Attribute from RNAstructure.
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

Check out :ref:`data processing sample code <data_processing_sample_code>`.



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
            |-- sequence
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



.. _data_processing_sample_code:

Sample code
***********

    *"Un bon croquis vaut mieux qu'un long discours."* (*A good sketch is worth more than a long speech.*) - Napoléon Bonaparte

For this example, we will use the example shown in `the getting started branch <https://github.com/yvesmartindestaillades/dreem_nap/tree/getting_started>`_ 

::

    >>> import pandas as pd
    >>> from dreem_nap import data_wrangler
    >>> import json
    >>> 
    >>> ## DREEM
    >>> # List the files that you want to process and create your pickles dict
    >>> samples_list = ['A1', 'A2','B3']
    >>> pickles = {sample: f"data/DREEM/{sample}/mutation_histos.p" for sample in samples_list}
    >>> 
    >>> ## RNA-STRUCTURE
    >>> # Indicate where is your RNAstructure file
    >>> RNAstructureFile = 'data/RNAstructureFile.csv'
    >>> 
    >>> ## DATA-WRANGLER
    >>> # Define what is the min base coverage values that you tolerate
    >>> min_bases_cov = 1000
    >>> 
    >>> ## DATABASE
    >>> # Select your root folder for the database 
    >>> folder = 'my_project_1/tutorial'
    >>> 
    >>> # Load Firebase credentials file 
    >>> firebase_credentials_file = 'data/credentials_firebase.json'
    >>> with open(firebase_credentials_file) as file:
    ...     firebase_credentials = json.load(file)
    ... 
    >>> ## PROCESS DATA
    >>> # Process your pickles files and push them to Firebase!
    >>> data_wrangler.push_samples_to_firebase(pickles = pickles,
    ...                     RNAstructureFile = RNAstructureFile,
    ...                     firebase_credentials = firebase_credentials,
    ...                     min_bases_cov = min_bases_cov, 
    ...                     folder=folder)
    Push pickles to firebase!
    A1 A2 B3




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




