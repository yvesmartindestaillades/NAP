.. _loading_data:

================
How to load data
================

Create a study
==============

First, you have to create a study. This is a python object that contains data that you want to analyze. A study analyzes a set of samples. A specific experimental condition ``label`` varies between samples, and is associated with the values of ``conditions``.

.. code-block:: python

    salt = Study.from_dict({'name': 'salt',
                            'description': 'Change the Na concentration', 
                            'samples': ['A6', 'B6', 'C6', 'D6', 'E6'], 
                            'label': 'Na quantity [M]', 
                            'conditions': [0.15, 0.3, 0.6, 1.0, 1.2]})

.. note::
    
    If you want to use several sets of samples, create several studies.


You also can create a study from a csv file.

::

    studies = Study.load_studies('studies.csv')
    study = studies['salt']


.. list-table:: studies.csv
   :widths: 25 75 25 25 50
   :header-rows: 1

   * - name
     - description
     - samples
     - conditions
     - label
   * - salt
     - Change the Na concentration
     - A6
     - 0.15
     - Na quantity [M]
   * - salt
     - Change the Na concentration
     - B6
     - 0.3
     - Na quantity [M]
   * - salt
     - Change the Na concentration
     - C6
     - 0.6
     - Na quantity [M]
   * - salt
     - Change the Na concentration
     - D6
     - 1
     - Na quantity [M]
   * - salt
     - Change the Na concentration
     - E6
     - 1.2
     - Na quantity [M]


Load and filter the data
========================

Load
****


Then, you want to load your data. 
Your data has to be organized like this:

::

    path/to/data/
        |-- A6.p
        |-- B6.p
        |-- C6.p
        |-- D6.p
        |-- E6.p

.. note::

    To get the output of DREEM in this format, use Herschlag lab's ``dreem_herschlag`` python package to run DREEM.

Filter by study or by sample
****************************

When loading your data, you need to give a **minimal base-coverage**.
If a construct has any of its bases covered by less than this value, this construct will be discarded.

.. code-block:: python

    salt.load_df_from_local_files(path_to_data= '/path/to/data/', 
                                  min_cov_bases=1000, 
                                  filter_by = 'study')


If you use the argument ``filter_by='study'``, NAP will filter out constructs that didn't reach the right amount of coverage **in all samples**.
If you use the argument ``filter_by='sample'``, NAP will filter out constructs that didn't reach the right amount of coverage **in the sample**.

Advanced filtering
******************

You can use more advanced filtering by specifying which bases you want to reach the minimum coverage.
See :ref:`filtering_data` for more details.

.. code-block:: python

    salt.load_df_from_local_files(path_to_data= '/path/to/data/', 
        min_cov_bases=1000, 
        filter_by = 'study',
        index = 'roi', # only filter by these bases. Could be a list of indexes (ex: [13, 14, .., 34]), a unique sub-sequence (ex: 'ATCTAGGTTAC') or 'all' (default).
        base_type = ['A', 'C'], # only filter by A, C
        base_paired = True, # only filter by paired bases. Can be False (only non-paired bases) or None (all bases, default).
        structure = 'structure_DMS') # the structure prediction to use if you want to filter by base-pairing. Default is None.

.. note::

  Advanced filtering technique takes approximately 50 times as long as the basic filtering technique.


Use a config file (optional)
============================

Write a file ``config.yaml`` that contains the parameters you want to use.

.. code-block:: yaml
    :caption: config.yaml

    path_to_data: /Users/ymdt/src/data/Gabe
    path_to_studies: /Users/ymdt/src/data/Gabe/studies.csv
    min_cov_bases: 1000
    filter_by: study
    index: all
    base_type: ['A', 'C']
    base_paired: True
    structure: structure_DMS

Load your config file using:

.. code-block:: python

    with open(path+'config.yml', 'r') as ymlfile:
        cfg = yaml.safe_load(ymlfile)


Summary of the code (using a config file)
=========================================

.. code-block:: python

    from dreem_nap.study import Study

    with open(path+'config.yml', 'r') as ymlfile:
      cfg = yaml.safe_load(ymlfile)

    salt = Study.load_studies(cfg['path_to_studies'])['salt']

    salt.load_df_from_local_files(path_to_data= '/path/to/data/',
                                  min_cov_bases= cfg['min_cov_bases'],
                                  filter_by = cfg['filter_by'],
                                  index = cfg['index'],
                                  base_type = cfg['base_type'],
                                  base_paired = cfg['base_paired'],
                                  structure = cfg['structure'])

    # Show the dataframe
    salt.get_df().head()


.. note::    

    NAP loads 1-indexed data from DREEM and returns 0-indexed data for arrays such as ``cov_bases``, ``mut_bases``, and more.


.. note::    

    If the mut_histogram object loaded from DREEM doesn't contain a cluster attribute (i.e doesn't use Expectation-Maximization algorithm), NAP will define cluster=0 and use it by default.



