
=========
Examples
=========


Load data
=========

.. code-block:: python
    :caption: Example code

    from dreem_nap.study import Study

    salt = Study.from_dict({'name': 'salt',
                            'description': 'Just an example study',
                            'samples': ['A1', 'B2', 'B3'], 
                            'label': "Salt concentration [nM]", 
                            'conditions': [60, 120, 180]})

    salt.load_df_from_local_files(path_to_data= '/path/to/data/',
                                  min_cov_bases= 1000)
    # Show the dataframe
    salt.get_df().head()


.. note::

    Part :ref:`loading_data` teaches you how to:

    * read studies from a csv file.
    * filter constructs based on their base coverage on specific bases.
    
    



Plot data
=========

A single plot
*************

.. code-block:: python

    salt.plot.mut_histogram(samp='B6',
                            construct='3114',
                            cluster=0,
                            index=list(range(19,80)),
                            base_paired=True,
                            structure='structure')

.. image:: img/mut_histogram.png
    :align: center

Multiple plots
**************

.. code-block:: python

    from dreem_nap import util
    sub_lib = 'MS2' # say that you only want to plot this sub-library
    mpl.use('agg')  # use this to avoid display issues
    df = salt.get_df()
    for samp in salt.samples:
        for construct in salt.constructs:
            for cluster in df[(df['samp']==samp)&(df['construct']==construct)]['cluster'].unique():
                if df[(df['samp']==samp)&(df['construct']==construct)&(df['cluster']==cluster)]['sub-library'].iloc[0] == sub_lib:
                    salt.plot.mut_histogram(samp=samp,construct=construct, cluster=cluster)
                    # save the figure and closes it
                    util.save_fig(path_to_figs+'/'+salt.name+'/'+samp+'_'+construct+'_mut_histogram.png') 


.. note::

    Take a look at the plots gallery in :ref:`plots`.            

Download data
=============

This part shows how to filter data from dataframes and to turn it into a csv file. 

.. note::

    Learn more about downloading to csv functions at :ref:`to_csv`


Download several columns of a single sample-construct-cluster
*****************************************************************

.. code-block:: python

    df = salt.mani.get_SCC(samp='C6',
                           construct='9572', 
                           cols=['mut_rates','sequence','structure','cov_bases'],
                           base_type=['A','C'], 
                           index=list(range(40,50))) 
    df.to_csv('example.csv')

===== ======================= ======= ============ ========= 
 .    mut_rates               base    cov_bases    paired   
===== ======================= ======= ============ ========= 
41    0.008445106805762544    C       1991.0       False    
43    0.06855439642324888     C       1988.0       False    
45    0.007948335817188276    C       1955.0       True     
47    0.007451564828614009    A       1897.0       True     
===== ======================= ======= ============ ========= 


Download a single column of several sample-construct-clusters
**************************************************************

.. code-block:: python

    df = study.mani.get_col_across_constructs(samp=470, 
                                              col='mut_rates',
                                              index=list(range(40,50))) 
    df.to_csv('example.csv')



====== ======================= ======================= ======================= ====================== ======================== ======================= ======================= ======================= ======================= ======================= 
.       40                      41                      42                      43                     44                       45                      46                      47                      48                      49                     
====== ======================= ======================= ======================= ====================== ======================== ======================= ======================= ======================= ======================= ======================= 
323    0.001987083954297069    0.008445106805762544    0.003974167908594138    0.06855439642324888    0.00894187779433681      0.007948335817188276    0.003477396920019871    0.007451564828614009    0.006951340615690168    0.011420059582919563   
478    0.009218163195629908    0.016729259132809832    0.0013656538067599864   0.048822123591669514   0.0027313076135199728    0.05769887333560942     0.04848071013997952     0.0013656538067599864   0.006828269033799932    0.006145442130419939   
619    0.0028622540250447226   0.008586762075134168    0.006797853309481216    0.0611587982832618     0.00536480686695279      0.010014306151645207    0.006437768240343348    0.009298998569384835    0.002861230329041488    0.004291845493562232   
834    0.0007651109410864575   0.008416220351951033    0.0007651109410864575   0.06006120887528692    0.14957918898240244      0.010328997704667177    0.061208875286916604    0.011859219586840091    0.020275439938791124    0.0971690895179801     
====== ======================= ======================= ======================= ====================== ======================== ======================= ======================= ======================= ======================= ======================= 



