=========================
Dataframe columns content
=========================


DREEM's output
===============

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


Library.csv content
===================

See the attributes of the library at `library.csv attributes from dreem_herschlag <https://github.com/yvesmartindestaillades/dreem_herschlag/blob/main/DREEM_Herschlag/resources/library_attributes.yml>`_.


Samples.csv content
===================

See the attributes per-sample at `samples.csv attributes from dreem_herschlag <https://github.com/yvesmartindestaillades/dreem_herschlag/blob/main/DREEM_Herschlag/resources/sample_attributes.yml>`_.


RNAstructure predictions
========================

**Structure predictions**:

    * ``structure``: structure prediction of the sequence only
    * ``structure_DMS``: structure prediction of the sequence using the DMS signal
    * ``structure_ROI``: structure prediction of the ROI sub-sequence
    * ``structure_ROI_DMS``: structure prediction of the ROI sub-sequence using the DMS signal

**DeltaG predictions**:

    * ``deltaG_min``: minimum energy for the full sequence structure prediction.
    * ``deltaG_min_DMS``: minimum energy for the full sequence structure prediction using the DMS signal.
    * ``deltaG_min_ROI``: minimum energy for the region of interest sub-sequence structure prediction.
    * ``deltaG_min_ROI_DMS``: minimum energy for the region of interest sub-sequence structure prediction using the DMS signal.
    * ``deltaG_ens``: ensemble energy for the full sequence structure prediction.
    * ``deltaG_ens_ROI``: ensemble energy for the region of interest sub-sequence structure prediction.

**Base-pairing prediction**

    * ``base_pairing_prob``: probability of mutation for each residue. Computed using the full-sequence only.

.. note::

    The ensemble energy is the average of the deltaG prediction values, weighted by the alternative structures partition shares.


Poisson confidence intervals
============================

    * ``poisson_low``: lower bound of the 95% confidence interval for the real mutation rate.
    * ``poisson_high``: upper bound of the 95% confidence interval for the real mutation rate.

