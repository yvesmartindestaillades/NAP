=========
Data reference
=========


DREEM's columns
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



