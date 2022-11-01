
.. _terminology:

===========
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

:Cluster:

    A **cluster** is one of the outputs of DREEM's Expectation-Maximization clustering algorithm applied to the mutation rates of a construct in a sample is clustered.
    Clusters are numbered 0, 1, 2 and so on.


:Sample-construct-cluster (SCC):

    A **SCC** corresponds to a cluster of a construct in a sample.
    Each SCC has a set of **attributes**.


:Attributes:

    **Attributes** are a set of experimental results and properties related to a **sample-construct**, read from dreem package's mutation_fraction.  
    
    Examples: ``flank``, ``mut_bases``, ``user``. 

:Study:

    A **study** is a set of **samples** that are relevant to be studied together, typically because of their experimental conditions.
    For example, a study can be a set of replicates, or the same library with temperature variating, etc.

    The study object contains methods to load and plot the data for its samples.

:ROI:

    **ROI** is the acronym for Region of Interest.
    It corresponds to a sub-sequence in the RNA sequence that we want to study.

:DREEM pickle file:

    **DREEM pickle file** is the output of DREEM pipeline.
    It consists in a class call ``MutationHistogram`` that contains one **sample**.
    
    ``MutationHistogram`` is compressed using the Python module ``pickle``.

