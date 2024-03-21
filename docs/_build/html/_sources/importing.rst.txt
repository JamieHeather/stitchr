.. _importing-label:

Importing ``stitchr`` for use in other Python scripts
===================================================

The underlying core functionality of ``stitchr`` can be imported for use in other Python scripts, as is done to power the ``thimble`` and ``gui_stitchr`` tools.

The main things to remember are that:

* Whatever script is importing ``stitchr`` will need to initialise the necessary data using the ``get_imgt_data`` function in the ``stitchrfunctions`` script

    * This must be done *per chain* that needs to be stitched
    * Therefore if users wish to stitch heterodimers or multiple loci in one session, they must generate the necessary data as separate variables or dictionary entries
    * The underlying data to be stitched must be installed in the ``stitchr`` data directory (see the :ref:`input-data-label` section)

* The core ``stitch`` function works on **single rearrangements**

    * Heterodimers are therefore produced by calling ``stitch`` twice, once per locus
    * Each call of the ``stitch`` function requires certain values

Here's a simple example of how to import and call ``stitchr``:

.. code:: python

    # import stitchr
    from Stitchr import stitchrfunctions as fxn
    from Stitchr import stitchr as st

    # specify details about the locus to be stitched
    chain = 'TRB'
    species = 'HUMAN'

    # initialise the necessary data
    tcr_dat, functionality, partial = fxn.get_imgt_data(chain, st.gene_types, species)
    codons = fxn.get_optimal_codons('', species)

    # provide details of the rearrangement to be stitched
    tcr_bits = {'v': 'TRBV7-3*01', 'j': 'TRBJ1-1*01', 'cdr3': 'CASSYLQAQYTEAFF',
                'l': 'TRBV7-3*01', 'c': 'TRBC1*01',
                'skip_c_checks': False, 'species': species, 'seamless': False,
                '5_prime_seq': '', '3_prime_seq': '', 'name': 'TCR'}

    # then run stitchr on that rearrangement
    stitched = st.stitch(tcr_bits, tcr_dat, functionality, partial, codons, 3, '')

    print(stitched)
    # Which produces
    (['TCR', 'TRBV7-3*01', 'TRBJ1-1*01', 'TRBC1*01', 'CASSYLQAQYTEAFF', 'TRBV7-3*01(L)'],
     'ATGGGCACCAGGCTCCTCTGCTGGGCAGCCCTGTGCCTCCTGGGGGCAGATCACACAGGTGCTGGAGTCTCCCAGACCCCCAGTAACAAGGTCACAGAGAAGGGAAAATATGTAGAGCTCAGGTGTGATCCAATTTCAGGTCATACTGCCCTTTACTGGTACCGACAAAGCCTGGGGCAGGGCCCAGAGTTTCTAATTTACTTCCAAGGCACGGGTGCGGCAGATGACTCAGGGCTGCCCAACGATCGGTTCTTTGCAGTCAGGCCTGAGGGATCCGTCTCTACTCTGAAGATCCAGCGCACAGAGCGGGGGGACTCAGCCGTGTATCTCTGTGCCAGCAGCTACCTGCAGGCCCAGTACACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAGAGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTTCCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACGGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTTACCTCGGTGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATGGTCAAGAGAAAGGATTTC',
     0)

Note that the ``stitch`` function outputs three items:

* A list of the relevant TCR features used to stitch that rearrangement, which default ``stitchr`` uses to compile a FASTA header
* A string detailing the nucleotide sequence of the stitched TCR
* An integer (0-2) detailing the translation offset (only relevant if a user-supplied 5 prime sequence has been included)

If users elect to do import ``stitchr`` into their own pipelines we recommend that they familiarise themselves with the code, particularly how ``stitchr`` is called in ``thimble`` and ``gui_stitchr``. The should also validate their results using vanilla ``stitchr`` and/or known TCR sequence controls. It is also recommended that users pay attention to the warnings produced, which can be instructive even for properly stitched sequences.





