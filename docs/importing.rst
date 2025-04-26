.. _importing-label:

Importing ``stitchr`` for use in other Python scripts
=====================================================

The underlying core functionality of ``stitchr`` can be imported for use in other Python scripts, as is done to power the ``thimble`` and ``gui_stitchr`` tools. Note that **as of ``stitchr`` version 1.2.0, the output format of the ``stitch()`` function has changed, so pipelines calling it will need to be updated!**

The main things to remember are that:

* Whatever script is importing ``stitchr`` will need to initialise the necessary data using the ``get_ref_data`` function in the ``stitchrfunctions`` script

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
    tcr_dat, functionality, partial = fxn.get_ref_data(chain, st.gene_types, species)
    codons = fxn.get_optimal_codons('', species)
    j_res, low_conf_js = fxn.get_j_motifs(species)
    c_res = fxn.get_c_motifs(species)

    # provide details of the rearrangement to be stitched
    tcr_bits = {'v': 'TRBV7-3*01', 'j': 'TRBJ1-1*01', 'cdr3': 'CASSYLQAQYTEAFF',
                'l': 'TRBV7-3*01', 'c': 'TRBC1*01', 'mode': '',
                'skip_c_checks': False, 'skip_n_checks': False, 'no_leader': False,
                'species': species, 'seamless': False,
                '5_prime_seq': '', '3_prime_seq': '', 'name': 'my-cool-TCR'}

    # then run stitchr on that rearrangement
    stitched = st.stitch(tcr_bits, tcr_dat, functionality, partial, codons, 3, '',
                         c_res, j_res, low_conf_js)

    print(stitched)
    {'in': {'v': 'TRBV7-3*01', 'j': 'TRBJ1-1*01', 'cdr3': 'CASSYLQAQYTEAFF', 'l': 'TRBV7-3*01', 'c': 'TRBC1*01', 'mode': '', 'skip_c_checks': False, 'skip_n_checks': False, 'no_leader': False, 'species': 'HUMAN', 'seamless': False, '5_prime_seq': '', '3_prime_seq': '', 'name': 'my-cool-TCR'}, 'used': {'l': 'TRBV7-3*01(L)', 'v': 'TRBV7-3*01', 'j': 'TRBJ1-1*01', 'c': 'TRBC1*01', 'cdr3': 'CASSYLQAQYTEAFF', 'non_templated': 'non_templated'}, 'seqs': {'l': 'ATGGGCACCAGGCTCCTCTGCTGGGCAGCCCTGTGCCTCCTGGGGGCAGATCACACA', 'v': 'GGTGCTGGAGTCTCCCAGACCCCCAGTAACAAGGTCACAGAGAAGGGAAAATATGTAGAGCTCAGGTGTGATCCAATTTCAGGTCATACTGCCCTTTACTGGTACCGACAAAGCCTGGGGCAGGGCCCAGAGTTTCTAATTTACTTCCAAGGCACGGGTGCGGCAGATGACTCAGGGCTGCCCAACGATCGGTTCTTTGCAGTCAGGCCTGAGGGATCCGTCTCTACTCTGAAGATCCAGCGCACAGAGCGGGGGGACTCAGCCGTGTATCTCTGTGCCAGCAGC', 'j': 'ACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAG', 'c': 'AGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTTCCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACGGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTTACCTCGGTGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATGGTCAAGAGAAAGGATTTC', 'non_templated': 'TACCTGCAGGCCCAGTAC'}, 'input_type': 'aa', 'translation_offset': 0, 'out_list': ['my-cool-TCR', 'TRBV7-3*01', 'TRBJ1-1*01', 'TRBC1*01', 'CASSYLQAQYTEAFF', 'TRBV7-3*01(L)'], 'stitched_nt': 'ATGGGCACCAGGCTCCTCTGCTGGGCAGCCCTGTGCCTCCTGGGGGCAGATCACACAGGTGCTGGAGTCTCCCAGACCCCCAGTAACAAGGTCACAGAGAAGGGAAAATATGTAGAGCTCAGGTGTGATCCAATTTCAGGTCATACTGCCCTTTACTGGTACCGACAAAGCCTGGGGCAGGGCCCAGAGTTTCTAATTTACTTCCAAGGCACGGGTGCGGCAGATGACTCAGGGCTGCCCAACGATCGGTTCTTTGCAGTCAGGCCTGAGGGATCCGTCTCTACTCTGAAGATCCAGCGCACAGAGCGGGGGGACTCAGCCGTGTATCTCTGTGCCAGCAGCTACCTGCAGGCCCAGTACACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAGAGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTTCCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACGGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTTACCTCGGTGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATGGTCAAGAGAAAGGATTTC'}


Note that the ``stitch`` function returns a nested dictionary, containing various levels that relate to information about the rearrangement in question at different stages of its handling. These fields are described in the :ref:`output-label` page (although note that the raw ``stitch()`` function output is only a subset of the fields described there, as the individual scripts each add additional data post-stitching.

If users elect to do import ``stitchr`` into their own pipelines we recommend that they familiarise themselves with the code, particularly how ``stitchr`` is called in ``thimble`` and ``gui_stitchr``. The should also validate their results using vanilla ``stitchr`` and/or known TCR sequence controls. It is also recommended that users pay attention to the warnings produced (which most ``stitchr`` scripts keep track of using the ``warnings.catch_warnings()`` function), which can be instructive even for properly stitched sequences.
