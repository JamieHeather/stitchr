.. _thimble-label:

Thimble
=======

Run ``stitchr`` high-throughput on multiple and paired TCRs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Instead of running ``stitchr`` on rearrangements one by one, you can fill out the necessary details into a tab separated file (.tsv) and submit it to ``thimble``. The format of the input data can be found in the empty and example templates located in the ``templates/`` directory.

Note that there are two kinds templates, one for alpha/beta TCRs, and another for gamma/delta TCRs, with the only difference being the gene names in the header fields. Users can only use ``thimble`` to stitch TCRs of one type per operation, and thus cannot mix a/b and g/d TCRs in the same input files.

You can tell ``thimble`` what flavour of TCR you’re making directly, using the ``-r / --receptor`` field followed by a single or double digit string (case-insensitive), e.g.:

.. code:: bash

   # Alpha/beta TCRs
   thimble -in somefile.tsv -r a
   thimble -in somefile.tsv -r AB
   thimble -in somefile.tsv -r b
   # Gamma/delta TCRs
   thimble -in somefile.tsv -r g
   thimble -in somefile.tsv -r GD
   thimble -in somefile.tsv -r dg

Alternatively if you don’t use the ``-r`` flag, ``thimble`` will automatically infer the TCR loci from the header line of the input file. While the ‘TRA-TRB’ and ‘TRG-TRD’ labels are not explicitly used by
``thimble``, they are used by the ``gui-stitchr`` script described below, and help make it clearer what’s in which files.

The species can be explicitly set via the ``-s / --species`` flag, or by including the common name in the input file somewhere. Note that using either the receptor or species flag explicitly will take precedence over any details inferred from input file.

All of the recombination-specific fields that can ordinarily be specified at the command line in ``stitchr`` can be applied per row using ``thimble``, with the exception of species (which must be kept the same for all TCRs in a given ``thimble`` run).

Note that the input to ``thimble`` can also be used to generate rearrangements for both chains of a given TCR (a/b or g/d) on one row, with additional options to link those sequences together (e.g. for gene synthesis). A number of x2A potential linkers are provided in the Data/linkers.tsv file. If custom linkers are desired, you can either edit that linkers file or just enter the nucleotide sequence of the desired linker into the Linker column of the input tsv. ``thimble`` will allow linkers that disrupt the frame (i.e. have a length not divisible by 3) but will throw a warning, so use carefully. 5’ and 3’ sequences can be added to both ends of either chain in a heterodimer, again allowing but throwing a warning if a custom sequence could potentiallydisrupt the frame.

By default, ``thimble`` produces linked TCRs in the order 5’ - beta chain - linker - alpha chain - 3’, as `this has been shown to increase the surface presentation of ectopic TCRs <https://doi.org/10.1038/mtna.2012.52>`_. However this can still be specified with the ‘Link_order’ column in the input template file, using ‘AB’ or ‘BA’ to refer to ‘alpha-beta’ or ‘beta-alpha’ orders respectively. Link order is ignored if no linker is provided. The same holds true for gamma-deltas (defaulting to DG over GD).

Any warnings and errors generated on a per-TCR basis are recorded in the final output file; it is recommended that users check this information, to ensure they understand the potential limitations of a specific sequence.

.. _example-usage-1:

Example usage
-------------

.. code:: bash

   thimble -in [input tsv] -o [output tsv]

   thimble -in thimble_input_example_TRA-TRB.tsv -o testing

The test input example is available in the ``templates/`` directory of the Github repo.

Stitch multiple TCRs per line with list and wildcard fields
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some applications may require stitching multiple variants of a particular rearrangement. When fields are specified with multiple options, a new line will be generated for each specified option, or one for each combination of options if more than one field with multiple fields are specified. Thimble allows you to specify multiple options for a particular field in several ways:

-  For almost every field, multiple options can be provided by separating options with commas.

   -  E.g. ``TRAV1-1\*01,TRAV1-2\*01,TRAV1-2\*02`` for TRAV or``CAVLDSNYQLIW,CAVLMSNYQLIW`` for TRA_CDR3.

-  For every field that relates to a specific germline TCR region (L/V/J/C), every allele for a given gene can be specified by using the ``%`` wildcard in place of an allele number.

   -  E.g. using default IMGT-GENE/DB, ``TRBV19*%`` will stitch the rearrangement using each of ``TRBV19\*01``, ``TRBV19\*02``, and ``TRBV19\*03``

-  For the same fields, users can opt to systemically try every available gene and allele for a region by *only* entering ``%`` in that field.

   -  E.g. placing ``%`` in the TRBJ field will stitch that rearrangement in combination with every different beta J gene.
   -  Note that this option is expected to generate a large number of warning messages, as not all genes are likely or able to take part in certain rearrangements.

Optional arguments
~~~~~~~~~~~~~~~~~~

-  ``-h`` - see a help menu, containing all the command line options
-  ``-s`` - specify a species, as with ``stitchr``
-  ``-cu`` - use an alternative codon usage file, from which to generate the sequences for the non-templated residues (see below)
-  ``-p`` - specify a path containing gene allele preferences (see below)
-  ``-r`` - specify the kind of TCR, i.e. a/b or g/d
-  ``-xg`` - toggle providing additional/custom genes to be stitched into TCR transcripts in the Data/additiona-genes.fasta file
-  ``-jt`` - length of J gene substring that has to be matched to avoid throwing a warning (decrease to get fewer notices about short J matches), default = 3
