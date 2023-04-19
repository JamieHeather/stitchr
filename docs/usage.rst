.. _usage-label:

Using ``stitchr``
=================

This script can take either amino acid or nucleotide sequences for the CDR3 junction. However when submitting an amino acid CDR3 sequence, ``stitchr`` will in most cases **not** produce the actual recombined sequences that encoded the original TCRs, apart from a few edge cases (such as particularly germ-like like alpha chain rearrangements). In these cases, ``stitchr`` recreates an equivalent full length DNA sequence that will encode the same protein sequence. It aims to produce a sequence as close to germline as possible, so all CDR3 residues that *can* be germline encoded by the V and J genes are. Non-templated residues in the CDR3 (or those templated by the D gene, which is treated as non-templated for the purpose of stitching) are chosen from taking the most commonly used codon per residue.

When provided with a nucleotide CDR3 sequence ``stitchr`` can simply line up the edges of the V and the J and pop it in. (The exception that might still produce slightly different nucleotide sequences is when during non-templated deletion and addition a long stretch of V or J gene nucleotides were removed and then a different sequence coincidentally encoding the same amino acids was introduced.)

Care must be taken to ensure that the correct TCR informaton is input, ensuring that:

* **You’re using proper IMGT gene nomenclature**

    * Older/deprecated gene names will not work

* **You have the correct and full CDR3 junction sequence**

    * `I.e. running inclusively from the conserved cysteine to the conserved phenylalanine (or rarely, tryptophan) residues <http://www.imgt.org/FAQ/#question39>`_
    * Either as amino acid or DNA sequences

* **You are using the right alleles for the TCR genes in question** if known

    * I.e. the bit after the asterisk in the gene name
    * There are many known non-synonymous polymorphisms (and undoubtedly many more unknown ones) which could be impacting on antigen recognition, surface expression, and other aspects of TCR biology
    * For best results, try to get long read TCR sequence data and process it through V/J/CDR3 annotation software which is capable of a) providing allele-level resolution, and b) can take up-to-date germline reference files.

The script produces a TCR from the information given, trying to provide warnings or errors if it detects an improbable or implausible combination, yet it’s possible that the script might produce output that *looks* OK yet which does not reproduce a coding sequence for the intended TCR.

Gene/allele default behaviour
-----------------------------

If you request an allele for which there isn’t complete sequence data, the script will attempt to default to the prototypical allele (\*01) of that gene, or a preferred default allele if the ``-p`` flag is set (see below). Similarly it will attempt to use the correct leader sequences (L-PART1+L-PART2) for the specified allele, but if it can’t find one it’ll default back to the prototype’s. In both cases, if it cannot find sequence for that allele then it will attempt to use an alternative allele for which data exists in the reference (see `this comment for more details <https://github.com/JamieHeather/stitchr/issues/25#issuecomment-1146626463>`_). Note that IMGT-provided gene sequences which are ‘partial’ at either end of their sequence are discounted entirely, as full length sequences are needed. If the script is needed to stitch TCRs that make use of genes that are partial at their recombination-distal ends then you can modify the FASTA header for these entries in the installed Data directory.

For human and mouse TCRs, the script will use the TRBC gene located in the same cluster as the J gene (i.e. TRBJ1-1 through TRBJ1-6 will get TRBC1, while TRBJ2-1 through TRBJ2-7 will get TRBC2). This can be overriden (see optional arguments). Unfortunately we are not experts in TCR loci architecture of all species, so we haven’t hard-wired any other constant region assumptions, so for all other species you’ll need to explicitly state which constant region you want used.

Translation sequences
---------------------

By default ``stitchr`` does not include stop codons at the end of the coding sequence; if desired, this must be specified using the 3’ flag (``-3p``), i.e. ``-3p TAA``, ``-3p TAG``, or ``-3p TGA``. Similarly, no sequence is included before that of the IMGT-recorded L1 leader sequence. If desired, this can be added using the 5’ flag (``-5p``), e.g. to add the pre-start codon section of an optimal Kozak sequence: ``-5p GCCGCCACC``. Note that translated sequence in the output is the *whole* stitched sequence, including any added 5’/3’ sequences: addition of 5’ sequences may cause the introduction of underscores (’\_‘) to appear in the translated output, representing incomplete codons that could not be translated. Also note that the translated sequence of an individual chain may differ from the corresponding section of a linked heterodimer for this reason, depending on the length/frame of the 5’ rearrangement.

Seamless mode
-------------

If users care about accurately replicating the exact nucleotide sequence of specific V(D)J rearrangements, and they have additional nucleotide sequences beyond the edges of the CDR3 junction, they can make use of the optional ``-sl`` ‘seamless’ flag to stitch together the complete recombined sequence as faithfully as possible.

E.g. instead of these first two options:

.. code:: bash

   stitchr -v TRBV7-6 -j TRBJ1-4 -cdr3 CASSSGQGLGEKLFF
   stitchr -v TRBV7-6 -j TRBJ1-4 -cdr3 TGTGCCAGCAGTTCCGGACAGGGCTTGGGAGAAAAACTGTTTTTT

… you would run (NB non-CDR3 nucleotides shown in lower case for display purposes):

.. code:: bash

   stitchr -sl -v TRBV7-6 -j TRBJ1-4 -cdr3 catgtatcgcTGTGCCAGCAGTTCCGGACAGGGCTTGGGAGAAAAACTGTTTTTTggcagtggaa

In this example aligning the results shows that the second serine in the CDR3 was actually encoded by ‘AGT’ in the rearrangement: the ‘AGC’ codon present in the germline gene must have been deleted and this alternative ‘S’ codon added or completed by Tdt. Thus while all options should produce the same amino acid sequence, the seamless option allows for truer generation of the sequence as was present in the clonotype. Note that the seamless option adds significantly to the time it takes to run ``stitchr`` (which only really matters when running it on high-throughput datasets using ``thimble``).

In order to best use the seamless option, please ensure that:

* You have sufficient nucleotide context on either side of the CDR3 (especially the V) - ideally 20-30 nucleotides.
* Do not include any leader or constant region nucleotides - this may involve trimming nucleotide sequences.
* Ensure your V gene and allele calling is accurate, or at the very least that the contextual sequence lacks polymorphisms or errors in its very 5’.

    * ``stitchr`` will attempt to detect and deal with single nucleotide mismatches with the stated allele, but more complex polymorphisms will result in a failure.

Other optional arguments
------------------------

-  ``-h`` - see a help menu, containing all the command line options
-  ``-c`` - specify a particular constant region gene (in the case of TRBC) or allele
-  ``-s`` - specify a species: ‘HUMAN’ is the default, see :ref:`species-covered-label` section for all options (which must be downloaded with `stitchrdl` or manually produced prior to use)
-  ``-aa`` - provide an incomplete amino acid sequence (spanning at least the CDR3, with some padding on either side), to assess the accuracy of the stitched TCR sequence. Must be a single string,unbroken by spaces or linebreaks
-  ``-cu`` - specify the path to an alternative codon usage file, from which to generate the sequences for the non-templated residues (see the :ref:`codon-files-label` section)
-  ``-p`` - specify a path containing gene allele preferences (see the :ref:`preferred-allele-label` section)
-  ``-l`` - use a different leader region to that present with the given V
-  ``-n`` - provide a name for the TCR chain, which will be included in the FASTA file header
-  ``-3p`` - provide a sequence to come immediately after the end of the constant region (e.g. a stop codon)
-  ``-5p`` - provide a sequence to come immediately before the start of the L1 leader sequence (e.g. a Kozak sequence)
-  ``-m`` - define an output mode, to define which sequences get printed to the terminal
-  ``-xg`` - toggle providing additional/custom genes to be stitched into TCR transcripts in the Data/additional-genes.fasta file
-  ``-sc`` - toggle skipping the constant region gene check (for genes not present in the C-region-motifs.tsv file)
-  ``-sw`` - suppress warning text, which may be especially useful in conjunction with some of the alternative output modes (see the :ref:`output-modes-label` section)

.. _output-modes-label:

Output modes
~~~~~~~~~~~~

``stitchr`` can output the TCR sequences it generates in a number of different formats, which `may help integrate its output into certain pipelines <https://github.com/JamieHeather/stitchr/issues/22>`_. These modes can be specified using the ``-m / --mode`` flag, using one of the following options:

-  ``-m BOTH_FA``

   -  Default option
   -  Outputs a horizontal line, followed by the full, formatted, descriptive FASTA sequence of the stitched TCR, both nucleotide and translated amino acid sequence

-  ``-m NT_FA``

   -  Outputs a horizontal line and the FASTA nucleotide sequence of the stitched TCR

-  ``-m AA_FA``

   -  Outputs a horizontal line and the FASTA translated amino acid sequence of the stitched TCR

-  ``-m NT``

   -  Outputs just the nucleotide sequence of the stitched TCR (no lines, no linebreaks, no FASTA header)

-  ``-m AA``

   -  Outputs just the translated amino acid sequence of the stitched TCR (no lines, no linebreaks, no FASTA header)

Providing a partial amino acid sequence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you provide a partial amino acid sequence using the ``-aa`` flag, ``stitchr`` will perform a rudimentary pairwise alignment, just to give a quick visual assessment of the quality of the sequence generation.

As an example, let’s take the example of the well described A2-NLV restricted `C25 TCR from the 5D2N PDB structure <https://www.rcsb.org/structure/5d2n>`_. We can take the amino acid sequence straight from the PDB FASTA file:

::

   >5D2N:E|PDBID|CHAIN|SEQUENCE
   MGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTI
   QRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVEL
   SWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSA
   EAWGRAD

We can then pull out the V, J, and CDR3 information. There’s lots of ways to do this, but the easiest manual way is to find the CDR3 and then search the immediately neighbouring sequences against V/J amino acid sequences (obtainable via IMGT/GENE-DB). This gives:

``TRBV7-6 / TRBJ1-4 / CASSLAPGTTNEKLFF``

Then we can run the code like this:

.. code:: bash

   stitchr -v TRBV7-6 -j TRBJ1-4 -cdr3 CASSLAPGTTNEKLFF -n C25 -aa MGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD

   # Produces the following output

   >nt-C25-TRBV7-6*01-TRBJ1-4*01-TRBC1*01-CASSLAPGTTNEKLFF-leader-TRBV7-6*01
   ATGGGCACCAGTCTCCTATGCTGGGTGGTCCTGGGTTTCCTAGGGACAGATCACACAGGTGCTGGAGTCTCCCAGTCTCCCAGGTACAAAGTCACAAAGAGGGGACAGGATGTAGCTCTCAGGTGTGATCCAATTTCGGGTCATGTATCCCTTTATTGGTACCGACAGGCCCTGGGGCAGGGCCCAGAGTTTCTGACTTACTTCAATTATGAAGCCCAACAAGACAAATCAGGGCTGCCCAATGATCGGTTCTCTGCAGAGAGGCCTGAGGGATCCATCTCCACTCTGACGATCCAGCGCACAGAGCAGCGGGACTCGGCCATGTATCGCTGTGCCAGCAGCCTGGCCCCCGGCACCACTAATGAAAAACTGTTTTTTGGCAGTGGAACCCAGCTCTCTGTCTTGGAGGACCTGAACAAGGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTTCCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACGGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCCGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTTACCTCGGTGTCCTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCCTGCTAGGGAAGGCCACCCTGTATGCTGTGCTGGTCAGCGCCCTTGTGTTGATGGCCATGGTCAAGAGAAAGGATTTC

   >aa-C25-TRBV7-6*01-TRBJ1-4*01-TRBC1*01-CASSLAPGTTNEKLFF-leader-TRBV7-6*01
   MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDF

   MG------------------AGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQ
   ||                  ||||||||||||||||||||||||||||||||||||||||
   MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQ

   GPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTT
   ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   GPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTT

   NEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFY-PDHVELS
   |||||||||||||||||||||||||||||||||||||||||||||||||||  |||||||
   NEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGF-FPDHVELS

   WWVNGKEVHSGVC-TDPQPLKEQPALNDSRYA-LSSRLRVSATFWQNPRNHFRCQVQFYG
   ||||||||||||  |||||||||||||||||  |||||||||||||||||||||||||||
   WWVNGKEVHSGV-STDPQPLKEQPALNDSRY-CLSSRLRVSATFWQNPRNHFRCQVQFYG

   LSENDEWTQDRAKPVTQIVSAEAWGRAD--------------------------------
   ||||||||||||||||||||||||||||
   LSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVL

   ---------------

   VSALVLMAMVKRKDF

We can see that there’s a few mismatches in the latter half of the stitched sequence, so perhaps this crystal actually used the other TRBC gene. We can swap that in:

.. code:: bash

   stitchr -v TRBV7-6 -j TRBJ1-4 -cdr3 CASSLAPGTTNEKLFF -n C25 -c TRBC2 -aa MGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD

   # Produces:

   >nt-C25-TRBV7-6*01-TRBJ1-4*01-TRBC2*01-CASSLAPGTTNEKLFF-leader-TRBV7-6*01
   ATGGGCACCAGTCTCCTATGCTGGGTGGTCCTGGGTTTCCTAGGGACAGATCACACAGGTGCTGGAGTCTCCCAGTCTCCCAGGTACAAAGTCACAAAGAGGGGACAGGATGTAGCTCTCAGGTGTGATCCAATTTCGGGTCATGTATCCCTTTATTGGTACCGACAGGCCCTGGGGCAGGGCCCAGAGTTTCTGACTTACTTCAATTATGAAGCCCAACAAGACAAATCAGGGCTGCCCAATGATCGGTTCTCTGCAGAGAGGCCTGAGGGATCCATCTCCACTCTGACGATCCAGCGCACAGAGCAGCGGGACTCGGCCATGTATCGCTGTGCCAGCAGCCTGGCCCCCGGCACCACTAATGAAAAACTGTTTTTTGGCAGTGGAACCCAGCTCTCTGTCTTGGAGGACCTGAAAAACGTGTTCCCACCCGAGGTCGCTGTGTTTGAGCCATCAGAAGCAGAGATCTCCCACACCCAAAAGGCCACACTGGTGTGCCTGGCCACAGGCTTCTACCCCGACCACGTGGAGCTGAGCTGGTGGGTGAATGGGAAGGAGGTGCACAGTGGGGTCAGCACAGACCCGCAGCCCCTCAAGGAGCAGCCCGCCCTCAATGACTCCAGATACTGCCTGAGCAGCCGCCTGAGGGTCTCGGCCACCTTCTGGCAGAACCCCCGCAACCACTTCCGCTGTCAAGTCCAGTTCTACGGGCTCTCGGAGAATGACGAGTGGACCCAGGATAGGGCCAAACCTGTCACCCAGATCGTCAGCGCCGAGGCCTGGGGTAGAGCAGACTGTGGCTTCACCTCCGAGTCTTACCAGCAAGGGGTCCTGTCTGCCACCATCCTCTATGAGATCTTGCTAGGGAAGGCCACCTTGTATGCCGTGCTGGTCAGTGCCCTCGTGCTGATGGCCATGGTCAAGAGAAAGGATTCCAGAGGC

   >aa-C25-TRBV7-6*01-TRBJ1-4*01-TRBC2*01-CASSLAPGTTNEKLFF-leader-TRBV7-6*01
   MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQGPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTTNEKLFFGSGTQLSVLEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDSRG

   MG------------------AGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQ
   ||                  ||||||||||||||||||||||||||||||||||||||||
   MGTSLLCWVVLGFLGTDHTGAGVSQSPRYKVTKRGQDVALRCDPISGHVSLYWYRQALGQ

   GPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTT
   ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   GPEFLTYFNYEAQQDKSGLPNDRFSAERPEGSISTLTIQRTEQRDSAMYRCASSLAPGTT

   NEKLFFGSGTQLSVLEDLNK-VFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELS
   |||||||||||||||||| | |||||||||||||||||||||||||||||||||||||||
   NEKLFFGSGTQLSVLEDL-KNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELS

   WWVNGKEVHSGVC-TDPQPLKEQPALNDSRYA-LSSRLRVSATFWQNPRNHFRCQVQFYG
   ||||||||||||  |||||||||||||||||  |||||||||||||||||||||||||||
   WWVNGKEVHSGV-STDPQPLKEQPALNDSRY-CLSSRLRVSATFWQNPRNHFRCQVQFYG

   LSENDEWTQDRAKPVTQIVSAEAWGRAD--------------------------------
   ||||||||||||||||||||||||||||
   LSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSESYQQGVLSATILYEILLGKATLYAVL

   -----------------

   VSALVLMAMVKRKDSRG

This produces even more mismatches! This is an instance where the constant region used in the crystal has been altered for expression/crystallization purposes.

A note on CDR3 C-terminal residues
----------------------------------

``stitchr`` assumes that the J gene will not undergo deletion past the C-terminal residue of the CDR3 junction (which occurs approximately in the middle of the J). Thus the code looks for the appropriate residue at the end of the CDR3, which in the majority of cases will be a phenylalanine (F). However in some cases it might be something else, like a W (not uncommon in human TRAJ/mice genes) or even something more exotic like a C, L or H (which occur in certain mouse J genes). Note that most of these non-F/W residues are found in J genes with a predicted `‘ORF’ IMGT status <http://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html>`_, and thus might not contribute to functioning TCRs, but ``stitchr`` will still let you generate a plausible sequence using them.
