.. stitchr documentation master file, first created by
   sphinx-quickstart on Tue Apr 18 16:30:31 2023.

stitchr
=======

Stitch together TCR coding nucleotide sequences from V/J/CDR3 information
-------------------------------------------------------------------------

.. image:: ../images/stitchr-logo.png
   :scale: 40 %
   :alt: the stitchr logo (which is pretty neat, if I do say so myself), in which the dot over the eye is a small pair of scissors, and the T, C, and R characters are a different shade and threaded by a needle
   :align: center

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   inputdata
   thimble
   guistitchr
   output
   immunoglobulins
   importing
   troubleshooting
   referencing

Sometimes you need a complete TCR nucleotide or amino acid sequence, but all you have is limited information. This script aims to generate a coding nucleotide sequence for a given rearrangement (e.g. for use when generating TCR expression vectors) in those situations.

The script takes the known V/J/CDR3 information, and uses that to pull out the relevant germline TCR nucleotide sequences and stitch them together. Its modular approach can be used for the automated generation of TCR sequences for gene synthesis and functional testing, or for TCR engineering through supplying modified germline sequences.

Out of the box, ``stitchr`` works on all common jawed vertebrate TCR loci (alpha/beta/gamma/delta), for all species for which there is currently data available in IMGT.

What's new in v1.3
==================

(2025-04-25)

A couple of minor case use capabilities have been added in this version.

* A ``-nl / --no_leader`` flag has been added, so that TCRs can be generated without a leader region (as `had been requested in the Issues page <https://github.com/JamieHeather/stitchr/issues/45>`_). Note that this option should be avoided for construction TCR sequences for functional expression, as these regions encode the signal peptide and are thus required for proper polypeptide trafficking.

* A ``-sn / --skip_n_checks`` flag has been added, allowing the production of `TCRs with CDR3 junctions which do not begin with the classical conserved second cysteine <https://github.com/JamieHeather/stitchr/issues/46>`_ when providing exact CDR3 junctions. Note that this option will only work for rearrangements which have deleted this CYS104 residue; in the absence of any N-terminal overlap (which previously was the Cys as a minimum), it now places the provided junction sequence with the first CDR3 residue in place of CYS104. (The only way to get ``stitchr`` to reliably produce rearrangements which have deleted beyond this position is to provide nucleotide sequence which extends further into the V region, in conjunction with the ``-sl / --seamless`` flag.)

* In order to simplify tracking of different ``stitchr`` updates, script-specific versions have been removed and updates are tracked wholly through package-level versioning. This can be easily obtained by running the ``--version`` option in any of the command scripts, or ``--cite`` for a fuller description.

What's new in v1.2
==================

(2025-02-24)

There are a number of minor tweaks and quality of life improvements in the update from version 1.1.3 to 1.2.0, with the major changes aiming to improve repeatability and user convenience, including:

* ``stitchr`` has changed how it stores and handles data, as laid out in the :ref:`output-label` section...
* ... Which has enabled outputting both ``stitchr`` and ``thimble`` results as either JSON or GenBank files!
* ``stitchrdl`` has also been updated to allow easier addition of FASTA reads to the ``additional-genes.fasta`` file, as per the :ref:`input-data-label` page.
* Finally some streamlining has increased the speed of the scripts, making ``thimble`` about a third faster.

   * See the below barplot, in which one million paired TCRs were processed with the old v1.1.3 and new v1.2.0.
   * Each version was run three times, error bars show 95% confidence intervals.

.. image:: ../images/speed-v1-2-0.png
   :width: 300
   :alt: barplot showing two groups: on the left 'v1.1.3' (the previous stitchr version) takes just over 60 seconds to process one million sequences, on the right 'v1.2.0' (the new version) takes about 40
   :align: center


Links
=====

* `stitchr on GitHub <https://github.com/JamieHeather/stitchr>`_

* `stitchr on PyPI <https://pypi.org/project/stitchr/>`_

* `stitchr publication (NAR 2022) <https://doi.org/10.1093/nar/gkac190>`_

* Related tools:
   * `autoDCR (for flexible TCR annotation) <https://github.com/JamieHeather/autoDCR>`_
   * `hladl (for automatically downloading or inferring HLA gene/protein sequences) <https://github.com/JamieHeather/hladl>`_
