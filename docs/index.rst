.. stitchr documentation master file, created by
   sphinx-quickstart on Tue Apr 18 16:30:31 2023.

stitchr
=======

Stitch together TCR coding nucleotide sequences from V/J/CDR3 information
-------------------------------------------------------------------------

.. image:: ../images/stitchr-logo.png
   :scale: 40 %
   :alt: the stitchr logo (which is pretty neat, if I do say so myself)
   :align: center

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   installation
   usage
   inputdata
   thimble
   guistitchr
   immunoglobulins
   referencing

Sometimes you need a complete TCR nucleotide or amino acid sequence, but all you have is limited information. This script aims to generate a coding nucleotide sequence for a given rearrangement (e.g. for use when generating TCR expression vectors) in those situations.

The script takes the known V/J/CDR3 information, and uses that to pull out the relevant germline TCR nucleotide sequences and stitch them together. Its modular approach can be used for the automated generation of TCR sequences for gene synthesis and functional testing, or for TCR engineering through supplying modified germline sequences.

Out of the box, ``stitchr`` works on all common jawed vertebrate TCR loci (alpha/beta/gamma/delta), for all species for which there is currently data available in IMGT.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
