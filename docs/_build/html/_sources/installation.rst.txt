Getting started
===============

Installation
------------

``stitchr`` runs on Python3, and can be installed via ``pip``:

``pip install stitchr``

In order to automatically download the necessary data for stitching, `IMGTgeneDL <https://github.com/JamieHeather/IMGTgeneDL>`_ is also required. If it's not automatically installed alongside ``stitchr``, it can be installed with:

``pip install IMGTgeneDL``

After installing ``stitchr`` via ``pip``, ``IMGTgeneDL`` can be used via the ``stitchrdl`` command to download suitably formatted data sets to the required directory like so:

``stitchrdl -s human``

See the :ref:`species-covered-label` section for details on the species for which data can be downloaded in this manner.

There are also some optional functions that require additional dependencies:

* The ``-aa`` alignment function of ``stitchr`` requires Biopython

    * ``pip install Bio``

* The graphical user interface version :ref:`gui-label` requires ``PySimpleGUI``

    * ``pip install PySimpleGUI``
    * This also might require the `installation of Tkinter <https://tkdocs.com/tutorial/install.html>`_

Quick start example
-------------------

The only required fields are the minimal components describing a single rearranged TCR chain: V gene name, J gene name, and CDR3 sequence (either DNA or amino acids). Constant regions must also be specified for all non-human/non-mouse species.

.. code:: bash

   stitchr -v [IMGT V gene] -j [IMGT J gene] -cdr3 [CDR3aa]

   stitchr -v TRBV7-3*01 -j TRBJ1-1*01 -cdr3 CASSYLQAQYTEAFF

   stitchr -v TRAV1-2 -j TRAJ33 -cdr3 TGTGCTGTGCTGGATAGCAACTATCAGTTAATCTGG

See the :ref:`usage-label` section for more detailed usage instructions. ``stitchr`` can also be run in a high-throughput manner (see :ref:`thimble-label`), or via a simple graphical user interface (see :ref:`gui-label`).