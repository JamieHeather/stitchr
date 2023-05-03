.. _input-data-label:

``stitchr`` input data
======================

.. _species-covered-label:

Species covered
~~~~~~~~~~~~~~~

``stitchr`` can stitch any TCR loci for which it has the necessary raw data, appropriately formatted in the ``Data`` directory (which can be located by running ``stitchr -dd`` or ``stitchr --data_dir``). IMGT currently contains enough TCR data for stitching (i.e. at least one leader, V gene, J gene, and constant region for a given loci) for the following species/loci:

.. csv-table:: stitchr species/loci table
    :header: "Common Name", "Genus species", "TRA", "TRB", "TRG", "TRD", "Kazusa species ID"
    :widths: auto

    "**CAT**", "*Felis catus*", "✔", "✔", "✔", "✔", "`9685 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9685>`_"
    "**COW**", "*Bos taurus*", "✔", "✔", "✔", "✔", "`9913 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9913>`_"
    "**CYNOMOLGUS\_MONKEY**", "*Macaca fascicularis*", , "✔", , , "`9541 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9541>`_"
    "**DOG**", "*Canis lupus familiaris*", "✔", "✔", "✔", "✔", "`9615 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9615>`_"
    "**DOLPHIN**", "*Tursiops truncatus*", "✔", , "✔", "✔", "`9739 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9739>`_"
    "**DROMEDARY**", "*Camelus dromedarius*", , "✔", "✔", , "`9838 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9838>`_"
    "**FERRET**", "*Mustela putorius furo*", , "✔", , , "`9669 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9669>`_"
    "**HUMAN**", "*Homo sapiens*", "✔", "✔", "✔", "✔", "`9606 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606>`_"
    "**MOUSE**", "*Mus musculus*", "✔", "✔", "✔", "✔", "`10090 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=10090>`_"
    "**NAKED\_MOLE-RAT**", "*Heterocephalus glaber*", "✔", "✔", "✔", "✔", "-"
    "**PIG**", "*Sus scrofa*", , "✔", , , "`9823 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9823>`_"
    "**RABBIT**", "*Oryctolagus cuniculus*", , "✔", "✔", "✔", "`9986 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9986>`_"
    "**RHESUS\_MONKEY**", "*Macaca mulatta*", "✔", "✔", "✔", "✔", "`9544 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9544>`_"
    "**SHEEP**", "*Ovis aries*", "✔", "✔", , "✔", "`9940 <https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9940>`_"


The species can be specified using the ``-s / --species`` command line flag when running your ``stitchr`` command (default = 'human'). E.g., here's an example using everyone's favourite mouse TCR, OT-I (sequences inferred from `this plasmid on AddGene <https://www.addgene.org/52111/>`_:

.. code:: bash

   stitchr -s mouse -v TRBV12-1 -j TRBJ2-7 -cdr3 CASSRANYEQYF
   stitchr -s mouse -v TRAV14D-1 -j TRAJ33 -cdr3 CAASDNYQLIW

It must be noted that many of the less well studied species have poorer gene annotations and germline variation covered by IMGT, so TCRs produced using these datasets should be treated with more caution than say for humans. E.g. different mouse strains will have different alleles (and different numbers of gene family members), so accuracy of stitched TCRs will depend both on the quality of both germline gene information and TCR clonotyping.

Generating new IMGT input files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may wish to update your raw TCR data files after an IMGT update, as sequences can be added (or even changed), or as more species become available. Assuming IMGT maintains its current noming conventions and webhosting scheme, these tasks can be undertaken automatically using the ``stitchrdl`` command, which is a wrapper for the `IMGTgeneDL script <https://github.com/JamieHeather/IMGTgeneDL>`_ (details in its own repo).

Users should note that when this is run, it creates a ``data-production-date.tsv`` file in the directory of that species, which contains the IMGT release number and date of downloaded sequences,
which should be included in any published reporting of the TCR sequences used. We also recommend that you update the germline TCR data for ``stitchr`` at the same time you update the database used in whatever TCR gene annotation software you use, to ensure that there's no discrepancy in allele nomenclature between the two.

.. _data-formatting-label:

Stitchr data formatting
~~~~~~~~~~~~~~~~~~~~~~~

Each species you wish to stitch TCRs for must have its own folder in the installed ``Data/`` directory, named after whatever flag you wish you use when giving ``stitchr`` information through the ``-s / --species`` flag. Note that ``stitchr`` will assume every folder in ``Data/`` that isn't named 'kazusa' or 'GUI-Examples' is a potential TCR germline folder, so it's advised to not put any other folders there. It is recommended to use ``stitchrdl`` to generate and update these folders
where possible.

Inside that folder there should be various files:

* ``data-production-date.tsv``

    * Contains information about the IMGT and script versions used to generate this data
    * Technically not used by ``stitchr`` as it runs, but contains important information for recording or relaying the products of ``stitchr``

* ``imgt-data.fasta``

    * Contains all of the FASTA reads that were successfully downloaded for this species
    * Again, this file isn't used during stitching but it's a useful reference to have

* ``J-region-motifs.tsv``

    * Contains automatically inferred CDR3 junction ending motifs and residues (using the process established in `the autoDCR TCR assignation tool <https://github.com/JamieHeather/autoDCR>`_), for use in finding the ends of junctions in ``stitchr``

* ``C-region-motifs.tsv``

    * Contains automatically inferred in-frame constant region peptide sequences, for use in finding the correct frame of stitched sequences

* ``TR[A/B/G/D].fasta``

    * FASTA files of the individual loci's genes
    * FASTA headers require a flag specifying what type of gene they are at the very end, after a tilde (~) character, being one of LEADER/VARIABLE/JOINING/CONSTANT

Note that the method used to automatically download TCR data (`IMGTgeneDL <https://github.com/JamieHeather/IMGTgeneDL>`_) seems to struggle for constant regions in certain species, when trying to download the fully spliced sequences. As such, the script will instead download each of the individual exons and splice them together.

This is particularly important for users who wish to stitch gamma chain TCRs with constant regions which may have multiple possible exon configurations. Any constant region which uses non-conserved or duplicated exons has an additional suffix to the allele field, in which the non-standard exon labels are appended after an underscore. E.g. the human ``TRGC2*05`` gene has the arrangement EX1+EX2\ **T**\ +EX2\ **R**\ +EX2+EX3 (instead of the usual EX1+EX2+EX3), so it is labelled ``TRGC2*05_TR``. This allows users to have multiple isoforms of the same allele with different sequences.

.. _codon-files-label:

Codon usage files
^^^^^^^^^^^^^^^^^

Non-templated based are assigned by taking the most common nucleotide triplet for a given amino acid, in a provided codon usage file.

Codon usage files are provided for all species for which data is available on the the `Kazusa website <https://www.kazusa.or.jp/>`_, and can be found in the ``Data/kazusa/`` directory. Alternative files can be provided, but must be in the same format (e.g. those provided by `HIVE <https://hive.biochemistry.gwu.edu/dna.cgi?cmd=refseq_processor&id=569942>`_), and named according to the common species name used for the rest of the data and placed in that directory if not specified using the ``-cu / --codon_usage`` flag. U/T can be used interchangeably, as all U bases will be replaced with T anyway.

If no species-specific codon usage file is found the script will default to using the human file.

.. _preferred-allele-label:

Preferred allele files
^^^^^^^^^^^^^^^^^^^^^^

``stitchr`` requires an exact allele to know which sequence to pull out of the database. By default, it always prioritises using the exact allele specified, but if the user just gives a gene identifier without an allele specified (e.g. TRAV1-1) then ``stitchr`` will use the prototypical '01' allele (e.g. TRAV1-1\*01), as this is the only allele which every IMGT-provided gene should theoretically have. It will similarly default to \*01 if explicitly given an allele which it can't find in the
input data.

There are occasions when this is not the biologically appropriate allele to choose. While users can of course specify the allele explicitly when providing the gene, they may alternatively wish to make use of the ``-p / --preferred_alleles_path`` command line option, which allows them to point to a tab-delimited file detailing specific alleles which should be used. Here users can specify four fields:

-  **Gene**: the IMGT gene name
-  **Allele**: the preferred allele (i.e. the text after the '\*' in a complete name)
-  **Region**: one of LEADER/VARIABLE/JOINING/CONSTANT, which tells ``stitchr`` explicitly what kind of sequence it is
-  **Loci**: specify which locus or loci this preferred allele covers using three digit codes, comma-delimiting if >1 (e.g. “TRB” or “TRA,TRD”)
-  **Source**: this field is not used by ``stitchr``, but can be useful for keeping track of the origin of/reason for including each preferred allele

This feature is particularly of use when generating large numbers of stitched sequences from a particular individual or strain where non-prototypical alleles are known. Note that if you are specifying adifferent allele for a variable gene that has variants in its leader sequence as well, make sure you add entries for both VARIABLE and LEADER alleles.

A template and two example common mouse strain files are included in the ``templates/preferred-alleles/`` directory. These examples are for the common mouse strains C56/Bl6 and Balb/c, and were produced by using the subspecies/strain field of the IMGT headers. Note that even then users should take care, as some genes have multiple alleles associated with them, despite being from inbred mice – e.g. TRAV9D-2 has two alleles associated with it for Balb/c (01 and 03). I've tried to pick the ones that are more likely to be functional (F > ORF > P, e.g. choosing ``TRBV24*03`` over 02 for Balb/c, or ``TRAV9D-4*04`` over 02 for C57/BL6), or are from better inferred data (e.g. taking one with functionality 'F' over '(F)').

Providing additional gene sequences
-----------------------------------

Sometimes you may wish to generate TCRs using additional gene sequences which won't be provided by IMGT (at least in the context of a given species). This can be used to introduce sequences from other loci/species, and modified or otherwise non-naturally occurring gene combinations.

Genes to be included can be added to the Data/additional-genes.fasta file, and then when ``stitchr`` or ``thimble`` is run these sequences will be read in by use of the ``-xg`` flag. As constant region gene switching is a common modification used in TCR expression and engineering studies, human alpha/beta/gamma/delta and mouse alpha/beta constant regions have been preloaded into this file. Genes added to thisfasta must have a FASTA header in the format:

::

   >accession/ID|gene*allele|species|functionality|sequence type
   e.g.
   >X02883|hTRAC*01|Homo sapiens|F|EX1+EX2+EX3+EX|anything else...

Only the second and fifth fields are important for these additional genes, and all other fields can be left empty (with empty functionality calls being presumed functional). The second field contains gene name and allele information: the gene name can be any alphanumeric string (that doesn't contain an asterisk), while the allele should be a zero-padded two (or more) digit integer (e.g. '01'). Any case can be used in gene names, but bear in mind all will be made upper case when running. The fifth field corresponds to the relevant portion of a final TCR transcript, again drawing on IMGT nomenclature. There are four valid options: V-REGION, J-REGION, EX1+EX2+EX3+EX4 (constant region), or
L-PART1+L-PART2 (leader sequence). Some things to remember when using custom sequences:

-  Functional leader sequences usually have lengths that are multiples of 3. They don't need to be, but if they're not the V gene will need to account for it to maintain the reading frame.
-  The 3' nucleotide of the J gene is the first nucleotide of the first codon of the constant region.
-  Constant regions in default settings are trimmed by the script to run up to the codon just before the first stop codon (as occur in EX4UTR exons of TRAC and TRDC). This is not required, and stop codons can be left in if desired, but care must be taken if the intention is to use ``thimble`` or ``gui-stitchr`` with these genes to make bicistronic expression constructs. It's recommended to leave stop codons off any constant regions added to additional-genes.fasta, and then provide them in ``thimble`` instead as needed.
-  Most of the gene sequence and format checks cannot be applied, so extra care must be taken to ensure input genes are valid. For instance, using the ``-xg`` flag automatically sets the ``-sc`` flag,
   which skips the usual constant region frame check (as ``stitchr`` doesn't know what frame is intended, see below).
-  Extra genes added via the additional-genes.fasta file are supplemented to the working dictionaries in ``stitcher`` *after* IMGT gene sequences are read in; any extra genes with the same gene name/allele combination as one already in the IMGT dataset will    overwrite the default sequence. If you wish to use both in the same rearrangement or ``thimble`` run, use novel naming in the input FASTA file - e.g. the example constant regions added have 'm' and 'h' prefixes, denoting their human or mouse origin, but any chance to ensure unique names will work.

Skipping constant region checks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the default loci covered, ``stitchr`` has a constant region frame-checking function that uses known correctly-translated sequences to infer the right frame (and where appropriate, placement of endogenous stop codons). If you wish to override these checks for some reason (most likely if you're manually creating your own non-standard or engineered constant region sequences) then you can get the ``-sc / --skip_c_checks`` flag in the command line. Under these circumstances, ``stitchr`` will instead determine the correct frame of the C terminal domain by finding the one with the longest stretch of amino acids before hitting a stop codon. Note that this is less reliable and slower than using the pre-computed motif files. This feature will also **only activate if the gene name of the relevant constant region is not found in the C-region-motifs.tsv file for that species**.

If for some reason users which to skip the C region checks (using the automatically inferred translation frame) for a gene that already is covered in the pre-generated motifs file, they should add a renamed variant of that gene to the ``additional-genes.fasta`` file, and use the ``-xg`` extra genes flag. Note that using the ``-xg`` flag will automatically set the ``-sc / --skip_c_checks`` on.

