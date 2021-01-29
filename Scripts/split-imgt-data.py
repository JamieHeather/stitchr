#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
split-imgt-data.py

Take IMGT-data.fasta file, and split out the separate alpha/beta chain information to separate files.

Data acquired by:
* Going to IMGT/GENE-DB (http://www.imgt.org/genedb/)
* Selecting 'Species' = 'Homo sapiens' and 'Molecular component' = 'TR'
* Selecting all genes and extracting and combining the following into one file:
    * L-PART1+L-PART2
    * V-REGION
    * J-REGION
    * EX1+EX2+EX3+EX4

Note that EX1+EX2+EX3+EX4 annotations may not exist for all constant regions in all species.
Instead an entry can be manually produced by combining the individual entries (EX1, EX2, EX3, and EX4).
This is what was done for the Mus musculus TRBC genes.
"""

import functions as fxn
import os
import time
from datetime import datetime
import sys
import collections as coll

__version__ = '0.4.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


for species in ['HUMAN', 'MOUSE']:

    # Check whether the appropriate files are present
    species_dir = fxn.data_dir + '/' + species + '/'
    all_files = os.listdir(species_dir)
    if 'imgt-data.fasta' not in all_files:
        print("Error: imgt-data.fasta file not detected for\'", species + \
                "'. Please generate and place it in the appropriate Data subdirectory.")
        sys.exit()

    # If so, check the modification time for the imgt-data.fasta file, assuming that's the last download time
    input_imgt_file = species_dir + 'imgt-data.fasta'
    mod_date = datetime.fromtimestamp(os.path.getmtime(input_imgt_file)).strftime('%Y-%m-%d')

    # Then read through the FASTA and sort into the appropriate chains
    with open(input_imgt_file, 'rU') as in_file, \
            open(species_dir + 'TRA.fasta', 'w') as TRA, \
            open(species_dir + 'TRB.fasta', 'w') as TRB:

        prot = coll.defaultdict(coll.defaultdict)

        for fasta_id, seq, blank in fxn.read_fa(in_file):
            gene, allele = fasta_id.split('|')[1].split('*')

            # NB: TRDV included with TRA genes due to the evidence that even non 'TRAV/DV' genes can recombine with TRAJ
            if 'TRA' in gene or 'TRDV' in gene:
                TRA.write(fxn.fastafy(fasta_id, seq))
            elif 'TRB' in gene:
                TRB.write(fxn.fastafy(fasta_id, seq))

    # Finally log the dates
    log_txt = 'imgt-data.fasta_last_modified ' + mod_date + '\nsplit-imgt-data.py_last_run ' + fxn.today()
    with open(species_dir + 'data-production-date.txt', 'w') as log_file:
        log_file.write(log_txt)
