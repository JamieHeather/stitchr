#!/usr/bin/python
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

"""

import functions as fxn
import os
import sys

__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


all_files = os.listdir(fxn.imgt_data_dir)
if 'imgt-data.fasta' not in all_files:
    print "Error: imgt-data.fasta file not detected. Please generate and place in the Data directory."
    sys.exit()

with open(fxn.imgt_data_dir + 'imgt-data.fasta', 'rU') as in_file, \
        open(fxn.imgt_data_dir + 'TRA.fasta', 'w') as TRA, \
        open(fxn.imgt_data_dir + 'TRB.fasta', 'w') as TRB:
    for fasta_id, seq, blank in fxn.read_fa(in_file):
        bits = fasta_id.split('|')
        if 'TRA' in bits[1]:
            TRA.write(fxn.fastafy(fasta_id, seq))
        elif 'TRB' in bits[1]:
            TRB.write(fxn.fastafy(fasta_id, seq))
