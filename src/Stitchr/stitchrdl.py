# -*- coding: utf-8 -*-

"""
stitchrdl.py

Download IMGT/GENE-DB data for stitching, for a given species
"""
from . import stitchrfunctions as fxn
import argparse
import os
import subprocess
import sys
import shutil

# Make sure IMGTgeneDL is installed
try:
    import IMGTgeneDL
except (ImportError, ModuleNotFoundError) as err:
    print("IMGTgeneDL is required to run this script (pip install IMGTgeneDL)")
    sys.exit()


__version__ = '0.2.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

#sys.tracebacklimit = 0  # comment when debugging
# TODO

def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description="stitchrdl v" + str(__version__) + '\n' + ": Process input data for use with stitchr. \n"
                    "See https://github.com/JamieHeather/stitchr and https://doi.org/10.1093/nar/gkac190.")

    # Input and output options
    parser.add_argument('-s', '--species', required=False, type=str, default='HUMAN',
                        help="Species (common name). Optional. Default = HUMAN")

    parser.add_argument('--version', action='version', version=__version__,
                        help="Print current stitchr version.")

    parser.add_argument('-fa', '--fasta', required=False, type=str, default='',
                        help="Provide a path to a FASTA file to add sequences to the additional-genes.fasta reference.")

    parser.add_argument('--cite', action=fxn.GetCitation, help="Print citation details.", nargs=0)

    return parser.parse_args()


region_key = {'L-': 'LEADER', 'V-': 'VARIABLE', 'J-': 'JOINING', 'EX': 'CONSTANT', 'CH': 'CONSTANT',
              'V': 'VARIABLE', 'J': 'JOINING', 'C': 'CONSTANT'}


def infer_genetype_genename(gene_name, gene_seq):
    """
    :param gene_name: str of gene name as supplied
    :param gene_seq: str of DNA sequence
    :return: inferred gene region
    NB for 'V' genes sequences <100 nt are presumed to be leaders, >100 nt presumed to be V-REGIONS
    """
    if gene_name[:2] in ['TR', 'IG'] and gene_name[3] in ['V', 'J', 'C']:
        if gene_name[3] == 'V':
            if len(gene_seq) < 100:
                print("Based off FASTA header details and gene length, '" + gene_name +
                      "' was presumed to be a leader sequence.")
                return 'LEADER'
            else:
                print("Based off FASTA header details and gene length, '" + gene_name +
                      "' was presumed to be a variable gene sequence.")
                return 'VARIABLE'
        else:
            return region_key[gene_name[3]]
    else:
        raise IOError("Cannot infer gene type for gene '" + gene_name + "' - please use full IMGT style format. ")


def infer_genetype_whole(fasta_header_bits, fasta_seq):
    """
    :param fasta_header_bits: 16-item list of strs, detailing the contents of an IMGT-style FASTA header
    :param fasta_seq: string of fasta read, to help distinguish V-REGIONs from leader sequences
    :return: str of 'L', 'V', 'J', or 'C' if a gene type can be inferred, or a null value if not
    """
    gene, allele = fasta_header_bits[1].split('*')

    # Preferentially go off the specific region field (as this is the most precise)
    region = fasta_header_bits[4][:2]
    if region:
        if region in region_key:
            return region_key[region]
    # Otherwise try to use the gene name (and sequence
    return infer_genetype_genename(gene, fasta_seq)


def gene_allele_check(given_id):
    if '*' in given_id:
        gene, allele = given_id.split('*')
        if gene and allele:
            return True
        else:
            return False


def add_from_fasta(path_to_fasta):
    """
    :param path_to_fasta: str detailing path to FASTA file containing sequences to add to the additional gene file
    Saves the additional gene FASTA file with new reads from the provided FASTA
    NB: exact ID matches to an existing entry will over-write
    """
    # Read in both existing additional gene file, and the specified sequences to add
    current_fa = fxn.fasta_to_dict(fxn.additional_genes_file)
    to_add_fa = fxn.fasta_to_dict(path_to_fasta)

    # Then go through the new file and add the reads in
    for new_entry in to_add_fa:
        print(new_entry)
        # Check if entry meets the required criteria
        bits = new_entry.split('|')
        imgt_format, gene_allele, inferrable_type, stated_type = [False] * 4
        if len(bits) == 16:
            imgt_format = True
            if gene_allele_check(bits[1]):
                gene, allele = bits[1].split('*')
                gene_allele = True
                if gene[3] in ['V', 'J', 'C'] or bits[4][:2] in ['L-', 'V-', 'J-', 'EX', 'CH']:
                    inferrable_type = True
                if bits[15] in ['~LEADER', '~VARIABLE', '~JOINING', '~CONSTANT']:
                    stated_type = True
        elif len(bits) > 1:
            raise IOError("Cannot automatically read in sequence (need either just sequence ID or full IMGT format) for "
                          + new_entry)

        if not gene_allele and not gene_allele_check(new_entry):
            raise IOError("Failed to read in " + new_entry + ".\n"
                          "Input FASTA entries must have a proper asterisk-delimited 'gene*allele' field (e.g. "
                          "TRAV1-1*01), noting that any arbitrary string can be used in either section. ")

        if imgt_format and gene_allele and stated_type:
            current_fa[new_entry] = to_add_fa[new_entry]
        elif imgt_format and gene_allele and inferrable_type:
            inferred_type = infer_genetype_whole(bits[1], to_add_fa[new_entry])
            current_fa[new_entry[:new_entry.rfind('|')+1] + '~' + inferred_type] = to_add_fa[new_entry]
        elif not imgt_format:
            inferred_type = infer_genetype_genename(new_entry, to_add_fa[new_entry])
            new_id = '|'.join(['stitchrdl-added', new_entry] + [' '] * 14)
            new_id = new_id[:new_id.rfind('|')+1] + '~' + inferred_type
            current_fa[new_id] = to_add_fa[new_entry]

        else:
            raise IOError("Failed to automatically read in entry for '" + new_entry +
                          "' - please try full formatting and try again. ")

    # Over-write the additional-genes.fasta file in place
    with open(fxn.additional_genes_file, 'w') as out_file:
        for entry in current_fa:
            out_file.write(fxn.fastafy(entry, current_fa[entry]))


def main():

    in_args = vars(args())

    if in_args['fasta']:
        print("Reading in sequences from " + in_args['fasta'] + " to the additional-genes.fasta file, "
                                                                "in the Stitchr data directory. ")
        add_from_fasta(in_args['fasta'])

    else:

        species = in_args['species'].upper()

        try:
            cmd = 'IMGTgeneDL -m stitchr -n -s ' + species
            subprocess.call(cmd, shell=True)
        except Exception:
            print("Failed to download sequences from IMGTgeneDL.")

        # Assuming that made the folder appropriately
        if species in os.listdir(os.getcwd()):
            # Delete older entries
            if species in os.listdir(fxn.data_dir):
                print("Note: a directory for this species exists: overwriting.")
                shutil.rmtree(os.path.join(fxn.data_dir, species))
            shutil.move(os.path.join(os.getcwd(), species), os.path.join(fxn.data_dir, species))

        else:
            raise IOError("Failed to produce the expected species\' output directory. ")
