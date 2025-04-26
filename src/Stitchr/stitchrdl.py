# -*- coding: utf-8 -*-

"""
stitchrdl.py

Download IMGT/GENE-DB data for stitching, for a given species
"""

from . import stitchrfunctions as fxn
from . import __version__
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


__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

sys.tracebacklimit = 0  # comment when debugging

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

    parser.add_argument('-na', '--novel_alleles', required=False, action='store_true', default=False,
                        help="Option to download collated published novel TCR alleles "
                             "to the additional-genes.fasta reference.")

    parser.add_argument('-ns', '--novel_studies_threshold', required=False, type=int, default=2,
                        help="If using the -na/--novel_alleles tag, this optional field allows specifying the number of"
                             "studies an allele must be found in to be included. Default = 2.")

    parser.add_argument('-nd', '--novel_donors_threshold', required=False, type=int, default=3,
                        help="If using the -na/--novel_alleles tag, this optional field allows specifying the number of"
                             "donors (across all studies) an allele must be found in to be included. Default = 3.")

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
        # Check if entry meets the required criteria
        bits = new_entry.split('|')
        imgt_format, gene_allele, inferrable_type, stated_type = [False] * 4
        inferrable_type, stated_type = False, False
        if len(bits) == 16:
            imgt_format = True
            if gene_allele_check(bits[1]):
                gene, allele = bits[1].split('*')
                gene_allele = True
                if len(gene) > 3:
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

    # Then archive the current version of the additional-genes file,
    # to minimise the chance of losing sequences (and thus losing reproducibility)...
    shutil.move(fxn.additional_genes_file,
                os.path.join(fxn.data_dir, 'archived-' + fxn.today() + '-additional-genes.fasta'))

    # ... before writing out the modified additional-genes.fasta file in its place
    with open(fxn.additional_genes_file, 'w') as out_file:
        for entry in current_fa:
            out_file.write(fxn.fastafy(entry, current_fa[entry]))


def download_novel_alleles(study_threshold, donor_threshold):
    """
    :param study_threshold: int of # of studies a given allele must be found in to be retained
    :param donor_threshold: int of # of donors a given allele must be found in to be retained, summed across all studies
    This function grabs novel alleles from the repo where I collate them, and reads them into the additional-genes file
    """

    import pandas as pd
    import requests

    novel_repo_url = 'https://api.github.com/repos/JamieHeather/novel-tcr-alleles/contents/'
    summary_prefix = 'novel-TCR-alleles-'

    # Identify the current novel-tcr-allele file
    response = requests.get(novel_repo_url, headers={})
    novel_file_name = ''
    if response.status_code == 200:
        files = response.json()
        matching_files = [x for x in files if x['name'].startswith(summary_prefix)]
        if len(matching_files) == 1:
            novel_file_name = matching_files[0]['name']
            novel_file_url = matching_files[0]['download_url']

    if not novel_file_name:
        raise IOError("Unable to locate a suitable summary novel TCR allele file name in the GitHub repository. ")

    # If found, download it
    tsv_path = os.path.join(fxn.data_dir, novel_file_name)
    response = requests.get(novel_file_url, headers={})

    if response.status_code == 200:
        with open(tsv_path, 'wb') as out_file:
            out_file.write(response.content)
    else:
        raise IOError("Failed to download the novel TCR allele data. ")

    # Then read in and whittle down to those entries that meet select criteria
    novel = pd.read_csv(tsv_path, sep='\t')
    novel_out_fasta = os.path.join(fxn.data_dir, novel_file_name.replace('.tsv', '.fasta'))
    with open(novel_out_fasta, 'w') as out_file:
        for row in novel.index:
            row_bits = novel.loc[row]
            if pd.isna(row_bits['Notes']):
                notes = ''
            else:
                notes = str(row_bits['Notes'])
            func = '?'

            if 'Stop codon' in notes:
                func = 'P'

            # Disregard alleles found in fewer than the threshold # of studies
            if row_bits['Number-Datasets-In'] < study_threshold:
                continue
            # Disregard alleles found in fewer than the threshold # of donors
            elif row_bits['Number-Donors-In'] < donor_threshold:
                continue
            # If it's already added to IMGT, skip
            elif isinstance(row_bits['IMGT-ID'], str):
                continue
            # Same for if it's a shorter version
            elif 'Shorter version of IMGT' in notes:
                continue

            # Determine whether the allele has a valid name
            allele_id = ''
            if row_bits['Standard-ID'].startswith('TR'):
                allele_id = row_bits['Standard-ID']
            else:
                # If not, borrow one of the names from the paper(s) where it was discovered
                studies = [x for x in novel if '-Name' in x]
                for study in studies:
                    if not pd.isna(row_bits[study]):
                        allele_id = row_bits[study] + '-' + study.replace('-Name', '')
                        if allele_id[:2] != 'TR' or '*' not in allele_id:
                            allele_id = row_bits['Gene'] + '*' + allele_id

            if not allele_id:
                raise IOError("Unable to determine a possible allele ID for sequence with data:\n" + str(row_bits))

            else:
                header = '|'.join([novel_file_name.replace('.tsv', ''), allele_id, 'Homo sapiens', func,
                                   '', '', str(len(row_bits['Ungapped-Sequence'])) + ' nt',
                                   '', '', '', '', '', '', '', notes, '~' + region_key[allele_id[3]]])
                out_file.write(fxn.fastafy(header, row_bits['Ungapped-Sequence']))

    # Having written out the selected novel alleles FASTA, add that to the additional-genes file
    add_from_fasta(novel_out_fasta)
    print("Successfully added novel alleles from " + novel_file_name.replace('.tsv', ''))


def main():

    in_args = vars(args())
    alt_text_suffix = " to the additional-genes.fasta file, in the Stitchr data directory. "

    if in_args['fasta']:
        print("Reading in sequences from " + in_args['fasta'] + alt_text_suffix)
        add_from_fasta(in_args['fasta'])

    elif in_args['novel_alleles']:
        print("Downloading published novel human TCR alleles from GitHub" + alt_text_suffix)
        download_novel_alleles(in_args['novel_studies_threshold'], in_args['novel_donors_threshold'])

    else:
        # Default behaviour, download from IMGT/GENE-DB via IMGTgeneDL
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
