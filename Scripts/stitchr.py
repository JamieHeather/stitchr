#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
stitchr.py

Takes command line V, J and (amino acid) CDR3 information
to generate a full length coding nucleotide TCR sequence.
Can be used for TCR vector design, and other purposes.

"""

import functions as fxn
import argparse
import sys


__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description="stiTChR v" + str(__version__) + '\n' + \
                    ": Stitch together a coding TCR nucleotide sequence from V, J, and CDR3 info.\n" + \
                    "Use IMGT gene names, and include terminal CDR3 residues (C/F).\n\n" + \
                    "E.g. \'python stitchr.py -v TRBV20-1 -j TRVJ1-2 -cdr3 CASWHATEVERF\'")

    # Input and output options
    parser.add_argument('-v', '--v', required=True, type=str,
                help='V gene name. Required. Specific allele not required, will default to prototypical (\*01)')

    parser.add_argument('-j', '--j', required=True, type=str,
                help='J gene name. Required. Specific allele not required, will default to prototypical (\*01)')

    parser.add_argument('-cdr3', '--cdr3', required=True, type=str,
                help='CDR3 amino acid sequence. Required. Must include terminal residues (C/F)')

    parser.add_argument('-c', '--c', required=False, type=str,
                help='Constant gene. Optional. Specific allele not required, will default to prototypical (\*01)/TRBC1')

    parser.add_argument('-l', '--l', required=False, type=str,
                help='Leader region. Optional. Will default to the appropriate V gene.')

    parser.add_argument('-cu', '--codon_usage', required=False, type=str, default='../Data/kuzusa-human.txt',
                help='Path to a file of Kuzusa-formatted codon usage frequencies. Optional.')

    parser.add_argument('-aa', '--aa', required=False, type=str,
                help='Partial amino acid sequence, if known. Optional. Can be used to check stitching success')

    parser.add_argument('-n', '--name', required=False, type=str,
                help='Name for TCR sequence. Optional. Will be added to output FASTA header.')

    return parser.parse_args()


if __name__ == '__main__':

    # TODO move all this to one large bracketing function?
    # Get input arguments, determine the TCR chain in use, then load the IMGT data in
    fxn.check_scripts_dir()
    input_args, chain = fxn.sort_input(vars(args()))

    regions = {'v': 'V-REGION', 'j': 'J-REGION', 'c': 'EX1+EX2+EX3+EX4', 'l': 'L-PART1+L-PART2'}
    gene_types = regions.values()
    imgt_dat, functionality = fxn.get_imgt_data(chain, gene_types)

    # Then find each of the appropriate sequences
    done = {}
    for r in regions:

        if '*' in input_args[r]:
            gene, allele = input_args[r].split('*')
            if allele not in imgt_dat[regions[r]][gene]:
                print "\tCannot find", r.upper(), "gene", input_args[r] + \
                                                          ": attempting prototypical allele (" + gene + "*01)"
                allele = '01'
        else:
            gene = input_args[r]
            allele = '01'

        if functionality[gene][allele] != 'F':
            print "Warning: gene", gene + '*' + allele, "has a IMGT-assigned functionality of \'" + \
                functionality[gene][allele] + "\', and thus may not express or function correctly."

        if allele in imgt_dat[regions[r]][gene]:
            done[r] = imgt_dat[regions[r]][gene][allele]
            vars()[r + '_used'] = gene + '*' + allele

        else:
            print "Cannot find TCR sequence data for", r, "gene:", gene + '*' + allele + ". Aborting run"
            sys.exit()

    # Get the germline encoded bits
    n_term_nt, n_term_aa = fxn.tidy_n_term(done['l'] + done['v'])
    c_term_nt, c_term_aa = fxn.tidy_c_term(done['j'] + done['c'], chain)

    # Then figure out where the CDR3 will slot in - look at the CDR3 edges to see how much overlap needs to be removed
    # Start with 4 residues chunks, move from end of V gene up to 10 residues in (very generous deletion allowance)
    n_term_nt_trimmed, cdr3_n_offset = fxn.determine_v_interface(input_args['cdr3'], n_term_nt, n_term_aa)
    c_term_nt_trimmed, cdr3_c_end = fxn.determine_j_interface(input_args['cdr3'], c_term_nt, c_term_aa)

    # Get the most frequent codons for each residue and generate the non-templated sequence
    codons = fxn.get_optimal_codons(input_args['codon_usage'])
    non_templated_aa = input_args['cdr3'][cdr3_n_offset:cdr3_c_end]
    non_templated_nt = ''.join([codons[x] for x in non_templated_aa])

    # Then finally stitch all that info together and output!
    stitched = n_term_nt_trimmed + non_templated_nt + c_term_nt_trimmed
    out_str = '-'.join([input_args['name'], v_used, j_used, c_used, input_args['cdr3'], 'leader', l_used])

    # TODO add a write to file option?
    print '----------------------------------------------------------------------------------------------'
    print fxn.fastafy('nt-' + out_str, stitched)
    print fxn.fastafy('aa-' + out_str, fxn.translate_nt(stitched))

    # If a known/partial amino acid sequence provided, ensure they match up with a quick printed alignment
    if input_args['aa']:
        from Bio import pairwise2
        from Bio.pairwise2 import format_alignment
        alignments = pairwise2.align.globalxx(input_args['aa'], fxn.translate_nt(stitched))
        for i in range(0, 600, 60):
            print '\n'
            if i > len(alignments[0][0]):
                break
            for y in [x[i:i+60] for x in format_alignment(*alignments[0]).split('\n')[:3]]:
                print y
