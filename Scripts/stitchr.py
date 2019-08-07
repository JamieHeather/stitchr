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


__version__ = '0.2.2'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

sys.tracebacklimit = 0


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description="stiTChR v" + str(__version__) + '\n' +
                    ": Stitch together a coding TCR nucleotide sequence from V, J, and CDR3 info.\n" +
                    "Use IMGT gene names, and include terminal CDR3 residues (C/F).\n\n" +
                    "E.g. \'python stitchr.py -v TRBV20-1 -j TRVJ1-2 -cdr3 CASWHATEVERF\'")

    # Input and output options
    parser.add_argument('-v', '--v', required=True, type=str,
                help='V gene name. Required. Specific allele not required, will default to prototypical (\*01)')

    parser.add_argument('-j', '--j', required=True, type=str,
                help='J gene name. Required. Specific allele not required, will default to prototypical (\*01)')

    parser.add_argument('-cdr3', '--cdr3', required=True, type=str,
                help='CDR3 amino acid sequence. Required. Must include terminal residues (C/F)')

    parser.add_argument('-s', '--species', required=False, type=str, default='human',
                help='Species. Optional. Only available default options are \'human\' or \'mouse\'. Default = human')

    parser.add_argument('-c', '--c', required=False, type=str,
                help='Constant gene. Optional. Specific allele not required, will default to prototypical (\*01)/TRBC1')

    parser.add_argument('-l', '--l', required=False, type=str,
                help='Leader region. Optional. Will default to the appropriate V gene.')

    parser.add_argument('-cu', '--codon_usage', required=False, type=str,
                help='Path to a file of Kazusa-formatted codon usage frequencies. Optional.')

    parser.add_argument('-aa', '--aa', required=False, type=str,
                help='Partial amino acid sequence, if known. Optional. Can be used to check stitching success')

    parser.add_argument('-n', '--name', required=False, type=str,
                help='Name for TCR sequence. Optional. Will be added to output FASTA header.')

    return parser.parse_args()


if __name__ == '__main__':

    # TODO move all this to one large bracketing function?
    # Get input arguments, determine the TCR chain in use, get codon table, then load the IMGT data in
    fxn.check_scripts_dir()
    input_args, chain, codons = fxn.sort_input(vars(args()))

    regions = {'v': 'V-REGION', 'j': 'J-REGION', 'c': 'EX1+EX2+EX3+EX4', 'l': 'L-PART1+L-PART2'}
    gene_types = regions.values()
    imgt_dat, functionality = fxn.get_imgt_data(chain, gene_types, input_args['species'])

    # Then find each of the appropriate sequences
    done = {}
    for r in regions:

        if '*' in input_args[r]:
            gene, allele = input_args[r].split('*')

        else:
            gene = input_args[r]
            allele = '01'

        # Check this gene exists
        if gene not in imgt_dat[regions[r]]:
            raise ValueError(gene + " is not found in the IMGT data for this species. Please check your gene name.")

        # And if it does, check whether or not the listed allele has a present value
        if allele not in imgt_dat[regions[r]][gene]:
            print "\tCannot find", r.upper(), "gene", input_args[r] + \
                                                      ": attempting prototypical allele (" + gene + "*01)"
            allele = '01'

        # Check functionality
        if functionality[gene][allele] != 'F':
            print "Warning: gene", gene + '*' + allele, "has a IMGT-assigned functionality of \'" + \
                functionality[gene][allele] + "\', and thus may not express or function correctly."

        if allele in imgt_dat[regions[r]][gene]:
            done[r] = imgt_dat[regions[r]][gene][allele]
            vars()[r + '_used'] = gene + '*' + allele

        else:
            raise ValueError("Cannot find TCR sequence data for " + r + " gene: " + gene + '*' + allele +
                             ". Aborting run")

    # Get information about the C-terminal residue of the CDR3 *should* be, given that J gene
    j_residue_exceptions, low_confidence_js = fxn.get_j_exception_residues(input_args['species'])

    # Throw a warning if the J gene is one in which the C-terminal residue cannot be confidently identified
    if j_used in low_confidence_js:
        print "Warning:", j_used, "has a \'low confidence\' CDR3-ending motif"

    # TODO allow users to force ignore the J CDR3 terminal residue check?
    # Then check the C-terminus of the CDR3 has an appropriate residue (putting the default F in the dict if not there)
    if j_used not in j_residue_exceptions:
        j_residue_exceptions[j_used] = 'F'

    if input_args['cdr3'][-1] != j_residue_exceptions[j_used]:
        raise ValueError("CDR3 provided does not end with the expected residue for this J gene (" +
        j_residue_exceptions[j_used] + ")\nDeletion this far in to the J is extremely unlikely. Aborting.")

    # Get the germline encoded bits
    n_term_nt, n_term_aa = fxn.tidy_n_term(done['l'] + done['v'])
    c_term_nt, c_term_aa = fxn.tidy_c_term(done['j'] + done['c'], chain, input_args['species'])

    # Then figure out where the CDR3 will slot in - look at the CDR3 edges to see how much overlap needs to be removed
    # Start with 4 residues chunks, move from end of V gene up to 10 residues in (very generous deletion allowance)
    n_term_nt_trimmed, cdr3_n_offset = fxn.determine_v_interface(input_args['cdr3'], n_term_nt, n_term_aa)
    c_term_nt_trimmed, cdr3_c_end = fxn.determine_j_interface(input_args['cdr3'], c_term_nt, c_term_aa)

    # Generate the non-templated sequences using common codons established earlier
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
    if 'aa' in input_args:
        from Bio import pairwise2
        from Bio.pairwise2 import format_alignment
        alignments = pairwise2.align.globalxx(input_args['aa'], fxn.translate_nt(stitched))
        for i in range(0, 600, 60):
            print '\n'
            if i > len(alignments[0][0]):
                break
            for y in [x[i:i+60] for x in format_alignment(*alignments[0]).split('\n')[:3]]:
                print y

    # TODO incorporate table input/output? With options for 2A/stop/kozak
    # TODO add 'strain' option for mice? Could allow automatic allele selection
    # TODO incorporate a DCR like amino acid check? Could input partial protein sequence, infer genes/CDR3, then get nt
