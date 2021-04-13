#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
stitchr.py

Takes command line V, J and (amino acid) CDR3 information
to generate a full length coding nucleotide TCR sequence.
Can be used for TCR vector design, and other purposes.

"""

import functions as fxn
import argparse
import warnings
import sys

__version__ = '0.6.1'
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

    parser.add_argument('-jt', '--j_warning_threshold', required=False, type=int, default=3,
                help='J gene substring length warning threshold. Default = 3. '
                     'Decrease to get fewer notes on short J matches.')

    parser.add_argument('-5p', '--5_prime_seq', required=False, type=str, default='',
                help='Optional sequence to add to the 5\' of the output sequence (e.g. a Kozak sequence).')

    parser.add_argument('-3p', '--3_prime_seq', required=False, type=str, default='',
                help='Optional sequence to add to the 3\' out the output sequence (e.g. a stop codon).')

    parser.add_argument('-xg', '--extra_genes', action='store_true', required=False, default=False,
                help='Optional flag to use additional (non-database/-natural) '
                     'sequences in the \'additional-genes.fasta\' file.')

    parser.add_argument('-sc', '--skip_c_checks', action='store_true', required=False, default=False,
                help='Optional flag to skip usual constant region gene checks.')

    return parser.parse_args()


def stitch(specific_args, locus, tcr_info, functionality, codon_dict, j_warning_threshold):
    """skip_c_checks
    Core function, that performs the actual TCR stitching
    :param specific_args: basic input arguments of a given rearrangement (e.g. V/J/CDR3)
    :param locus: which chain is this looking at, i.e. TRA or TRB
    :param tcr_info: sequence data for the alleles of a specific locus read in from IMGT data
    :param functionality: predicted functionality of different TCR genes, as according to IMGT
    :param codon_dict: dictionary of which codons to use for which amino acids
    :param j_warning_threshold: int threshold value, if a J substring length match is shorter it will throw a warning
    :return: list of details of the TCR as constructed, plus the stitched together nucleotide sequence
    """

    # Find each of the appropriate sequences
    done = {}
    used_alleles = {}

    for r in regions:

        if '*' in specific_args[r]:
            gene, allele = specific_args[r].split('*')

        else:
            gene = specific_args[r]
            allele = '01'
            # TODO if  there were a 'strain' option for mice, appropriate allele selection would happen here

        # Check this gene exists
        if gene not in tcr_info[regions[r]]:

            # If it's a leader sequence, it might be a user-defined DNA sequence
            if r == 'l' and fxn.dna_check(specific_args['l']):
                # If it is, add that info in...
                used_alleles[r] = 'UserSpecifiedLeader*' + specific_args['l']
                done[r] = specific_args['l']

                # Check it's likely to translate in frame
                if len(specific_args['l']) % 3 != 0:
                    warnings.warn("User specified leader sequence is not evenly divisible by 3. "
                                  "Stitched TCR frame will likely be wrong.")

                # ...and jump ahead to the stitching (skipping irrelevant TCR gene checks)
                continue

            raise ValueError("Error: " + gene +
                             " is not found in the IMGT data for this chain/species. Please check your gene name.")

        # And if it does, check whether or not the listed allele has a present value
        if allele not in tcr_info[regions[r]][gene]:
            warnings.warn("Cannot find " + r.upper() + " gene " + specific_args[r] +
                          ": attempting prototypical allele (" + gene + "*01). ")
            allele = '01'

        # TODO add feature: if gene has >1 allele (especially if with coding diffs) add a warning to the output

        # Check functionality
        if functionality[gene][allele] != 'F':
            warnings.warn("Warning: gene " + gene + '*' + allele + " has a IMGT-assigned functionality of \'" +
                          functionality[gene][allele] + "\', and thus may not express or function correctly. ")

        if allele in tcr_info[regions[r]][gene]:
            done[r] = tcr_info[regions[r]][gene][allele]
            used_alleles[r] = gene + '*' + allele

        else:
            raise ValueError("Error: Cannot find TCR sequence data for "
                             + r.upper() + " gene: " + gene + '*' + allele + ". ")

    # Get information about the C-terminal residue of the CDR3 *should* be, given that J gene
    j_residue_exceptions, low_confidence_js = fxn.get_j_exception_residues(specific_args['species'])

    # Throw a warning if the J gene is one in which the C-terminal residue cannot be confidently identified
    if used_alleles['j'] in low_confidence_js:
        warnings.warn("Warning: " + used_alleles['j'] + " has a \'low confidence\' CDR3-ending motif. ")

    # TODO allow users to force ignore the J CDR3 terminal residue check?
    # Then check the C-terminus of the CDR3 has an appropriate residue (putting the default F in the dict if not there)
    if used_alleles['j'] not in j_residue_exceptions:
        j_residue_exceptions[used_alleles['j']] = 'F'

    if specific_args['cdr3'][-1] != j_residue_exceptions[used_alleles['j']]:
        raise ValueError("Error: CDR3 provided does not end with the expected residue for this J gene (" +
                    j_residue_exceptions[used_alleles['j']] + "). Deletion this far in to the J is extremely unlikely. ")

    # Get the germline encoded bits
    n_term_nt, n_term_aa = fxn.tidy_n_term(done['l'] + done['v'])
    c_term_nt, c_term_aa = fxn.tidy_c_term(done['j'] + done['c'], locus,
                                           specific_args['species'], specific_args['skip_c_checks'])

    # Then figure out where the CDR3 will slot in - look at the CDR3 edges to see how much overlap needs to be removed
    # Start with 4 residues chunks, move from end of V gene up to 10 residues in (very generous deletion allowance)
    n_term_nt_trimmed, cdr3_n_offset = fxn.determine_v_interface(specific_args['cdr3'], n_term_nt, n_term_aa)
    c_term_nt_trimmed, cdr3_c_end = fxn.determine_j_interface(specific_args['cdr3'], c_term_nt, c_term_aa,
                                                              j_warning_threshold)

    # Generate the non-templated sequences using common codons established earlier
    non_templated_aa = specific_args['cdr3'][cdr3_n_offset:cdr3_c_end]
    non_templated_nt = fxn.rev_translate(non_templated_aa, codon_dict)

    # If optional 5'/3' sequences are specified, add them to the relevant place
    if '5_prime_seq' in specific_args:
        n_term_nt_trimmed = specific_args['5_prime_seq'] + n_term_nt_trimmed
        # Translation offset allows simple translation of output NT without having to figure out the frame
        transl_offset = 3 - (len(specific_args['5_prime_seq']) % 3)
    else:
        transl_offset = 0

    if '3_prime_seq' in specific_args:
        c_term_nt_trimmed += specific_args['3_prime_seq']

    # Then finally stitch all that info together and output!
    stitched_nt = n_term_nt_trimmed + non_templated_nt + c_term_nt_trimmed
    out_bits = [specific_args['name'], used_alleles['v'], used_alleles['j'],
                used_alleles['c'], specific_args['cdr3'], used_alleles['l'] + '(L)']

    # TODO add information to output header if additional 5'/3' sequences specified?
    return out_bits, stitched_nt, transl_offset


regions = {'v': 'V-REGION', 'j': 'J-REGION', 'c': 'EX1+EX2+EX3+EX4', 'l': 'L-PART1+L-PART2'}
gene_types = list(regions.values())


if __name__ == '__main__':

    # Get input arguments, determine the TCR chain in use, get codon table, then load the IMGT data in
    fxn.check_scripts_dir()
    input_args, chain, codons = fxn.sort_input(vars(args()))

    imgt_dat, tcr_functionality = fxn.get_imgt_data(chain, gene_types, input_args['species'])

    if 'extra_genes' in input_args:
        if input_args['extra_genes']:
            imgt_dat, tcr_functionality = fxn.get_additional_genes(imgt_dat, tcr_functionality)
            input_args['skip_c_checks'] = True
    else:
        input_args['skip_c_checks'] = False

    out_list, stitched, offset = stitch(input_args, chain, imgt_dat, tcr_functionality, codons,
                                input_args['j_warning_threshold'])
    out_str = '|'.join(out_list)

    print('----------------------------------------------------------------------------------------------')
    print(fxn.fastafy('nt|' + out_str, stitched))
    # Use the offset to 5' pad the stitched sequence with 'N's to make up for non-codon length 5' added sequences
    print(fxn.fastafy('aa|' + out_str, fxn.translate_nt('N' * offset + stitched)))

    # If a known/partial amino acid sequence provided, ensure they match up with a quick printed alignment
    if 'aa' in input_args:
        from Bio import pairwise2
        from Bio.pairwise2 import format_alignment
        alignments = pairwise2.align.globalxx(input_args['aa'], fxn.translate_nt('N' * offset + stitched))
        for i in range(0, 600, 60):
            print('\n')
            if i > len(alignments[0][0]):
                break
            for y in [x[i:i+60] for x in format_alignment(*alignments[0]).split('\n')[:3]]:
                print(y)

    # TODO add 'strain' option for mice? Could allow automatic allele selection
    # TODO incorporate a DCR like amino acid check? Could input partial protein sequence, infer genes/CDR3, then get nt
