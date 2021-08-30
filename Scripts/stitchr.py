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


__version__ = '0.8.1'
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

    parser.add_argument('-aa', '--aa', required=False, type=str,
                help='Partial amino acid sequence, if known. Optional. Can be used to check stitching success')

    parser.add_argument('-n', '--name', required=False, type=str,
                help='Name for TCR sequence. Optional. Will be added to output FASTA header.')

    parser.add_argument('-sl', '--seamless', action='store_true', required=False, default=False,
                help='Optional flag to integrate known nucleotide sequences seamlessly. \n'
                     'NB: nucleotide sequences covering the CDR3 junction with additional V gene context required.')

    parser.add_argument('-5p', '--5_prime_seq', required=False, type=str, default='',
                help='Optional sequence to add to the 5\' of the output sequence (e.g. a Kozak sequence).')

    parser.add_argument('-3p', '--3_prime_seq', required=False, type=str, default='',
                help='Optional sequence to add to the 3\' out the output sequence (e.g. a stop codon).')

    parser.add_argument('-xg', '--extra_genes', action='store_true', required=False, default=False,
                help='Optional flag to use additional (non-database/-natural) '
                     'sequences in the \'additional-genes.fasta\' file.')

    parser.add_argument('-cu', '--codon_usage', required=False, type=str,
                help='Path to a file of Kazusa-formatted codon usage frequencies. Optional.')

    parser.add_argument('-jt', '--j_warning_threshold', required=False, type=int, default=3,
                help='J gene substring length warning threshold. Default = 3. '
                     'Decrease to get fewer notes on short J matches.')

    parser.add_argument('-sc', '--skip_c_checks', action='store_true', required=False, default=False,
                help='Optional flag to skip usual constant region gene checks.')

    return parser.parse_args()


def stitch(specific_args, locus, tcr_info, functionality, partial_info, codon_dict, j_warning_threshold):
    """
    Core function, that performs the actual TCR stitching
    :param specific_args: basic input arguments of a given rearrangement (e.g. V/J/CDR3)
    :param locus: which chain is this looking at, i.e. TRA or TRB
    :param tcr_info: sequence data for the alleles of a specific locus read in from IMGT data
    :param functionality: predicted functionality of different TCR genes, as according to IMGT
    :param partial_info: genes filtered out from input TCR data on account of being incomplete in the database
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
            warnings.warn("No allele specified for " + gene + " - defaulting to *01. ")
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
                    warnings.warn("User specified leader sequence is not evenly divisible by 3 - "
                                  "stitched TCR frame will likely be wrong. ")

                # ...and jump ahead to the stitching (skipping irrelevant TCR gene checks)
                continue

            raise ValueError("Error: " + gene +
                             " is not found in the IMGT data for this chain/species. Please check your gene name. ")

        # And if it does, check whether or not the listed allele has a present value
        if allele not in tcr_info[regions[r]][gene]:

            if partial_info[gene][allele]:
                warnings.warn("Cannot use " + r.upper() + " gene " + specific_args[r] + " as it is classed as \""
                              + partial_info[gene][allele] + "\": attempting prototypical allele (" + gene + "*01). ")
            else:
                warnings.warn("Cannot find " + r.upper() + " gene " + specific_args[r] +
                              ": attempting prototypical allele (" + gene + "*01). ")
            allele = '01'

        # TODO add feature: if gene has >1 allele (especially if with coding diffs) add a warning to the output?

        func_err_base = "Warning: gene " + gene + '*' + allele + " has a IMGT-assigned functionality of \'" \
                          + functionality[gene][allele] + "\', "

        # Check functionality
        if fxn.strip_functionality(functionality[gene][allele]) != 'F':
            warnings.warn(func_err_base + "and thus may not express or function correctly. ")

        if allele in tcr_info[regions[r]][gene]:
            done[r] = tcr_info[regions[r]][gene][allele]
            used_alleles[r] = gene + '*' + allele

            # Special check to account for IMGT not including the 3' terminal J residue if allele from cDNA
            if functionality[gene][allele] == '(F)' and r == 'j':

                cdna_j_err = func_err_base + "meaning it was only detected in cDNA, and thus IMGT doesn't record its" \
                                             " 3\' terminal nucleotide - "

                # Currently no (F) prototypical alleles for human/mouse, but let's plan ahead
                if allele == '01':
                    raise IOError(cdna_j_err + " unable to fix as this is the prototypical allele (*01). ")
                # Otherwise, use the 01 terminal residue
                else:
                    done[r] = done[r] + tcr_info[regions[r]][gene]['01'][-1]
                    warnings.warn(cdna_j_err + "substituting the *01 allele terminal base to maintain reading frame. ")

        else:
            raise ValueError("Cannot find TCR sequence data for "
                             + r.upper() + " gene: " + gene + '*' + allele + ". ")

    # Get information about the C-terminal residue of the CDR3 *should* be, given that J gene
    j_residue_exceptions, low_confidence_js = fxn.get_j_exception_residues(specific_args['species'])

    # Throw a warning if the J gene is one in which the C-terminal residue cannot be confidently identified
    if used_alleles['j'] in low_confidence_js:
        warnings.warn("Warning: " + used_alleles['j'] + " has a \'low confidence\' CDR3-ending motif. ")

    # Then determine whether CDR3 has been provided in amino or nucleic acid form
    if fxn.dna_check(specific_args['cdr3']):
        input_type = 'nt'
        specific_args['cdr3_nt'] = specific_args['cdr3']
        warnings.warn("CDR3 junction provided as DNA sequence: \'" + specific_args['cdr3_nt'] + '\'. ')

    else:
        input_type = 'aa'

    # Determine whether seamless integrated is requested (in platform independent way)
    seamless = False
    if 'seamless' in specific_args:
        if specific_args['seamless']:
            seamless = True

    # Get the germline encoded bits
    n_term_nt_raw = done['l'] + done['v']
    c_term_nt_raw = done['j'] + done['c']

    # Run the appropriate form of non-templated integration
    # First test eligibility for seamless integration
    if input_type == 'nt' and seamless:

        warnings.warn(
            "Seamless option selected: stitched sequence may not be accurate is nucleotide sequence provided is too "
            "short or contains polymorphisms or errors relative to the chosen genes/alleles near the edges. ")

        # Optimistic warning to ensure additional sequence provided for seamless stitching
        if fxn.translate_nt(specific_args['cdr3_nt'][:3]) == 'C':
            warnings.warn("Cys residue detected in first codon of provided seamless CDR3 junction: "
                          "note that seamless stitching requires additional sequence 5\' of the start of the CDR3 "
                          "(disregard this warning if this is not the CDR3 starting residue). ")

        n_term_nt_trimmed, v_overlap = fxn.find_v_overlap(n_term_nt_raw, specific_args['cdr3_nt'])

        # Check for suspiciously short V gene overlaps
        if len(v_overlap) < 10:
            warnings.warn("Only short V gene overlap detected (" + v_overlap + ") for seamless stitching. " )

            # Most common cause = unexpected polymorphism (e.g. SNP or PCR error) in the 5' of the padding sequence
            # Try the overlap search again, starting from one position upstream of the previous match
            n_term_nt_trimmed_2, v_overlap_2 = fxn.find_v_overlap(n_term_nt_raw,
                                                                  specific_args['cdr3_nt'][len(v_overlap)+1:])

            if len(v_overlap_2) > 10:
                warnings.warn("A longer (" + str(len(v_overlap_2)) + " nt) overlap was found after trimming the "
                              "short overlap +1 off, and stitching continued (NB: presumed SNP or PCR error). ")

                n_term_nt_trimmed = n_term_nt_trimmed_2[:-(len(v_overlap) + 1)]
                v_overlap = specific_args['cdr3_nt'][:len(v_overlap) + len(v_overlap_2) + 1]

            else:
                raise ValueError("No longer overlap was found even after trimming that short overlap + 1. Please check "
                                 "the V gene call and that the CDR3 padding sequence doesn't contain polymorphisms. ")

        c_term_nt_trimmed = fxn.find_j_overlap(specific_args['cdr3_nt'][len(v_overlap):], c_term_nt_raw)
        stitched_nt = n_term_nt_trimmed + specific_args['cdr3_nt'] + c_term_nt_trimmed

        # Use the tidy_c_term functionality to frame check/trim excess
        stitched_nt, stitched_trans = fxn.tidy_c_term(stitched_nt, locus, specific_args['species'], False)

        # Catch more 5' SNP errors: if there's a SNP in the edges of the contextual padding this can cause and indel,
        # ... which means tidy_c_term will trim some nt from the 5' of the gene, which we can use to trigger an IOError
        if not stitched_nt.startswith(done['l']):
            raise ValueError("An indel has been detected during seamless stitching, which is usually caused by "
                             "polymorphisms in the padding sequence relative to the genes selected: please either "
                             "ensure selected alleles are correct or provide more context beyond any polymorphisms. ")

    # Otherwise run regular amino-acid based germline determination
    else:

        # If an exact nucleotide junction is provided, first translate
        if input_type == 'nt':
            specific_args['cdr3'] = fxn.translate_nt(specific_args['cdr3_nt'])

            # Frame check (if not providing extra context for seamless integration)
            if len(specific_args['cdr3_nt']) % 3 != 0 and not seamless:
                warnings.warn("Warning: length of CDR3 DNA sequence provided is not evenly divisible by 3 "
                              "and seamless stitching not selected: stitched TCR frame will likely be wrong. ")

        # And check that users haven't asked for amino acid/seamless, which won't work
        elif seamless:
            raise IOError("The seamless option has been selected, yet provided CDR3-containing sequence is not DNA. ")

        # Get codon data, and use to check that there's no unexpected characters in the CDR3
        if len([x for x in list(set([x for x in specific_args['cdr3']])) if x not in list(codon_dict.keys())]) > 0:
            raise ValueError("Unexpected character in CDR3 string. "
                             "Please use only one-letter standard amino acid designations. ")

        # Then check the C-terminus of the CDR3 has an appropriate residue
        # (putting the default F in the dict if not there)
        if used_alleles['j'] not in j_residue_exceptions:
            j_residue_exceptions[used_alleles['j']] = 'F'

        if specific_args['cdr3'][-1] != j_residue_exceptions[used_alleles['j']]:
            warnings.warn("CDR3 provided does not end with the expected residue for this J gene (" +
                j_residue_exceptions[used_alleles['j']] + "). Deletion this far in to the J is extremely unlikely. ")

        # Tidy up the germline edges to be coding for whole codons without any remainders
        n_term_nt_inframe, n_term_aa = fxn.tidy_n_term(n_term_nt_raw)
        c_term_nt_inframe, c_term_aa = fxn.tidy_c_term(c_term_nt_raw, locus,
                                                       specific_args['species'], specific_args['skip_c_checks'])

        # Figure out where the AA CDR3 will slot in: look at the CDR3 edges & see how much overlap needs to be removed
        # Start with 4 residues chunks, move from end of V gene up to 10 residues in (very generous deletion allowance)
        n_term_nt_trimmed, cdr3_n_offset = fxn.determine_v_interface(specific_args['cdr3'],
                                                                     n_term_nt_inframe, n_term_aa)
        c_term_nt_trimmed, cdr3_c_end = fxn.determine_j_interface(specific_args['cdr3'][cdr3_n_offset:],
                                                                  c_term_nt_inframe, c_term_aa,
                                                                  len(done['j']), j_warning_threshold)

        # Generate the non-templated sequences using either supplied nucleotides or common codons established earlier
        if input_type == 'nt':
            non_templated_nt = specific_args['cdr3_nt'][cdr3_n_offset * 3:(cdr3_n_offset+cdr3_c_end) * 3]
        else:
            non_templated_aa = specific_args['cdr3'][cdr3_n_offset:cdr3_n_offset+cdr3_c_end]
            non_templated_nt = fxn.rev_translate(non_templated_aa, codon_dict)

        # Then finally stitch all that info together and output!
        stitched_nt = n_term_nt_trimmed + non_templated_nt + c_term_nt_trimmed

    # If optional 5'/3' sequences are specified, add them to the relevant place
    if '5_prime_seq' in specific_args:
        stitched_nt = specific_args['5_prime_seq'] + stitched_nt
        # Translation offset allows simple translation of output NT without having to figure out the frame
        transl_offset = 3 - (len(specific_args['5_prime_seq']) % 3)
    else:
        transl_offset = 0

    if '3_prime_seq' in specific_args:
        stitched_nt += specific_args['3_prime_seq']

    # Then finally stitch all that info together and output!
    out_bits = [specific_args['name'], used_alleles['v'], used_alleles['j'],
                used_alleles['c'], specific_args['cdr3'], used_alleles['l'] + '(L)']

    # TODO add information to output header if additional 5'/3' sequences specified?
    return out_bits, stitched_nt, transl_offset


regions = {'v': 'V-REGION', 'j': 'J-REGION', 'c': 'EX1+EX2+EX3+EX4', 'l': 'L-PART1+L-PART2'}
gene_types = list(regions.values())

if __name__ == '__main__':

    # Get input arguments, determine the TCR chain in use, get codon table, then load the IMGT data in
    fxn.check_scripts_dir()
    input_args, chain = fxn.sort_input(vars(args()))
    codons = fxn.get_optimal_codons(input_args['codon_usage'], input_args['species'])
    imgt_dat, tcr_functionality, partial = fxn.get_imgt_data(chain, gene_types, input_args['species'])

    if 'extra_genes' in input_args:
        if input_args['extra_genes']:
            imgt_dat, tcr_functionality = fxn.get_additional_genes(imgt_dat, tcr_functionality)
            input_args['skip_c_checks'] = True
    else:
        input_args['skip_c_checks'] = False

    out_list, stitched, offset = stitch(input_args, chain, imgt_dat, tcr_functionality, partial, codons,
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
