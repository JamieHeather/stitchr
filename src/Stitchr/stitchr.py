# -*- coding: utf-8 -*-

"""
stitchr.py

Takes command line V, J and (amino acid) CDR3 information
to generate a full length coding nucleotide TCR sequence.
Can be used for TCR vector design, and other purposes.

"""

from . import stitchrfunctions as fxn
from . import __version__
import argparse
import sys
import warnings


__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

sys.tracebacklimit = 0  # comment when debugging
warnings.formatwarning = fxn.custom_formatwarning

def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(
        description="stiTChR v" + str(__version__) + '\n' +
                    ": Stitch together a coding TCR nucleotide sequence from V, J, and CDR3 info.\n"
                    "Use IMGT gene names, and include terminal CDR3 residues (C/F).\n\n"
                    "E.g. \'stitchr -v TRBV20-1*01 -j TRBJ1-2*01 -cdr3 CASWHATEVERF\'.\n"
                    "See https://github.com/JamieHeather/stitchr and https://doi.org/10.1093/nar/gkac190.")

    # Input and output options
    parser.add_argument('-v', '--v', required=True, type=str, default='', help="V gene name. Required. "
                        "Specific allele not required, will default to prototypical (*01)")

    parser.add_argument('-j', '--j', required=True, type=str, default='', help="J gene name. Required. "
                        "Specific allele not required, will default to prototypical (*01)")

    parser.add_argument('-cdr3', '--cdr3', required=True, type=str, default='',
                        help="CDR3 amino acid sequence. Required. Must include terminal residues (e.g. C/F)")

    parser.add_argument('-s', '--species', required=False, type=str, default='HUMAN',
                        help="Species (common name). Optional: see data directory for all possible options. "
                             "Default = HUMAN")

    parser.add_argument('-c', '--c', required=False, type=str, default='', help="Constant gene. Optional. "
                        "Specific allele not required, will default to prototypical (*01). "
                        "See README re: alternative TRGC exon configurations.")

    parser.add_argument('-l', '--l', required=False, type=str, default='', help="Leader region. Optional. "
                        "Will default to match the appropriate V gene.")

    parser.add_argument('-aa', '--aa', required=False, type=str, default='',
                        help="Partial amino acid sequence, if known. Optional. Can be used to check stitching success.")

    parser.add_argument('-n', '--name', required=False, type=str, default='',
                        help="Name for TCR sequence. Optional. Will be added to output FASTA header.")

    parser.add_argument('-sl', '--seamless', action='store_true', required=False, default=False,
                        help="Optional flag to integrate known nucleotide sequences seamlessly. \n NB: "
                        "nucleotide sequences covering the CDR3 junction with additional V gene context required.")

    parser.add_argument('-5p', '--5_prime_seq', required=False, type=str, default='',
                        help="Optional sequence to add to the 5' of the output sequence (e.g. a Kozak sequence).")

    parser.add_argument('-3p', '--3_prime_seq', required=False, type=str, default='',
                        help="Optional sequence to add to the 3' out the output sequence (e.g. a stop codon).")

    parser.add_argument('-xg', '--extra_genes', action='store_true', required=False, default=False,
                        help="Optional flag to use additional (non-database/-natural) sequences in the "
                             "'additional-genes.fasta' file.")

    parser.add_argument('-p', '--preferred_alleles_path', required=False, type=str, default='',
                        help="Path to a file of preferred alleles to use when no allele specified (instead of *01). "
                             "Optional.")

    parser.add_argument('-m', '--mode', required=False, type=str, default='both_fa',
                        help="Standard out output mode. "
                             "Options are 'BOTH_FA' (default), 'AA_FA', 'NT_FA', 'AA', 'NT', 'GB', 'JSON'.")

    parser.add_argument('-cu', '--codon_usage_path', required=False, type=str, default='',
                        help="Path to a file of Kazusa-formatted codon usage frequencies. Optional.")

    parser.add_argument('-jt', '--j_warning_threshold', required=False, type=int, default=3,
                        help="J gene substring length warning threshold. Default = 3. "
                             "Decrease to get fewer notes on short J matches.")

    parser.add_argument('-sc', '--skip_c_checks', action='store_true', required=False, default=False,
                        help="Optional flag to skip usual constant region gene checks.")

    parser.add_argument('-sn', '--skip_n_checks', action='store_true', required=False, default=False,
                        help="Optional flag to skip the usual CDR3 N terminal (conserved 2nd Cys) checks.")

    parser.add_argument('-nl', '--no_leader', action='store_true', required=False, default=False,
                        help="Optional flag to ignore leader sequences (signal peptides).")

    parser.add_argument('-sw', '--suppress_warnings', action='store_true', required=False, default=False,
                        help="Optional flag to suppress warnings.")

    parser.add_argument('--version', action='version', version=__version__,
                        help="Print current stitchr package version.")

    parser.add_argument('--cite', action=fxn.GetCitation, help="Print citation details.", nargs=0)

    parser.add_argument('-dd', '--data_dir', action='version', version=fxn.data_dir,
                        help="Print installed stitchr data directory path.")

    return parser.parse_args()


def stitch(specific_args, tcr_info, functionality, partial_info, codon_dict, j_warning_threshold, preferences,
           c_motifs, j_residues, low_confidence_js):
    """
    Core function, that performs the actual TCR stitching
    :param specific_args: dict, basic input arguments of a given rearrangement (e.g. V/J/CDR3)
    :param tcr_info: dict, germline sequence data for the alleles of a specific locus
    :param functionality: dict, predicted functionality of different TCR genes
    :param partial_info: dict, genes filtered out from input TCR data on account of being incomplete in the database
    :param codon_dict: dictionary of which codons to use for which amino acids
    :param j_warning_threshold: int threshold value, if a J substring length match is shorter it will throw a warning
    :param preferences: nested dict of preferred alleles, like the tcr_info dict but one level shallower
    :param c_motifs: dict, product of get_c_motifs, detailing N terminal C gene translations to find correct ORFs
    :param j_residues: dict, product of get_j_motifs, detailing the conserved F/equivalent residue to look for
    :param low_confidence_js: list of str detailing 'low confidence' J regions which lack an obvious FGXG motif
    :return: stitch_dict of details of the in/in process/constructed TCR,
        including the stitched nt and translation offset (int, 0/1/2)
    """

    stitch_dict = {'in': specific_args, 'used': {}, 'seqs': {}}

    for r in fxn.regions:

        if r == 'l' and specific_args['no_leader']:
            continue

        # First establish what the input gene and allele values are
        gene_allele_in = specific_args[r]
        if '*' not in gene_allele_in:
            gene_allele_in += '*'

        in_gene, in_allele = gene_allele_in.split('*')
        gene, allele = '', ''

        # First check whether the gene exists in the input data for this species
        if in_gene in tcr_info[fxn.regions[r]]:
            gene = in_gene

        else:
            # If it's a leader sequence, it might be a user-defined DNA sequence
            if r == 'l' and fxn.dna_check(specific_args['l']):
                # If it is, add that info in...
                stitch_dict['used'][r] = 'UserSpecifiedLeader*' + specific_args['l']
                stitch_dict['seqs'][r] = specific_args['l']

                # Check it's likely to translate in frame
                if len(specific_args['l']) % 3 != 0:
                    warnings.warn("User specified leader sequence is not evenly divisible by 3 - "
                                  "stitched TCR frame will likely be wrong. ")

                # ...and jump ahead to the stitching (skipping irrelevant TCR gene checks)
                continue

            raise ValueError("Error: a " + fxn.regions[r] + " sequence region has not been found for gene " + gene +
                             " in the IMGT data for this chain/species. Please check your TCR and species data. ")

        # Having established the gene, then need to determine the allele
        if in_allele:

            # If an allele is provided, check whether it exists in the database for this gene and isn't partial
            if in_allele in tcr_info[fxn.regions[r]][gene]:

                if partial_info[gene][allele]:
                    warnings.warn("Cannot use " + gene + "*" + in_allele + " " + fxn.regions[r].lower() + "sequence, "
                                  "as it is classed as '" + partial_info[gene][in_allele] + "' (partial). ")
                else:
                    allele = in_allele

            else:
                warnings.warn("Cannot find the sequence of the requested allele, " + gene + "*" + in_allele +
                              ", for the " + fxn.regions[r].lower() + " sequence in the input FASTA data. ")

        # If no allele supplied (or is supplied but invalid) then use 1) a preferred allele or 2) the prototypical *01
        if not allele:

            warnings.warn("No valid " + fxn.regions[r].lower() + " region allele determined yet for " + gene + ". ")
            allele = '01'

            if preferences:
                if gene in preferences[fxn.regions[r]]:
                    allele = preferences[fxn.regions[r]][gene]
                    warnings.warn("Defaulting to allele *" + allele + " for the " + fxn.regions[r].lower() +
                                  " sequence, as specified in the preferred allele file. ")
                    # NB: we don't have to worry about partial alleles in the preference list, as those are filtered out
                else:
                    warnings.warn("Defaulting to *01, as " + gene + " isn't specified in the preferred allele file for"
                                                                    " the " + fxn.regions[r].lower() + " region. ")
            else:
                warnings.warn("Defaulting to *01 for the " + fxn.regions[r].lower() + " region, "
                              "in the absence of a preferred allele file being specified. ")

            # Check how feasible the defaulted allele is
            if allele == '01':
                # First check this allele actually exists - if not pick one from what's available
                if not tcr_info[fxn.regions[r]][gene][allele]:
                    new_allele = [x for x in tcr_info[fxn.regions[r]][gene] if tcr_info[fxn.regions[r]][gene][x]][0]
                    warnings.warn("No sequence found for " + gene + "*" + allele + " for the " + fxn.regions[r] +
                                  " region. " + gene + "*" + new_allele + " is being used instead - please double "
                                  "check either the intended allele or that this is a suitable replacement. ")
                    allele = new_allele

                # If the prototypical allele does exist, flag a warning if there are other options available
                elif len(tcr_info[fxn.regions[r]][gene]) > 1:
                    # Just check that at least one of those other alleles is not partial
                    for other_allele in tcr_info[fxn.regions[r]][gene]:
                        if other_allele != '01':
                            if not partial_info[gene][other_allele]:
                                warnings.warn("NB: the prototypical '*01' allele is being used for the " +
                                              fxn.regions[r].lower() + " region by default, but other alleles are "
                                                                       "available - consider double checking the right allele is asked for. ")
                                continue

        # Catchall double check both gene and allele are sorted
        if gene and allele:
            stitch_dict['seqs'][r] = tcr_info[fxn.regions[r]][gene][allele]
            stitch_dict['used'][r] = gene + '*' + allele

            # Add '(L)' to leader names
            if r == 'l':
                stitch_dict['used']['l'] = stitch_dict['used']['l'] + '(L)'

            func_err_base = "Warning: gene " + gene + '*' + allele + " has an assigned functionality of \'" \
                            + functionality[gene][allele] + "\', "

            # Check functionality
            if fxn.strip_functionality(functionality[gene][allele]) != 'F':
                warnings.warn(func_err_base + "and thus may not express or function correctly. ")

            # Special check to account for IMGT not including the 3' terminal J residue if allele from cDNA
            if functionality[gene][allele] == '(F)' and r == 'j':

                cdna_j_err = func_err_base + "meaning it was only detected in cDNA, and thus IMGT doesn't record its" \
                                             " 3\' terminal nucleotide - "

                # If this is an (F) call for a prototypical allele we have to skip, as there's no reference to go off
                if allele == '01':
                    raise IOError(cdna_j_err + " unable to fix as this is the prototypical allele (*01). ")
                # Otherwise, use the 01 terminal residue
                else:
                    stitch_dict['seqs'][r] = stitch_dict['seqs'][r] + tcr_info[fxn.regions[r]][gene]['01'][-1]
                    warnings.warn(cdna_j_err + "substituting the *01 allele terminal base to maintain reading frame. ")

        else:
            raise ValueError(
                "Cannot find TCR sequence data for " + r.upper() + " gene: " + in_gene + '*' + in_allele + ". ")

    # Throw a warning if the J gene is one in which the C-terminal residue cannot be confidently identified
    if stitch_dict['used']['j'] in low_confidence_js:
        warnings.warn("Warning: " + stitch_dict['used']['j'] + " has a \'low confidence\' CDR3-ending motif. ")

    # Then determine whether CDR3 has been provided in amino or nucleic acid form
    if fxn.dna_check(specific_args['cdr3']):
        stitch_dict['input_type'] = 'nt'
        stitch_dict['used']['cdr3'] = fxn.translate_nt(specific_args['cdr3'])
        warnings.warn("CDR3 junction provided as DNA sequence: \'" + specific_args['cdr3'] + '\'. ')
    else:
        stitch_dict['input_type'] = 'aa'
        stitch_dict['used']['cdr3'] = specific_args['cdr3']

    # Determine whether seamless integrated is requested (in platform independent way)
    seamless = False
    if 'seamless' in specific_args:
        if specific_args['seamless']:
            seamless = True
    else:
        specific_args['seamless'] = False
        stitch_dict['in']['seamless'] = False

    # Get the germline encoded bits
    if specific_args['no_leader']:
        n_term_nt_raw = stitch_dict['seqs']['v']
    else:
        n_term_nt_raw = stitch_dict['seqs']['l'] + stitch_dict['seqs']['v']
    c_term_nt_raw = stitch_dict['seqs']['j'] + stitch_dict['seqs']['c']

    # Run the appropriate form of non-templated integration
    # First test eligibility for seamless integration
    if stitch_dict['input_type'] == 'nt' and seamless:

        warnings.warn(
            "Seamless option selected: stitched sequence may not be accurate if provided NT sequence is too short"
            " or contains polymorphisms/errors near the recombining edges. ")

        # Optimistic warning to ensure additional sequence provided for seamless stitching
        if fxn.translate_nt(specific_args['cdr3'][:3]) == 'C':
            warnings.warn("Cys residue detected in first codon of provided seamless CDR3 junction: "
                          "note that seamless stitching requires additional sequence 5\' of the start of the CDR3 "
                          "(disregard this warning if this is not the CDR3 starting residue). ")

        n_term_nt_trimmed, v_overlap = fxn.find_v_overlap(n_term_nt_raw, specific_args['cdr3'])

        # Check for suspiciously short V gene overlaps
        if len(v_overlap) < 10:
            warnings.warn("Only short V gene overlap detected (" + v_overlap + ") for seamless stitching. ")

            # Most common cause = unexpected polymorphism (e.g. SNP or PCR error) in the 5' of the padding sequence
            # Try the overlap search again, starting from one position upstream of the previous match
            n_term_nt_trimmed_2, v_overlap_2 = fxn.find_v_overlap(n_term_nt_raw,
                                                                  specific_args['cdr3'][len(v_overlap)+1:])

            if len(v_overlap_2) > 10:
                warnings.warn("A longer (" + str(len(v_overlap_2)) + " nt) overlap was found after trimming the "
                                                                     "short overlap +1 off, and stitching continued (NB: presumed SNP or PCR error). ")

                n_term_nt_trimmed = n_term_nt_trimmed_2[:-(len(v_overlap) + 1)]
                v_overlap = specific_args['cdr3'][:len(v_overlap) + len(v_overlap_2) + 1]

            else:
                raise ValueError("No longer overlap was found even after trimming that short overlap + 1. Please check "
                                 "the V gene call and that the CDR3 padding sequence doesn't contain polymorphisms. ")

        c_term_nt_trimmed = fxn.find_j_overlap(specific_args['cdr3'][len(v_overlap):], c_term_nt_raw)
        stitched_nt = n_term_nt_trimmed + specific_args['cdr3'] + c_term_nt_trimmed

        stitch_dict['seqs']['seamless_insert'] = specific_args['cdr3']
        stitch_dict['used']['seamless_insert'] = 'seamless_insert'

        # Use the tidy_c_term functionality to frame check/trim excess
        stitched_nt, stitched_trans = fxn.tidy_c_term(stitched_nt, False, c_motifs, stitch_dict['used']['c'])

        # Catch more 5' SNP errors: if there's a SNP in the edges of the contextual padding this can cause and indel,
        # ... which means tidy_c_term will trim some nt from the 5' of the gene, which we can use to trigger an error
        if not stitched_nt.startswith(n_term_nt_trimmed[:30]):
            raise ValueError("An indel has been detected during seamless stitching, which is usually caused by "
                             "polymorphisms in the padding sequence relative to the genes selected: please either "
                             "ensure selected alleles are correct or provide more context beyond any polymorphisms. ")

    # Otherwise run regular amino-acid based germline determination
    else:

        # If an exact nucleotide junction is provided, work off the translated product
        if stitch_dict['input_type'] == 'nt':

            # Frame check (if not providing extra context for seamless integration)
            if len(stitch_dict['in']['cdr3']) % 3 != 0 and not seamless:
                # TODO this is the problem to fix
                warnings.warn("Warning: length of CDR3 DNA sequence provided is not evenly divisible by 3 "
                              "and seamless stitching not selected: stitched TCR frame will likely be wrong. ")

        # And check that users haven't asked for amino acid/seamless, which won't work
        elif seamless:
            raise IOError("The seamless option has been selected, yet provided CDR3-containing sequence is not DNA. ")

        # Get codon data, and use to check that there's no unexpected characters in the CDR3
        if len([x for x in list(set([x for x in stitch_dict['used']['cdr3']])) if x not in list(codon_dict.keys())]) > 0:
            raise ValueError("Unexpected character in CDR3 string. "
                             "Please use only one-letter standard amino acid designations. ")

        # Then check the C-terminus of the CDR3 has an appropriate residue
        # (putting the default F in the dict if not there)
        if stitch_dict['used']['j'] not in j_residues:
            j_residues[stitch_dict['used']['j']] = 'F'

        if stitch_dict['used']['cdr3'][-1] != j_residues[stitch_dict['used']['j']]:
            warnings.warn("CDR3 provided does not end with the expected residue for this J gene (" +
                          j_residues[stitch_dict['used']['j']] + "). Deletion this far in to the J is extremely unlikely. ")

        # Tidy up the germline edges to be coding for whole codons without any remainders
        n_term_nt_inframe, n_term_aa = fxn.tidy_n_term(n_term_nt_raw)

        c_term_nt_inframe, c_term_aa = fxn.tidy_c_term(c_term_nt_raw,
                                                       specific_args['skip_c_checks'],
                                                       c_motifs, stitch_dict['used']['c'])

        # Figure out where the AA CDR3 will slot in: look at the CDR3 edges & see how much overlap needs to be removed
        # Start with 4 residues chunks, move from end of V gene up to 10 residues in (very generous deletion allowance)
        n_term_nt_trimmed, cdr3_n_offset = fxn.determine_v_interface(stitch_dict['used']['cdr3'],
                                                                     n_term_nt_inframe, n_term_aa,
                                                                     specific_args['skip_n_checks'])

        c_term_nt_trimmed, cdr3_c_end = fxn.determine_j_interface(stitch_dict['used']['cdr3'][cdr3_n_offset:],
                                                                  c_term_nt_inframe, c_term_aa,
                                                                  len(stitch_dict['seqs']['j']), j_warning_threshold)

        # Generate the non-templated sequences using either supplied nucleotides or common codons established earlier
        if stitch_dict['input_type'] == 'nt':
            non_templated_nt = specific_args['cdr3'][cdr3_n_offset * 3:(cdr3_n_offset+cdr3_c_end) * 3]
        else:
            non_templated_aa = stitch_dict['used']['cdr3'][cdr3_n_offset:cdr3_n_offset+cdr3_c_end]
            non_templated_nt = fxn.rev_translate(non_templated_aa, codon_dict)

        # Then finally stitch all that info together and output!
        stitched_nt = n_term_nt_trimmed + non_templated_nt + c_term_nt_trimmed
        stitch_dict['seqs']['non_templated'] = non_templated_nt
        stitch_dict['used']['non_templated'] = 'non_templated'

    if specific_args['no_leader']:
        stitch_dict['seqs']['v'] = n_term_nt_trimmed
    else:
        stitch_dict['seqs']['v'] = n_term_nt_trimmed[len(stitch_dict['seqs']['l']):]
    stitch_dict['seqs']['c'] = stitch_dict['seqs']['c'][:stitch_dict['seqs']['c'].index(c_term_nt_trimmed[-20:]) + 20]
    stitch_dict['seqs']['j'] = c_term_nt_trimmed[:-len(stitch_dict['seqs']['c'])]
    # TODO potentially find a more elegant way to determine/retain the C start/end positions?

    # If optional 5'/3' sequences are specified, add them to the relevant place
    if specific_args['5_prime_seq']:

        stitched_nt = specific_args['5_prime_seq'] + stitched_nt
        stitch_dict['seqs']['5_prime_seq'] = specific_args['5_prime_seq']
        stitch_dict['used']['5_prime_seq'] = '5_prime_seq'

        # Translation offset allows simple translation of output NT without having to figure out the frame
        stitch_dict['translation_offset'] = 3 - (len(specific_args['5_prime_seq']) % 3)
        if stitch_dict['translation_offset'] == 3:
            stitch_dict['translation_offset'] = 0

    else:
        stitch_dict['translation_offset'] = 0

    if specific_args['3_prime_seq']:
        stitch_dict['used']['3_prime_seq'] = '3_prime_seq'
        stitched_nt += specific_args['3_prime_seq']
        stitch_dict['seqs']['3_prime_seq'] = specific_args['3_prime_seq']

    stitch_dict['out_list'] = [specific_args['name'], stitch_dict['used']['v'], stitch_dict['used']['j'],
                               stitch_dict['used']['c'], specific_args['cdr3']]
    if not specific_args['no_leader']:
        stitch_dict['out_list'].append(stitch_dict['used']['l'])
    stitch_dict['stitched_nt'] = stitched_nt

    # If this is to be output to GenBank format, also determine the location of each substr,
    # to ensure proper plotting of even non-unique string
    if specific_args['mode'] == 'GB':
        stitch_dict['indexes'] = {}
        region_order = ['5_prime_seq', 'l', 'v', 'cdr3', 'seamless_insert', 'j', 'c', '3_prime_seq']

        # Generate the nt sequence of the CDR3 (where possible)
        if not stitch_dict['in']['seamless']:
            translation = fxn.translate_nt('N' * stitch_dict['translation_offset'] + stitch_dict['stitched_nt'])
            cdr3_match = translation.find(stitch_dict['used']['cdr3']) * 3
            stitch_dict['seqs']['cdr3'] = stitch_dict['stitched_nt'][cdr3_match - stitch_dict['translation_offset']:
                                                                     cdr3_match + stitch_dict['translation_offset'] + len(stitch_dict['used']['cdr3']*3)]

        loc = 0
        for region in region_order:
            if region in stitch_dict['seqs'].keys():
                if region == '3_prime_seq':
                    stitch_dict['indexes'][region] = stitched_nt.rfind(stitch_dict['seqs'][region])
                else:
                    index = stitched_nt[loc:].index(stitch_dict['seqs'][region])
                    loc += index
                    stitch_dict['indexes'][region] = loc

        if 'non_templated' in stitch_dict['seqs'].keys():
            stitch_dict['indexes']['non_templated'] = (stitched_nt[:stitch_dict['indexes']['j']].rfind(
                stitch_dict['seqs']['non_templated']))

    return stitch_dict


gene_types = list(fxn.regions.values())


def main():

    # Get input arguments, determine the TCR chain in use, get codon table, then load the IMGT data in
    input_args, chain = fxn.sort_input(vars(args()))
    codons = fxn.get_optimal_codons(input_args['codon_usage_path'], input_args['species'])
    imgt_dat, tcr_functionality, partial = fxn.get_ref_data(chain, gene_types, input_args['species'])

    # Get information about the C-terminal residue of the CDR3 *should* be, given that J gene
    j_res, low_conf_js = fxn.get_j_motifs(input_args['species'])
    # And the motifs required for the correct frame inference and delineation of the constant region sequences
    c_res = fxn.get_c_motifs(input_args['species'])

    if input_args['suppress_warnings']:
        warnings.filterwarnings('ignore')

    if input_args['extra_genes']:
        imgt_dat, tcr_functionality = fxn.get_additional_genes(imgt_dat, tcr_functionality)
        input_args['skip_c_checks'] = True

    if input_args['preferred_alleles_path']:
        preferred_alleles = fxn.get_preferred_alleles(input_args['preferred_alleles_path'],
                                                      gene_types, imgt_dat, partial, chain)
    else:
        preferred_alleles = {}

    with warnings.catch_warnings(record=True) as stitch_warnings:
        stitched = stitch(input_args, imgt_dat, tcr_functionality, partial, codons, input_args['j_warning_threshold'],
                          preferred_alleles, c_res, j_res, low_conf_js)

    out_str = '|'.join(stitched['out_list'])

    # Output the appropriate strings to stdout
    if input_args['mode'] not in ['BOTH_FA', 'AA_FA', 'NT_FA', 'AA', 'NT', 'GB', 'JSON']:
        raise IOError("Unknown output mode detected: " + input_args['mode'] + ". \n"
                      "Should be one of 'BOTH_FA' (default), 'AA_FA', 'NT_FA', 'AA', 'NT', 'GB', or 'JSON'.")

    if '_FA' in input_args['mode']:
        print('----------------------------------------------------------------------------------------------')
        if input_args['mode'] == 'BOTH_FA' or input_args['mode'] == 'NT_FA':
            print(fxn.fastafy('nt|' + out_str, stitched['stitched_nt']))

        if input_args['mode'] == 'BOTH_FA' or input_args['mode'] == 'AA_FA':
            # Use the offset to 5' pad the stitched sequence with 'N's to make up for non-codon length 5' added seqs
            print(fxn.fastafy('aa|' + out_str,
                              fxn.translate_nt('N' * stitched['translation_offset']
                                               + stitched['stitched_nt'])))

    elif input_args['mode'] == 'NT':
        print(stitched['stitched_nt'])

    elif input_args['mode'] == 'AA':
        print(fxn.translate_nt('N' * stitched['translation_offset'] + stitched['stitched_nt']))

    elif input_args['mode'] == 'GB':

        gb_file_name = fxn.get_output_name(stitched)

        # Generate the DEFINITION field
        description = fxn.get_metadata_text(input_args['species'], 'stitchr', __version__)

        if input_args['extra_genes']:
            description += ("Note that this was run using the 'extra genes' flag, and thus additional non-reference "
                            "sequences may have been incorporated. ")

        description += ("\n" + fxn.spacer + "\n" + fxn.spacer +
                        "Produced using the stitchr command:\n" + fxn.spacer + "'" + ' '.join(sys.argv[:]) + "'. ")

        if stitch_warnings and not input_args['suppress_warnings']:
            description += ("\n" + fxn.spacer + "\n" + fxn.spacer +
                            "This rearrangement produced the following warning messages:\n" + fxn.spacer) + ' '.join(
                [str(stitch_warnings[x].message) for x in range(len(stitch_warnings))])

        # And define the individual features to be recorded
        gb_feats = []
        for feat in stitched['indexes']:
            gb_feats.append((stitched['used'][feat], '', [stitched['seqs'][feat]], [stitched['indexes'][feat]]))
        # TODO move stitchr/thimble-to-genbank formatting to own function?

        genbank_params = {'sequence_name': gb_file_name, 'full_sequence': stitched['stitched_nt'],
                          'description': description, 'topology': 'linear', 'features': gb_feats,
                          'save_dir_path': './', 'species': input_args['species'], 'numbering': False,
                          'plot_multi': True, 'division': 'SYN', 'journal': 'stitchr'}

        fxn.output_genbank(**genbank_params)

    elif input_args['mode'] == 'JSON':
        import json
        stitched['metadata'] = fxn.get_metadata_dict(input_args['species'], 'stitchr')
        stitched['Warnings/Errors'] = ''.join([str(stitch_warnings[x].message)
                                                               for x in range(len(stitch_warnings))])

        with open(fxn.get_output_name(stitched) + '.json', 'w') as out_file:
            json.dump(stitched, out_file)

    # If a known/partial amino acid sequence provided, ensure they match up with a quick printed alignment
    if input_args['aa']:
        from Bio import pairwise2
        from Bio.pairwise2 import format_alignment
        alignments = pairwise2.align.globalxx(input_args['aa'],
                                              fxn.translate_nt('N' * stitched['translation_offset']
                                                               + stitched['stitched_nt']))
        for i in range(0, 600, 60):
            print('\n')
            if i > len(alignments[0][0]):
                break
            for y in [x[i:i+60] for x in format_alignment(*alignments[0]).split('\n')[:3]]:
                print(y)

    # Finally print any warnings generated during the running
    if stitch_warnings:
        print('\n'.join([str(stitch_warnings[x].message) for x in range(len(stitch_warnings))]))
