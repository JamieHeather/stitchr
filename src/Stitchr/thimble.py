# -*- coding: utf-8 -*-

"""
thimble.py

Runs stiTChR in a high-throughput manner, on a tab separated input file.

No fancy backronym or contrived silent capitals - it just helps you stitch faster.

"""


from . import stitchrfunctions as fxn
from . import stitchr as st
from . import __version__
import warnings
import argparse
import os
import sys
import collections as coll
from time import time


__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

sys.tracebacklimit = 0  # comment when debugging


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(description="thimble:\nRun stiTChR on multiple and paired TCRs")

    # Input and output options
    parser.add_argument('-in', '--in_file', required=True, type=str,
                        help="Tab-separated input file. Must fit the format of the example file. Required.")

    parser.add_argument('-o', '--out_file', required=True, type=str,
                        help="Path/file name for output tsv file. Required.")

    parser.add_argument('-s', '--species', required=False, type=str, default='HUMAN',
                        help="Species. Optional. Default = 'HUMAN', options depend on directories in Data directory.")

    parser.add_argument('-m', '--mode', required=False, type=str, default='TSV',
                        help="Standard out output mode. Options are 'TSV' (default), 'GB', or 'JSON'.")

    parser.add_argument('-ol', '--only_linked', required=False, action='store_true', default=False,
                        help="Optional flag to speed up thimble by only translating linked/paired chain sequences.")

    parser.add_argument('-xg', '--extra_genes', action='store_true', required=False, default=False,
                        help="Optional flag to use additional (non-database/-natural) sequences in "
                             "the \'additional-genes.fasta\' file.")

    parser.add_argument('-sl', '--seamless', action='store_true', required=False, default=False,
                        help="Optional flag to integrate known nucleotide sequences seamlessly. \nNB: "
                             "nucleotide sequences covering the CDR3 junction with additional V gene context required.")

    parser.add_argument('-r', '--receptor', required=False, type=str,
                        help="Receptor to stitch, alpha/beta or gamma/delta TCR. Optional. Give double or single digits"
                             " (e.g. 'ab', 'gd', 'd', 'b').")

    parser.add_argument('-p', '--preferred_alleles_path', required=False, type=str, default='',
                        help="Path to a file of preferred alleles to use when no allele specified (instead of *01). "
                             "Optional.")

    parser.add_argument('-cu', '--codon_usage', required=False, type=str,
                        help="Path to a file of Kazusa-formatted codon usage frequencies. Optional.")

    parser.add_argument('-jt', '--j_warning_threshold', required=False, type=int, default=3,
                        help="J gene substring length warning threshold. Default = 3. "
                             "Decrease to get fewer notes on short J matches. Optional.")

    parser.add_argument('-sc', '--skip_c_checks', action='store_true', required=False, default=False,
                        help="Optional flag to skip usual constant region gene checks.")

    parser.add_argument('-sn', '--skip_n_checks', action='store_true', required=False, default=False,
                        help="Optional flag to skip the usual CDR3 N terminal (conserved 2nd Cys) checks.")

    parser.add_argument('-nl', '--no_leader', action='store_true', required=False, default=False,
                        help="Optional flag to ignore leader sequences (not recommended for paired/linked sequences!).")

    parser.add_argument('--version', action='version', version=__version__,
                        help="Print current stitchr package version.")

    parser.add_argument('--cite', action=fxn.GetCitation, help="Print citation details.", nargs=0)

    return parser.parse_args()


def locus_to_trx(trx_str):
    """
    :param trx_str: str containing (a) reference to a specific locus/loci (e.g. TRA/TRB or TRG/TRD)
    :return: same string, but with those loci specific IDs replaced with generic TR1/TR2 strings
    """
    return trx_str.replace('TRA', 'TR1').replace('TRB', 'TR2').replace('TRG', 'TR1').replace('TRD', 'TR2')


def populate_blanks(chain_dict, necessary_fields):
    """
    :param chain_dict: dictionary containing the stitchr fields for a given chain's recombination
    :param necessary_fields: list of fields that the tcr dictionary requires to pass the downstream steps
    :return: chain_dict with blank entries for each of the necessary_fields that it was missing
    """
    for nf in necessary_fields:
        if nf not in chain_dict:
            chain_dict[nf] = ''

    return chain_dict


in_headers = {'TRA/TRB': ['TCR_name', 'TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3',
                          'TRAC', 'TRBC', 'TRA_leader', 'TRB_leader', 'Linker', 'Link_order',
                          'TRA_5_prime_seq', 'TRA_3_prime_seq', 'TRB_5_prime_seq', 'TRB_3_prime_seq'],
              'TRG/TRD': ['TCR_name', 'TRGV', 'TRGJ', 'TRG_CDR3', 'TRDV', 'TRDJ', 'TRD_CDR3',
                          'TRGC', 'TRDC', 'TRG_leader', 'TRD_leader', 'Linker', 'Link_order',
                          'TRG_5_prime_seq', 'TRG_3_prime_seq', 'TRD_5_prime_seq', 'TRD_3_prime_seq']}

convert_fields = {'TRAV': 'v', 'TRAJ': 'j', 'TRA_CDR3': 'cdr3', 'TRBV': 'v', 'TRBJ': 'j', 'TRB_CDR3': 'cdr3',
                  'TRAC': 'c', 'TRBC': 'c', 'TRA_leader': 'l', 'TRB_leader': 'l', 'TCR_name': 'name',
                  'TRA_5_prime_seq': '5_prime_seq', 'TRA_3_prime_seq': '3_prime_seq',
                  'TRB_5_prime_seq': '5_prime_seq', 'TRB_3_prime_seq': '3_prime_seq',
                  'TRGV': 'v', 'TRGJ': 'j', 'TRG_CDR3': 'cdr3', 'TRDV': 'v', 'TRDJ': 'j', 'TRD_CDR3': 'cdr3',
                  'TRGC': 'c', 'TRDC': 'c', 'TRG_leader': 'l', 'TRD_leader': 'l',
                  'TRG_5_prime_seq': '5_prime_seq', 'TRG_3_prime_seq': '3_prime_seq',
                  'TRD_5_prime_seq': '5_prime_seq', 'TRD_3_prime_seq': '3_prime_seq'}

pre_stitch_list_fields = ['name', 'v', 'j', 'c', 'cdr3', 'l', '5_prime_seq', '3_prime_seq']
post_stitch_list_fields = ['_name', 'V', 'J', 'C', '_CDR3', '_leader', '_5_prime_seq', '_3_prime_seq']

out_headers = {'TRA/TRB': ['TCR_name', 'TRA_nt', 'TRB_nt', 'TRA_aa', 'TRB_aa',
                           'TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3',
                           'TRAC', 'TRBC', 'TRA_leader', 'TRB_leader', 'Linker', 'Link_order',
                           'TRA_5_prime_seq', 'TRA_3_prime_seq', 'TRB_5_prime_seq', 'TRB_3_prime_seq',
                           'Linked_nt', 'Linked_aa', 'Warnings/Errors'],
               'TRG/TRD': ['TCR_name', 'TRG_nt', 'TRD_nt', 'TRG_aa', 'TRD_aa',
                           'TRGV', 'TRGJ', 'TRG_CDR3', 'TRDV', 'TRDJ', 'TRD_CDR3',
                           'TRGC', 'TRDC', 'TRG_leader', 'TRD_leader', 'Linker', 'Link_order',
                           'TRG_5_prime_seq', 'TRG_3_prime_seq', 'TRD_5_prime_seq', 'TRD_3_prime_seq',
                           'Linked_nt', 'Linked_aa', 'Warnings/Errors']}


def main():

    # Get input arguments, get the required data read in
    input_args = vars(args())
    start = time()

    # Get species, in order of command line arg > inferred from input file name > default human
    if input_args['species']:
        if input_args['species'].upper() in fxn.find_species_covered():
            species = input_args['species'].upper()
        else:
            raise IOError("No data available for requested species: " + input_args['species'])
    else:
        species_inference = fxn.infer_species(input_args['in_file'])
        if species_inference:
            species = species_inference
        else:
            species = 'HUMAN'

    codons = fxn.get_optimal_codons(input_args['codon_usage'], species)
    linker_dict = fxn.get_linker_dict()

    # Fetch required motifs from the species folder in the data directory
    j_res, low_conf_js = fxn.get_j_motifs(species)
    c_res = fxn.get_c_motifs(species)

    tcr_dat = {}
    tcr_functionality = {}
    preferences = {}

    # Figure out whether a/b or g/d
    if input_args['receptor']:
        input_receptor = input_args['receptor'].upper()
        if ('A' in input_receptor or 'B' in input_receptor) and not ('G' in input_receptor or 'D' in input_receptor):
            receptor = 'TRA/TRB'
        elif ('G' in input_receptor or 'D' in input_receptor) and not ('A' in input_receptor or 'B' in input_receptor):
            receptor = 'TRG/TRD'
        else:
            raise IOError("Unable to determine receptor from '-r' command string: " + input_receptor + ". ")

    else:  # If not explicitly provided, infer from input TSV headers
        with fxn.opener(input_args['in_file']) as in_file:
            for line in in_file:
                if 'TRAV' in line and 'TRGV' not in line:
                    receptor = 'TRA/TRB'
                elif 'TRGV' in line and 'TRAV' not in line:
                    receptor = 'TRG/TRD'
                else:
                    raise IOError("Unable to determine receptor from input file header, please check template. ")
                break

    # Define the individual receptors (chains or loci, i.e. TRA and TRB or TRG and TRD) in play
    r1, r2 = receptor.split('/')
    for c in [r1, r2]:

        tmp_tcr_dat, tmp_functionality, partial = fxn.get_ref_data(c, st.gene_types, species)
        tcr_dat[c] = tmp_tcr_dat
        tcr_functionality[c] = tmp_functionality

        if 'extra_genes' in input_args:
            if input_args['extra_genes']:
                tcr_dat[c], tcr_functionality[c] = fxn.get_additional_genes(tcr_dat[c], tcr_functionality[c])
                input_args['skip_c_checks'] = True
            else:
                input_args['skip_c_checks'] = False

        # Allow for provision of preferred alleles
        if input_args['preferred_alleles_path']:
            preferences[c] = fxn.get_preferred_alleles(input_args['preferred_alleles_path'], list(fxn.regions.values()),
                                                       tcr_dat[c], partial, c)
        else:
            preferences[c] = ''

    if not input_args['out_file'].lower().endswith('.tsv'):
        input_args['out_file'] += '.tsv'

    # If using certain file output options, sort out the environment to cope with their specific outputs
    input_args['mode'] = input_args['mode'].upper()
    if input_args['mode'] in ['GB', 'JSON']:

        base_description = fxn.get_metadata_text(input_args['species'], 'thimble', __version__)

        base_description += ("Produced using the thimble command: '" + ' '.join(sys.argv[:]) + "'. ")

        if input_args['extra_genes']:
            base_description += ("Note that this was run using the 'extra genes' flag, and thus additional "
                                 "non-reference sequences may have been incorporated. ")

        input_args['out_folder'] = input_args['out_file'][:-4]
        if not os.path.exists(input_args['out_folder']):
            os.makedirs(input_args['out_folder'])

        if input_args['mode'] == 'GB':
            # Also create the individual subdirectories (unless 'only_linked' option is selected)
            subdirs = ['linked']
            if not input_args['only_linked']:
                subdirs += [r1, r2]
            for subdir in subdirs:
                subdir_to_make = os.path.join(input_args['out_folder'], subdir)
                if not os.path.exists(subdir_to_make):
                    os.makedirs(subdir_to_make)

        elif input_args['mode'] == 'JSON':
            import json
    else:
        input_args['out_folder'] = os.getcwd()

    # Then go through in file and stitch each TCR on each line
    if not os.path.isfile(input_args['in_file']):
        raise IOError(input_args['in_file'] + " not detected - please check and specify in file again.")

    with fxn.opener(input_args['in_file']) as in_file:

        line_count = 1
        tcr_count = 1
        out_data = ['\t'.join(out_headers[receptor])]

        for line in in_file:

            bits = line.replace('\n', '').replace('\r', '').replace('\"', '').split('\t')

            # Check there's the right number of columns (with the right headers)
            if line_count == 1:
                if bits != in_headers[receptor]:
                    raise ValueError("Headers in input file don't match the expectations - please check template.")

            else:
                # Generate a dict per line ...
                line_sorted_row_bits = {x: '' for x in out_headers[receptor]}
                entry_line_warnings = []

                if len(bits) != len(in_headers[receptor]):  # Pad entries from TSVs where whitespace on right is missing
                    bits += [''] * (len(in_headers[receptor]) - len(bits))

                # Allow for users to specify multiple options per line, to allow the quick stitching of related TCRs
                multiple_options = []
                for i in range(len(in_headers[receptor])):
                    line_sorted_row_bits[in_headers[receptor][i]] = bits[i]
                    if ',' in bits[i] or '%' in bits[i]:
                        multiple_options.append(i)

                # Make sure there's a (usable) name
                if line_sorted_row_bits['TCR_name']:
                    line_sorted_row_bits['TCR_name'] = line_sorted_row_bits['TCR_name'].strip()
                else:
                    line_sorted_row_bits['TCR_name'] = 'TCR-' + "{:05d}".format(tcr_count)

                tcr_lines = [line_sorted_row_bits]
                if multiple_options:
                    # Iterate over all possible combinations for multiply defined options
                    for multi_var in multiple_options:
                        updated_lines = []
                        var_name = in_headers[receptor][multi_var]
                        region = fxn.regions[var_name[3:].replace('_', '')[0].lower()]

                        # A manually entered list
                        if ',' in bits[multi_var]:
                            var_options = bits[multi_var].split(',')

                        # If not a list, it should be a wild card, but only for specified fields
                        elif var_name[3:] not in ['V', 'J', 'C', '_leader']:
                            entry_line_warnings.append("% wild cards only allowed in V/J/C/leader fields. ")
                            break

                        # Wild card just all alleles for a certain gene
                        elif '*%' in bits[multi_var]:
                            gene = bits[multi_var].split('*')[0]
                            if gene in tcr_dat[var_name[:3]][region]:
                                var_options = [gene + '*' + x for x in tcr_dat[var_name[:3]][region][gene].keys()]
                            else:
                                entry_line_warnings.append("Specified multiple entry option gene '" + gene + "' is not "
                                                           "found in the " + region + " region data. ")
                                break

                        # Or wild card all genes/allele combinations for a given region (AKA throw it all at the wall)
                        elif bits[multi_var] == '%':
                            var_options = []
                            for g in tcr_dat[var_name[:3]][region]:
                                for a in tcr_dat[var_name[:3]][region][g]:
                                    var_options.append(g + '*' + a)

                        else:
                            entry_line_warnings.append("Unexpected multi-option field entry detected: " +
                                                       bits[multi_var] + '. ')
                            break

                        # Then loop through those options filling out all specified combinations to stitch
                        for vo in var_options:
                            for out_line in tcr_lines:
                                updated_lines.append(dict(out_line))
                                updated_lines[-1][var_name] = vo
                                updated_lines[-1]['TCR_name'] = updated_lines[-1]['TCR_name'] + '-' + \
                                                                str(var_options.index(vo))
                        tcr_lines = updated_lines

                if entry_line_warnings:
                    tcr_lines[0]['Warnings/Errors'] = ''.join(entry_line_warnings)
                    continue

                for sorted_row_bits in tcr_lines:

                    thimble_dict = coll.defaultdict(str) # Nested dictionary for tracking all TCR-associated data
                    thimble_dict['linked'] = {'Warnings/Errors': '', 'linkable': True, 'Linker': ''}
                    thimble_dict['out'] = sorted_row_bits

                    # Need to run the stitch function on each individual submitted chain
                    for c in [r1, r2]:
                        chain_notes = ''

                        with warnings.catch_warnings(record=True) as chain_warnings:
                            warnings.simplefilter("always")

                            # Generate the necessary fields for stitching
                            tcr_bits = {'skip_c_checks': input_args['skip_c_checks'],
                                        'skip_n_checks': input_args['skip_n_checks'],
                                        'no_leader': input_args['no_leader'],
                                        'species': species,
                                        'mode': input_args['mode']}

                            if 'seamless' in input_args:
                                if input_args['seamless']:
                                    tcr_bits['seamless'] = True

                            # # Convert naming to output formats
                            for field in [x for x in sorted_row_bits if c in x]:
                                if field in convert_fields:
                                    if sorted_row_bits[field]:
                                        tcr_bits[convert_fields[field]] = sorted_row_bits[field]

                            # Determine whether there's supposed to be a TCR for this chain
                            # At the very least we need the relevant TCR info (V/J/CDR3)
                            featured_bits = [x for x in tcr_bits if x in ['v', 'j', 'cdr3']]

                            # If no TCR features present, just skip
                            if len(featured_bits) == 0:
                                thimble_dict['linked']['linkable'] = False
                                continue

                            elif len(featured_bits) in [1, 2]:
                                warnings.warn("Incomplete TCR information - need V/J/CDR3 as minimum.")
                                thimble_dict['linked']['linkable'] = False

                            else:
                                tcr_bits = fxn.autofill_input(populate_blanks(tcr_bits, pre_stitch_list_fields), c)
                                tcr_bits = fxn.tweak_thimble_input(tcr_bits)

                                try:
                                    # Actually run the stitching of this chain
                                    thimble_dict[c] = st.stitch(tcr_bits, tcr_dat[c], tcr_functionality[c], partial,
                                                                codons, input_args['j_warning_threshold'],
                                                                preferences[c], c_res, j_res, low_conf_js)

                                    thimble_dict['out'][c + '_nt'] = thimble_dict[c]['stitched_nt']

                                    # Skip translation of single chains if the 'only_linked' option is selected
                                    if not input_args['only_linked']:
                                        thimble_dict[c]['stitched_aa'] = fxn.pad_trans(thimble_dict[c]['stitched_nt'],
                                                                                       thimble_dict[c]['translation_offset'])
                                        thimble_dict['out'][c + '_aa'] = thimble_dict[c]['stitched_aa']
                                    else:
                                        thimble_dict[c]['stitched_aa'], thimble_dict['out'][c + '_aa'] = '', ''

                                    thimble_dict['out'].update(dict(list(zip(
                                        [c + x for x in post_stitch_list_fields], thimble_dict[c]['out_list']))))

                                except Exception as message:
                                    chain_notes += str(message) + 'Cannot stitch a sequence for ' + c + '. '
                                    thimble_dict['linked']['linkable'] = False

                            # Store all chain related warning messages (ignoring irrelevant errors)
                            # NB only save chain specific data if there's a rearranged chain (for single chain outputs)
                            # but all chain data will be stored in the final all encompassing entry
                            chain_notes += ' '.join(
                                [str(chain_warnings[x].message) for x in range(len(chain_warnings))
                                 if 'DeprecationWarning' not in str(chain_warnings[x].category)])
                            if c in thimble_dict.keys():
                                thimble_dict[c]['Warnings/Errors'] = chain_notes

                        if chain_notes:
                            thimble_dict['linked']['Warnings/Errors'] += chain_notes.replace('.', ' (' + c + ').')

                    with warnings.catch_warnings(record=True) as link_warnings:
                        warnings.simplefilter("always")

                        # If sequences are to be linked, determine the appropriate linker sequence and join together
                        if sorted_row_bits['Linker'] and thimble_dict['linked']['linkable'] == True:

                            if sorted_row_bits['Link_order']:
                                thimble_dict['linked']['Link_order'] = sorted_row_bits['Link_order'].upper()

                                # Only allow valid orders
                                if sorted_row_bits['Link_order'] not in [r2[2]+r1[2], r1[2]+r2[2]]:
                                    warnings.warn("Error: given link order not valid (not " + r1[2]+r2[2] + " or " +
                                                  r2[2]+r1[2] + ") - defaulting to " + r2[2]+r1[2] + ".")
                                    thimble_dict['linked']['Link_order'] = r2[2]+r1[2]

                            else:
                                # Default option is B
                                thimble_dict['linked']['Link_order'] = r2[2]+r1[2]

                            # Extract just the defining locus character (e.g. A/B) to infer proper link order
                            tr1, tr2 = thimble_dict['linked']['Link_order'][0], thimble_dict['linked']['Link_order'][1]

                            if thimble_dict[r1]['stitched_nt'] and thimble_dict[r2]['stitched_nt']:
                                try:
                                    thimble_dict['linked']['Linker_id'] = sorted_row_bits['Linker']
                                    thimble_dict['linked']['linker_nt'] = fxn.get_linker_seq(sorted_row_bits['Linker'],
                                                                                             linker_dict)

                                    thimble_dict['linked']['Linked_nt'] = thimble_dict['TR' + tr1]['stitched_nt'] + \
                                                                          thimble_dict['linked']['linker_nt'] + \
                                                                          thimble_dict['TR' + tr2]['stitched_nt']
                                    thimble_dict['linked']['Linked_aa'] = fxn.pad_trans(thimble_dict['linked']['Linked_nt'],
                                                                                        thimble_dict[r1]
                                                                                        ['translation_offset'])

                                    # Add warnings if sequences applied at maybe the wrong ends of things
                                    if (thimble_dict[r1]['in']['5_prime_seq']
                                        and thimble_dict['linked']['Link_order'] == r2[2]+r1[2]) or \
                                            (thimble_dict[r2]['in']['5_prime_seq']
                                             and thimble_dict['linked']['Link_order'] == r1[2]+r2[2]):
                                        warnings.warn("Warning: " + thimble_dict['linked']['Link_order'] + " order "
                                                    "specified, but 3' chain has an additional 5' sequence provided. ")

                                    if (thimble_dict[r1]['in']['3_prime_seq']
                                        and thimble_dict['linked']['Link_order'] == r1[2]+r2[2]) or \
                                            (thimble_dict[r2]['in']['3_prime_seq']
                                             and thimble_dict['linked']['Link_order'] == r2[2]+r1[2]):
                                        warnings.warn("Warning: " + thimble_dict['linked']['Link_order'] + " order "
                                                    "specified, but 5' chain has an additional 3' sequence provided. ")

                                # Store any error messages
                                except Exception as message:
                                    thimble_dict['linked']['Warnings/Errors'] += str(message)

                            else:
                                warnings.warn("Error: need both a " + r1 + " and " + r2 + " to link. ")

                        # Add any warnings/errors derived from the linkage, or either constituent chain
                        if thimble_dict['linked']['Warnings/Errors'] or link_warnings:
                            thimble_dict['linked']['Warnings/Errors'] += ' '.join(
                                [str(link_warnings[x].message) for x in range(len(link_warnings))
                                 if 'DeprecationWarning' not in str(link_warnings[x].category)])

                        if not thimble_dict['linked']['Warnings/Errors']:
                            thimble_dict['linked']['Warnings/Errors'] = "[None]"

                        # Move over the remaining fields from the 'linked' dict to the 'out' dict
                        for f in ['Linked_nt', 'Linked_aa']:
                            if f in thimble_dict['linked']:
                                thimble_dict['out'][f] = thimble_dict['linked'][f]
                            else:
                                thimble_dict['out'][f] = ''

                    thimble_dict['out']['Warnings/Errors'] = thimble_dict['linked']['Warnings/Errors']

                    # Store output as one long string, for a single write-out once input file is finished
                    for field in [x for x in out_headers[receptor] if x not in thimble_dict['out']]:
                        thimble_dict['out'][field] = thimble_dict['linked'][field]

                    out_data.append('\t'.join([thimble_dict['out'][x] for x in out_headers[receptor]]))

                    # If requested, also output additional more detailed file types
                    if input_args['mode'] == 'JSON':
                        thimble_dict['metadata'] = fxn.get_metadata_dict(input_args['species'], 'thimble')
                        with open(os.path.join(input_args['out_folder'], thimble_dict['out']['TCR_name'])
                                  + '.json', 'w') as out_file:
                            json.dump(thimble_dict, out_file)

                    elif input_args['mode'] == 'GB':

                        # Generate files for each chain, unless linked only option selected
                        subdirs = ['linked']
                        if not input_args['only_linked']:
                            subdirs += [r1, r2]
                        for subdir in subdirs:
                            if subdir not in thimble_dict:
                                continue
                            spec_description = base_description
                            if thimble_dict[subdir]['Warnings/Errors'] != "[None]":
                                spec_description +=  ("\n" + fxn.spacer + "\n" + fxn.spacer +
                                                      "This rearrangement produced the following warning messages:"
                                                      "\n" + fxn.spacer + thimble_dict[subdir]['Warnings/Errors'])

                            # TODO move stitchr/thimble-to-genbank formatting to own function?
                            gb_feats = []
                            if subdir != 'linked':
                                gb_feats = []
                                for feat in thimble_dict[subdir]['indexes']:
                                    gb_feats.append((thimble_dict[subdir]['used'][feat], '',
                                                     [thimble_dict[subdir]['seqs'][feat]],
                                                     [thimble_dict[subdir]['indexes'][feat]]))
                            else:
                                if not thimble_dict[subdir]['linkable']:
                                    continue
                                # Establish the order and re-engineer up the indexes to ensure correct linked placement
                                thimble_dict[subdir]['stitched_nt'] = thimble_dict[subdir]['Linked_nt']

                                sorted_chain1 = sorted(thimble_dict['TR' + tr1]['indexes'].items(), key=lambda x: x[1])
                                sorted_chain2 = sorted(thimble_dict['TR' + tr2]['indexes'].items(), key=lambda x: x[1])

                                for feat in sorted_chain1:
                                    gb_feats.append((thimble_dict['TR' + tr1]['used'][feat[0]], '',
                                                    [thimble_dict['TR' + tr1]['seqs'][feat[0]]],
                                                    [thimble_dict['TR' + tr1]['indexes'][feat[0]]]))

                                offset = len(thimble_dict['TR' + tr1]['stitched_nt'])
                                gb_feats.append((thimble_dict['linked']['Linker_id'], '',
                                                    [thimble_dict['linked']['linker_nt']],
                                                    [offset]))
                                offset += len(thimble_dict['linked']['linker_nt'])

                                for feat in sorted_chain2:
                                    gb_feats.append((thimble_dict['TR' + tr2]['used'][feat[0]], '',
                                                    [thimble_dict['TR' + tr2]['seqs'][feat[0]]],
                                                    [thimble_dict['TR' + tr2]['indexes'][feat[0]] + offset]))

                            genbank_params = {'sequence_name': thimble_dict['out']['TCR_name'],
                                              'full_sequence': thimble_dict[subdir]['stitched_nt'],
                                              'description': spec_description, 'topology': 'linear',
                                              'features': gb_feats,
                                              'save_dir_path': os.path.join(input_args['out_folder'], subdir),
                                              'species': input_args['species'],
                                              'numbering': False,
                                              'plot_multi': True, 'division': 'SYN', 'journal': 'thimble'}

                            fxn.output_genbank(**genbank_params)

                tcr_count += 1

            line_count += 1

    # TODO delete unused output folders for GB/JSON output?

    time_taken = time() - start
    print("Took", str(round(time_taken, 2)), "seconds")

    # Write out data
    with open(os.path.join(input_args['out_folder'], input_args['out_file']), 'w') as out_file:
        out_string = '\n'.join(out_data)
        out_file.write(out_string)
