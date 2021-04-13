#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
thimble.py

Runs stiTChR in a high-throughput manner, on a tab separated input file.

No fancy backronym or contrived silent capitals - it just helps you stitch faster.

"""


import functions as fxn
import stitchr as st
import warnings
import argparse
import sys
import os

__version__ = '0.4.1'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'

sys.tracebacklimit = 0


def args():
    """
    args(): Obtains command line arguments which dictate the script's behaviour
    """

    # Help flag
    parser = argparse.ArgumentParser(description="thimble v" + str(__version__) + '\n' +
                                                 ": Run stiTChR on multiple and paired TCRs")

    # Input and output options
    parser.add_argument('-in', '--in_file', required=True, type=str,
                help='Tab-separated input file. Must fit the format of the example file. Required.')

    parser.add_argument('-o', '--out_file', required=True, type=str,
                help='Path/file name for output tsv file. Required.')

    parser.add_argument('-s', '--species', required=False, type=str, default='human',
                help='Species. Optional. Only available default options are \'human\' or \'mouse\'. Default = human.')

    parser.add_argument('-cu', '--codon_usage', required=False, type=str,
                help='Path to a file of Kazusa-formatted codon usage frequencies. Optional.')

    parser.add_argument('-xg', '--extra_genes', action='store_true', required=False, default=False,
       help='Optional flag to use additional (non-database/-natural) sequences in the \'additional-genes.fasta\' file.')

    parser.add_argument('-jt', '--j_warning_threshold', required=False, type=int, default=3,
                help='J gene substring length warning threshold. Default = 3. '
                     'Decrease to get fewer notes on short J matches. Optional.')

    return parser.parse_args()


in_headers = ['TCR_name', 'TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3',
               'TRAC', 'TRBC', 'TRA_leader', 'TRB_leader', 'Linker', 'Link_order',
              'TRA_5_prime_seq', 'TRA_3_prime_seq', 'TRB_5_prime_seq', 'TRB_3_prime_seq']

convert_fields = {'TRAV': 'v', 'TRAJ': 'j', 'TRA_CDR3': 'cdr3', 'TRBV': 'v', 'TRBJ': 'j', 'TRB_CDR3': 'cdr3',
                  'TRAC': 'c', 'TRBC': 'c', 'TRA_leader': 'l', 'TRB_leader': 'l', 'TCR_name': 'name',
                  'TRA_5_prime_seq': '5_prime_seq', 'TRA_3_prime_seq': '3_prime_seq',
                 'TRB_5_prime_seq': '5_prime_seq', 'TRB_3_prime_seq': '3_prime_seq'}

stitch_list_fields = ['_name', 'V', 'J', 'C', '_CDR3', '_leader', '_5_prime_seq', '_3_prime_seq']

out_headers = ['TCR_name', 'TRA_nt', 'TRB_nt', 'TRA_aa', 'TRB_aa',
               'TRAV', 'TRAJ', 'TRA_CDR3', 'TRBV', 'TRBJ', 'TRB_CDR3',
               'TRAC', 'TRBC', 'TRA_leader', 'TRB_leader', 'Linker', 'Link_order',
              'TRA_5_prime_seq', 'TRA_3_prime_seq', 'TRB_5_prime_seq', 'TRB_3_prime_seq',
               'Linked_nt', 'Linked_aa', 'Warnings/Errors']

if __name__ == '__main__':

    # Get input arguments, get the required data read in
    fxn.check_scripts_dir()
    input_args = vars(args())

    codons = fxn.get_optimal_codons(input_args['codon_usage'], input_args['species'].upper())
    linker_dict = fxn.get_linker_dict()

    tcr_dat = {}
    tcr_functionality = {}

    for c in ['TRA', 'TRB']:

        tmp_tcr_dat, tmp_functionality = fxn.get_imgt_data(c, st.gene_types, input_args['species'].upper())
        tcr_dat[c] = tmp_tcr_dat
        tcr_functionality[c] = tmp_functionality

        if 'extra_genes' in input_args:
            if input_args['extra_genes']:
                tcr_dat[c], tcr_functionality[c] = fxn.get_additional_genes(tcr_dat[c], tcr_functionality[c])
                input_args['skip_c_checks'] = True
            else:
                input_args['skip_c_checks'] = False

    # Then go through in file and stitch each TCR on each line
    if not os.path.isfile(input_args['in_file']):
        raise IOError(input_args['in_file'] + " not detected - please check and specify in file again.")

    with fxn.opener(input_args['in_file']) as in_file:

        line_count = 0
        out_data = ['\t'.join(out_headers)]

        for line in in_file:

            bits = line.replace('\n', '').replace('\"', '').split('\t')

            # Check there's the right number of columns (with the right headers)
            if line_count == 0:
                if bits != in_headers:
                    raise ValueError("Headers in input file don't match the expectations - please check template.")

            else:
                # Generate a dict per line ...
                sorted_row_bits = {x: '' for x in out_headers}

                if len(bits) != len(in_headers):  # Pad entries from tsvs where whitespace on right is missing
                    bits += [''] * (len(in_headers) - len(bits))

                for i in range(len(in_headers)):
                    sorted_row_bits[in_headers[i]] = bits[i]

                # ... along with chain-specific dicts which are then folded into
                for c in ['TRA', 'TRB']:

                    with warnings.catch_warnings(record=True) as chain_warnings:
                        warnings.simplefilter("always")

                        tcr_bits = {'skip_c_checks': input_args['skip_c_checks']}

                        # Convert naming to output formats
                        for field in [x for x in sorted_row_bits if c in x]:
                            if sorted_row_bits[field]:
                                tcr_bits[convert_fields[field]] = sorted_row_bits[field]

                        # Determine whether or not there's supposed to be a TCR for this chain
                        if len(tcr_bits) > 1:

                            # At the very least we need the relevant TCR info (V/J/CDR3)
                            if not all(section in tcr_bits for section in ['v', 'j', 'cdr3']):
                                warnings.warn("Incomplete TCR information - need V/J/CDR3 as minimum.")

                            else:
                                tcr_bits = fxn.autofill_input(tcr_bits, c)
                                tcr_bits = fxn.tweak_thimble_input(tcr_bits, input_args)

                                try:
                                    out_list, stitched, offset = st.stitch(tcr_bits, c, tcr_dat[c],
                                                                   tcr_functionality[c], codons,
                                                                   input_args['j_warning_threshold'])
                                    sorted_row_bits[c + '_nt'] = stitched
                                    sorted_row_bits[c + '_aa'] = fxn.translate_nt('N' * offset + stitched)
                                    sorted_row_bits.update(dict(list(zip([c + x for x in stitch_list_fields],
                                                                         out_list))))

                                except Exception as message:
                                    sorted_row_bits['Warnings/Errors'] += '(' + c + ') ' + str(message)
                                    sorted_row_bits['Warnings/Errors'] += 'Cannot stitch a sequence for ' + c + '. '

                        # Store all chain related warning messages too, in same field, ignoring irrelevant errors
                        # (ignoring Biopython's len%3 != 0 error, which isn't very helpful given the circumstances)
                        sorted_row_bits['Warnings/Errors'] += ' '.join(
                            ['(' + c + ') ' + str(chain_warnings[x].message) for x in range(len(chain_warnings))
                             if 'Biopython' not in str(chain_warnings[x].category) and
                             'DeprecationWarning' not in str(chain_warnings[x].category)])

                with warnings.catch_warnings(record=True) as link_warnings:
                    warnings.simplefilter("always")

                    # If sequences are to be linked, determine the appropriate linker sequence and join together
                    if sorted_row_bits['Linker']:

                        if sorted_row_bits['Link_order']:
                            sorted_row_bits['Link_order'] = sorted_row_bits['Link_order'].upper()

                            # Only allow valid orders
                            if sorted_row_bits['Link_order'] not in ['BA', 'AB']:
                                warnings.warn("Error: given link order not valid (not AB or BA) - defaulting to BA.")
                                sorted_row_bits['Link_order'] = 'BA'

                        else:
                            # Default option is B
                            sorted_row_bits['Link_order'] = 'BA'

                        tr1, tr2 = sorted_row_bits['Link_order'][0], sorted_row_bits['Link_order'][1]

                        if sorted_row_bits['TRA_nt'] and sorted_row_bits['TRB_nt']:
                            try:
                                linker_seq = fxn.get_linker_seq(sorted_row_bits['Linker'], linker_dict)

                                linked = sorted_row_bits['TR' + tr1 + '_nt'] + \
                                         linker_seq + sorted_row_bits['TR' + tr2 + '_nt']
                                sorted_row_bits['Linked_nt'] = linked
                                sorted_row_bits['Linked_aa'] = fxn.translate_nt('N' * offset + linked)

                                # Add warnings if sequences applied at maybe the wrong ends of things
                                if (sorted_row_bits['TRA_5_prime_seq'] and ord == 'BA') or \
                                    (sorted_row_bits['TRB_5_prime_seq'] and ord == 'AB'):
                                    warnings.warn("Warning: " + ord + " order specified, but "
                                                                "3\' chain has an additional 5\' sequence provided. ")

                                if (sorted_row_bits['TRA_3_prime_seq'] and ord == 'AB') or \
                                    (sorted_row_bits['TRB_3_prime_seq'] and ord == 'BA'):
                                    warnings.warn("Warning: " + ord + " order specified, but "
                                                                "5\' chain has an additional 3\' sequence provided. ")

                            # Store any error messages
                            except Exception as message:
                                sorted_row_bits['Warnings/Errors'] += str(message)


                        else:
                            warnings.warn("Error: need both a TRA and TRB to link. ")

                    # And add any warnings/errors derived from the linkage
                    sorted_row_bits['Warnings/Errors'] += ' '.join(
                        ['(Link) ' + str(link_warnings[x].message) for x in range(len(link_warnings))
                         if 'Biopython' not in str(link_warnings[x].category) and
                         'DeprecationWarning' not in str(link_warnings[x].category)])

                # Store output as one long string, for a single write-out once input file is finished
                out_data.append('\t'.join([sorted_row_bits[x] for x in out_headers]))

            line_count += 1

    # Write out data
    if not input_args['out_file'].lower().endswith('.tsv'):
        input_args['out_file'] += '.tsv'

    with open(input_args['out_file'], 'w') as out_file:
        out_string = '\n'.join(out_data)
        out_file.write(out_string)
