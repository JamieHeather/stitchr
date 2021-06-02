#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
gui-stitchr.py

A graphical user interface for stitchr, powered by PySimpleGUI

"""


import PySimpleGUI as sg
import os
import functions as fxn
import stitchr as st
import thimble as th
import collections as coll
import warnings

__version__ = '0.1.0'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def check_species(value_data):
    """
    :param value_data: the 'values' dict, detailing the contents of the input fields
    :return: the inferred species, HUMAN or MOUSE
    """
    if value_data['rad_hs']:
        return 'HUMAN'
    elif value_data['rad_mm']:
        return 'MOUSE'
    else:
        raise IOError('Species cannot be determined')


def check_is_file(path_to_file):
    """
    :param path_to_file:
    """
    if not os.path.isfile(path_to_file):
        raise IOError("Cannot find file: " + path_to_file)


extra_gene_text = "(FASTA format, optional)"
box_width = 50
sz = (box_width, 1)
half_sz = (int(box_width / 2 - 2), 1)          # For half a column
third_sz = (int(box_width / 3 - 1.3), 1)          # For one third of a column

linkers = fxn.get_linker_dict()
link_choices = list(linkers.keys()) + ['Custom']

# General interface column
col1 = [

    [sg.Button('Example data', size=half_sz), sg.Button('Reset form', size=half_sz)],

    [sg.FileBrowse(key="uploaded-tcr", size=half_sz, button_text='Find TCR input file'),
     sg.Button("Upload TCR details", size=half_sz)],

    [sg.Text('Species', size=(30, 1), font=("Helvetica", 12))],
    [sg.Radio('Human', "RADIO1", key='rad_hs', default=True), sg.Radio('Mouse', "RADIO1", key='rad_mm')],

    [sg.Text('Additional genes', size=(30, 1), font=("Helvetica", 12))],
    [sg.MLine(default_text=extra_gene_text, size=(35, 3), key='additional_genes')],

    [sg.Checkbox('Link TRA/TRB chains', key='scheduler.enabled', enable_events=True),
     sg.InputOptionMenu(link_choices, key='linker_choice', default_value=link_choices[0])],

    [sg.Button('Run Stitchr', size=(box_width, 1))],

    [sg.InputText(key='Export stitched TCRs', do_not_clear=False, enable_events=True, visible=False,
                  size=(box_width, 1)),
     sg.FileSaveAs(file_types=(("FASTA files", "*.fasta"),), default_extension='.fasta',
                  size=(box_width, 1), button_text='Export stitched TCRs')],

    [sg.Button('Exit', size=(box_width, 1))]

]

# Alpha column
col2 = [

    [sg.Text('TRA')],

    [sg.Text('TRAV')], [sg.InputText('', key='TRAV', size=sz)],

    [sg.Text('TRAJ')], [sg.InputText('', key='TRAJ', size=sz)],

    [sg.Text('TRA CDR3')], [sg.InputText('', key='TRA_CDR3', size=sz)],

    [sg.Text('TRA name')], [sg.InputText('', key='TRA_name', size=sz)],

    [sg.Text('5\'/3\' sequences')],
    [sg.InputText('', key='TRA_5_prime_seq', size=(int(box_width/3), 1)),
     sg.InputText('', key='TRA_3_prime_seq', size=(int(box_width/3), 1))],

    [sg.Text('TRA out')], [sg.MLine(default_text='', size=(50, 20), key='TRA_out')]

]

# Beta column
col3 = [

    [sg.Text('TRB')],

    [sg.Text('TRBV')], [sg.InputText('', key='TRBV', size=sz)],

    [sg.Text('TRBJ')], [sg.InputText('', key='TRBJ', size=sz)],

    [sg.Text('TRB CDR3')], [sg.InputText('', key='TRB_CDR3', size=sz)],

    [sg.Text('TRB name')], [sg.InputText('', key='TRB_name', size=sz)],

    [sg.Text('5\'/3\' sequences')],
    [sg.InputText('', key='TRB_5_prime_seq', size=half_sz),
     sg.InputText('', key='TRB_3_prime_seq', size=half_sz)],

    [sg.Text('TRB out')], [sg.MLine(default_text='', size=(50, 20), key='TRB_out')]

]

layout = [[sg.Column(col1, element_justification='l'),
           sg.Column(col2, element_justification='l'),
           sg.Column(col3, element_justification='l')]]

window = sg.Window("stitchr", layout)

example_data = {
    'HUMAN': {
        'TRAV': 'TRAV12-2',
        'TRAJ': 'TRAJ23',
        'TRA_CDR3': 'CAVNFGGGKLIF',
        'TRA_out': '',
        'TRA_name': 'DMF5 (3QEU) TRA',
        'TRBV': 'TRBV6-4',
        'TRBJ': 'TRBJ1-1',
        'TRB_CDR3': 'CASSLSFGTEAFF',
        'TRB_name': 'DMF5 (3QEU) TR',
        'TRB_out': ''},
    'MOUSE': {
        'TRAV': 'TRAV14-1',
        'TRAJ': 'TRAJ33',
        'TRA_CDR3': 'CAASDNYQLIW',
        'TRA_out': '',
        'TRA_name': 'OT-I TRA',
        'TRBV': 'TRBV12-1',
        'TRBJ': 'TRBJ2-7',
        'TRB_CDR3': 'CASSRANYEQYF',
        'TRB_name': 'OT-I TRB',
        'TRB_out': ''}
}

while True:
    event, values = window.read()

    # Determine species
    species = check_species(values)

    if event == 'Example data':

        for field in example_data[species]:
            window[field].update(example_data[species][field])

    elif event == 'Reset form':

        for field in example_data[species]:
            window[field].update('')

        outputs = coll.defaultdict()

    elif event == "Upload TCR details":

        # This section uses the Thimble format input file to populate the GUI fields
        check_is_file(values['uploaded-tcr'])

        with open(values['uploaded-tcr'], 'r') as in_file:
            line_count = 0

            for line in in_file:
                bits = line.replace('\n', '').replace('\r', '').split('\t')

                # Use header line to check it's the right file format
                if line_count == 0:
                    if bits != th.in_headers:
                        raise IOError("Input TCR file doesn't have expected columns. Refer to template.")

                # Then use the data line to replace any matching entries on the webform
                elif line_count == 1:
                    for x in range(len(th.in_headers)):
                        if x == 0:  # Have to add TCR name field individually, as it's only provided once per line
                            window['TRA_name'].update(bits[x])
                            window['TRB_name'].update(bits[x])
                        elif th.in_headers[x] in example_data[species]:  # Then only update those fields that are shared
                            window[th.in_headers[x]].update(bits[x])

                elif line_count > 1:
                    warnings.warn("More than one data line detected in input TCR file. Ignoring lines after first.")
                    break

                line_count += 1

        print(values["uploaded-tcr"])

    elif event == 'Run Stitchr':

        # Disable stitchr button while code is running
        window['Run Stitchr'].update(disabled=True)

        # Loop through both chains, determine which are asked for, and read data in
        codons = fxn.get_optimal_codons('../Data/' + species + '/kazusa.txt', species)
        outputs = coll.defaultdict()

        for chain in ['TRA', 'TRB']:

            if values[chain + 'V'] and values[chain + 'J'] and values[chain + '_CDR3']:
                tcr_dat, functionality = fxn.get_imgt_data(chain, st.gene_types, species)

                # TODO here would be the place to add additional sequences
                # if values['additional_genes'] != extra_gene_text:

                tcr_bits = {'v': values[chain + 'V'], 'j': values[chain + 'J'], 'cdr3': values[chain + '_CDR3'],
                                'skip_c_checks': False, 'species': species,
                                'name': values[chain + '_name'].replace(' ', '_')}

                for end in ['5', '3']:
                    if values[chain + '_' + end + '_prime_seq']:
                        tcr_bits[end + '_prime_seq'] = values[chain + '_' + end + '_prime_seq']

                tcr_bits = fxn.autofill_input(tcr_bits, chain)

                # Run the stitching
                outputs[chain + '_out_list'], outputs[chain + '_stitched'], outputs[chain + '_offset'] = st.stitch(
                    tcr_bits, chain, tcr_dat, functionality, codons, 3)

                outputs[chain + '_out_str'] = '|'.join(outputs[chain + '_out_list'])
                outputs[chain + '_fasta'] = fxn.fastafy('nt|' + outputs[chain + '_out_str'],
                                                                outputs[chain + '_stitched'])

                window[chain + '_out'].update(outputs[chain + '_fasta'])

        # Re-enable stitchr button once completed
        window['Run Stitchr'].update(disabled=False)

    elif event == 'Export stitched TCRs':

        # Only need to bother saving something if there's a stitched TCR to save
        if 'outputs' in dir():
            if len(outputs) > 0:
                out_str = ''
                for chain in ['TRA', 'TRB']:
                    if outputs[chain + '_fasta']:
                        out_str += outputs[chain + '_fasta']

                # TODO also add linked FASTA saved here

                out_file = values['Export stitched TCRs']
                if out_file:
                    with open(out_file, 'w') as out_file:
                        out_file.write(out_str)

    elif event in ('Exit', None):
        break

window.close()

# TODO add error handling
# TODO add linking
# TODO incorporate provided additional genes
# TODO add logo? Might need to make it into a gif
