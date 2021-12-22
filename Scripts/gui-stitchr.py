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

__version__ = '0.3.2'
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


def read_fasta_box(fasta_text):
    """
    :param fasta_text: Contents of text box containing FASTA format text
    """
    header, seq = None, []
    for line in fasta_text:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                yield(header.replace('>', ''), ''.join(seq))
            header, seq = line, []
        else:
            seq.append(line)
    if header:
        yield(header, ''.join(seq))


extra_gene_text = ">TCRgenename*01\nATCG...\n"
box_width = 70
sz = (int(box_width * 0.9), 1)
half_sz = (int(box_width * 0.44), 1)          # For half a column
third_sz = (int(box_width / 3 - 1.3), 1)          # For one third of a column
quart_sz = (int(box_width / 4 - 1.1), 1)          # For one quarter of a column
out_box_font = ('Courier New', 10)

linkers = fxn.get_linker_dict()
link_choices = list(linkers.keys()) + ['Custom']

fnt = 'Arial'

# General interface column
col1 = [

    [sg.Button('Example data', size=quart_sz), sg.Button('Reset form', size=quart_sz)],

    [sg.FileBrowse(key="uploaded_tcr", size=quart_sz, button_text='Find TCR input file'),
     sg.Button("Upload TCR details", size=quart_sz)],

    [sg.Text('Species', size=quart_sz, font=(fnt, 12))],
    [sg.Radio('Human', "RADIO1", key='rad_hs', default=True), sg.Radio('Mouse', "RADIO1", key='rad_mm')],

    [sg.Text('Additional genes', size=half_sz, font=(fnt, 12))],
    [sg.MLine(default_text=extra_gene_text, size=(int(box_width / 2 - 1.3), 3), key='additional_genes')],

    [sg.Checkbox('Link TRA/TRB', key='chk_linker', enable_events=True, size=quart_sz),
     sg.Combo(link_choices, key='linker_choice', default_value=link_choices[0], size=quart_sz, enable_events=True)],

    [sg.InputText('', key='custom_linker', visible=False, size=third_sz)],

    [sg.Text('Link order', size=third_sz, font=(fnt, 12))],
    [sg.Radio('AB', "RADIO2", key='rad_AB'), sg.Radio('BA', "RADIO2", key='rad_BA', default=True)],

    [sg.Checkbox('Seamless stitching', key='chk_seamless', enable_events=True, font=(fnt, 12))],

    [sg.Button('Run Stitchr', size=(int(box_width/4), 2), font=(fnt, 20))],

    [sg.InputText(key='Export output', do_not_clear=False, enable_events=True, visible=False,
                  size=quart_sz),
    sg.FileSaveAs(file_types=(("FASTA files", "*.fasta"),), default_extension='.fasta',
                  size=quart_sz, button_text='Export output'),

    sg.Button('Exit', size=quart_sz)],

    [sg.Text('Linked out', key='linked_out_text')],
    [sg.MLine(default_text='', size=(int(box_width/2), 10), key='linked_out',font=out_box_font)],

    [sg.Text('Linked log', key='linked_log_text')],
    [sg.MLine(default_text='', size=(int(box_width / 2), 5), key='linked_log', font=out_box_font)],
]

# Alpha column
col2 = [

    [sg.Text('Alpha chain TCR parameters', font=(fnt, 16))],

    [sg.Text('TRAV')], [sg.InputText('', key='TRAV', size=sz)],

    [sg.Text('TRAJ')], [sg.InputText('', key='TRAJ', size=sz)],

    [sg.Text('TRA CDR3 junction')], [sg.InputText('', key='TRA_CDR3', size=sz)],

    [sg.Text('TRA name')], [sg.InputText('', key='TRA_name', size=sz)],

    [sg.Text('TRA leader\t\t\t TRAC')],


    [sg.InputText('', key='TRA_leader', size=half_sz),
     sg.InputText('', key='TRAC', size=half_sz)],

    [sg.Text('5\' sequence\t\t\t 3\' sequence')],
    [sg.InputText('', key='TRA_5_prime_seq', size=half_sz),
     sg.InputText('', key='TRA_3_prime_seq', size=half_sz)],

    [sg.Text('TRA out')],
    [sg.MLine(default_text='', size=(box_width-9, 20), key='TRA_out', font=out_box_font)],

    [sg.Text('TRA log', key='TRA_log_text')],
    [sg.MLine(default_text='', size=(box_width-9, 5), key='TRA_log', font=out_box_font)]

]

# Beta column
col3 = [

    [sg.Text('Beta chain TCR parameters', font=(fnt, 16))],

    [sg.Text('TRBV')], [sg.InputText('', key='TRBV', size=sz)],

    [sg.Text('TRBJ')], [sg.InputText('', key='TRBJ', size=sz)],

    [sg.Text('TRB CDR3 junction')], [sg.InputText('', key='TRB_CDR3', size=sz)],

    [sg.Text('TRB name')], [sg.InputText('', key='TRB_name', size=sz)],

    [sg.Text('TRB leader\t\t\t TRBC')],

    [sg.InputText('', key='TRB_leader', size=half_sz),
     sg.InputText('', key='TRBC', size=half_sz)],

    [sg.Text('5\' sequence\t\t\t 3\' sequence')],

    [sg.InputText('', key='TRB_5_prime_seq', size=half_sz),
     sg.InputText('', key='TRB_3_prime_seq', size=half_sz)],

    [sg.Text('TRB out')],
    [sg.MLine(default_text='', size=(box_width-9, 20), key='TRB_out', font=out_box_font)],

    [sg.Text('TRB log', key='TRB_log_text')],
    [sg.MLine(default_text='', size=(box_width-9, 5), key='TRB_log', font=out_box_font)]

]

layout = [[sg.Column(col1, element_justification='c'),
           sg.Column(col2, element_justification='l'),
           sg.Column(col3, element_justification='l')]]

window = sg.Window("stitchr", layout, finalize=True)

window.bind('<Escape>', 'Exit')

example_data = {
    'HUMAN': {
        'TRAV': 'TRAV12-2',
        'TRAJ': 'TRAJ23',
        'TRA_CDR3': 'CAVNFGGGKLIF',
        'TRA_out': '',
        'TRA_name': 'DMF5 (3QEU) TRA',
        'TRA_leader': '',
        'TRAC': '',
        'TRA_5_prime_seq': '',
        'TRA_3_prime_seq': '',
        'TRBV': 'TRBV6-4',
        'TRBJ': 'TRBJ1-1',
        'TRB_CDR3': 'CASSLSFGTEAFF',
        'TRB_name': 'DMF5 (3QEU) TRB',
        'TRB_leader': '',
        'TRBC': '',
        'TRB_5_prime_seq': '',
        'TRB_3_prime_seq': '',
        'TRB_out': ''},
    'MOUSE': {
        'TRAV': 'TRAV14-1',
        'TRAJ': 'TRAJ33',
        'TRA_CDR3': 'CAASDNYQLIW',
        'TRA_out': '',
        'TRA_name': 'OT-I TRA',
        'TRA_leader': '',
        'TRAC': '',
        'TRA_5_prime_seq': '',
        'TRA_3_prime_seq': '',
        'TRBV': 'TRBV12-1',
        'TRBJ': 'TRBJ2-7',
        'TRB_CDR3': 'CASSRANYEQYF',
        'TRB_name': 'OT-I TRB',
        'TRB_leader': '',
        'TRBC': '',
        'TRB_5_prime_seq': '',
        'TRB_3_prime_seq': '',
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

        for field in ['linked_out', 'linked_out_text',
                      'linked_log', 'linked_log_text',
                      'TRA_log', 'TRA_log_text',
                      'TRB_log', 'TRB_log_text']:
            if field.endswith('log') or field.endswith('out'):
                window[field].update('')

        outputs = coll.defaultdict()

        window['rad_AB'].update(value=False)
        window['rad_BA'].update(value=True)
        window['additional_genes'].update(extra_gene_text)

    elif event == "Upload TCR details":

        # This section uses the Thimble format input file to populate the GUI fields

        # First check file details
        if not os.path.isfile(values['uploaded_tcr']):
            sg.Popup('Please use \'Find TCR input file\' button to try again.', title = 'TCR file not found')
        else:

            with open(values['uploaded_tcr'], 'r') as in_file:
                line_count = 0

                for line in in_file:
                    bits = line.replace('\n', '').replace('\r', '').split('\t')

                    # Use header line to check it's the right file format
                    if line_count == 0:
                        if bits != th.in_headers:
                            sg.Popup("Input TCR file doesn't have expected columns.\n"
                                     "Please refer to template and try again.",
                                     title="TCR file error")
                            break

                    # Then use the data line to replace any matching entries on the webform
                    elif line_count == 1:
                        for x in range(len(th.in_headers)):

                            # Have to add TCR name field individually, as it's only provided once per line
                            if x == 0:
                                window['TRA_name'].update(bits[x])
                                window['TRB_name'].update(bits[x])

                            # Then only update those fields that are shared
                            elif th.in_headers[x] in example_data[species]:
                                window[th.in_headers[x]].update(bits[x])

                            # Update link order...
                            elif th.in_headers[x] == 'Link_order':
                                window['chk_linker'].update(value=True)

                                if bits[x]:
                                    if bits[x] == 'AB':
                                        window['rad_AB'].update(value=True)
                                        window['rad_BA'].update(value=False)
                                    elif bits[x] == 'BA':
                                        window['rad_AB'].update(value=False)
                                        window['rad_BA'].update(value=True)
                                    else:
                                        raise warnings.warn("Invalid link order input - must be AB or BA.")

                            # ... and linker sequence
                            elif th.in_headers[x] == 'Linker':
                                window['chk_linker'].update(value=True)

                                if bits[x]:
                                    if bits[x] in linkers.keys():
                                        window['linker_choice'].update(value=bits[x])
                                    else:
                                        window['linker_choice'].update('Custom')
                                        window['custom_linker'].update(bits[x], visible=True)

                    elif line_count > 1:
                        warnings.warn("More than one data line detected in input TCR file. Ignoring lines after first.")
                        break

                    line_count += 1

    elif event == 'Run Stitchr':
        
        warning_msgs = coll.defaultdict(str)

        window['linked_out'].update('')
        window['linked_log'].update('')

        # Disable stitchr button while code is running
        window['Run Stitchr'].update(disabled=True)

        # Loop through both chains, determine which are asked for, and read data in
        codons = fxn.get_optimal_codons('../Data/' + species + '/kazusa.txt', species)
        outputs = coll.defaultdict()

        # If additional genes provided, read in and run rudimentary checks
        if values['additional_genes'] != extra_gene_text + '\n':
            outputs['additional_fastas_raw'] = [x for x in read_fasta_box(
                values['additional_genes'].split('\n') + ['>\n'])][:-1]

            # Check no redundant gene names
            if len(list(set([x[0] for x in outputs['additional_fastas_raw']]))) != \
                    len(outputs['additional_fastas_raw']):
                window['additional_genes'].update("Multiple FASTAs detected with the same identifiername name.\n"
                                                  "Additional genes ignored; correct and retry")

            # Check they all have allele numbers to match the expected gene name format
            outputs['additional_fastas'] = []
            for extra_gene in outputs['additional_fastas_raw']:
                if '*' in extra_gene[0]:
                    outputs['additional_fastas'].append(extra_gene)
                else:
                    outputs['additional_fastas'].append((extra_gene[0] + '*01', extra_gene[1]))

                # Also throw in an alert if non-DNA characters used
                if not fxn.dna_check(extra_gene[1]):
                    warnings.warn("Warning: user-provided gene " + extra_gene[0] + " contains non-DNA sequences.")

        # Check if seamless stitching selected
        if values['chk_seamless']:
            seamless = True
        else:
            seamless = False

        # Then stitch each individual chain...
        for chain in ['TRA', 'TRB']:

            window[chain + '_out'].update('')
            window[chain + '_log'].update('')

            with warnings.catch_warnings(record=True) as chain_log:
                warnings.simplefilter("always")

                if values[chain + 'V'] and values[chain + 'J'] and values[chain + '_CDR3']:
    
                    try:
                        tcr_dat, functionality, partial = fxn.get_imgt_data(chain, st.gene_types, species)
    
                        # If additional genes provided, just add them to all possible gene segment types
                        if values['additional_genes'] != extra_gene_text + '\n':
                            for extra_gene in outputs['additional_fastas']:
                                gene, allele = extra_gene[0].split('*')
    
                                for gene_type in tcr_dat.keys():
    
                                    if gene not in tcr_dat[gene_type]:
                                        tcr_dat[gene_type][gene] = coll.defaultdict(list)
    
                                    if allele in tcr_dat[gene_type][gene]:
                                        raise warnings.warn("User provided gene/allele combination " + extra_gene[0] +
                                                      " already exists in TCR germline data. Please change and try.")
                                    else:
                                        tcr_dat[gene_type][gene][allele] = extra_gene[1].upper()
                                        functionality[gene][allele] = '?'
    
                        tcr_bits = {'v': values[chain + 'V'], 'j': values[chain + 'J'], 'cdr3': values[chain + '_CDR3'],
                                    'skip_c_checks': False, 'species': species, 'seamless': seamless,
                                    'name': values[chain + '_name'].replace(' ', '_')}
    
                        # Can't do C checks if user providing genes, as it may be a C
                        if values['additional_genes'] != extra_gene_text + '\n':
                            tcr_bits['skip_c_checks'] = True
    
                        for end in ['5', '3']:
                            if values[chain + '_' + end + '_prime_seq']:
                                tcr_bits[end + '_prime_seq'] = values[chain + '_' + end + '_prime_seq']
    
                        for section in ['_leader', 'C']:
                            if values[chain + section]:
                                tcr_bits[th.convert_fields[chain + section]] = values[chain + section]
    
                        tcr_bits = fxn.autofill_input(tcr_bits, chain)
    
                        # Run the stitching
                        outputs[chain + '_out_list'], outputs[chain + '_stitched'], outputs[chain + '_offset'] = st.stitch(
                            tcr_bits, chain, tcr_dat, functionality, partial, codons, 3)
    
                        outputs[chain + '_out_str'] = '|'.join(outputs[chain + '_out_list'])
                        outputs[chain + '_fasta'] = fxn.fastafy('nt|' + outputs[chain + '_out_str'],
                                                                outputs[chain + '_stitched'])

                        window[chain + '_out'].update(outputs[chain + '_fasta'])
    
                    except Exception as message:
                        warning_msgs[chain + '_out'] = str(message)

                elif values[chain + 'V'] or values[chain + 'J'] or values[chain + '_CDR3']:
                    warnings.warn('V gene, J gene, and CDR3 sequence are all required to stitch a TCR chain.')

            warning_msgs[chain + '_out'] += ''.join([str(chain_log[x].message) for x in range(len(chain_log))
                                                    if 'DeprecationWarning' not in str(chain_log[x].category)])

            window[chain + '_log'].update(warning_msgs[chain + '_out'])

        # ... and if asked for, link together
        if values['chk_linker']:

            with warnings.catch_warnings(record=True) as link_log:
                warnings.simplefilter("always")

                # Only link if both chains present
                try:
                    if 'TRA_out_str' in outputs and 'TRB_out_str' in outputs:

                        # Determine order
                        if values['rad_AB'] and not values['rad_BA']:
                            tr1, tr2 = 'A', 'B'
                        elif values['rad_BA'] and not values['rad_AB']:
                            tr1, tr2 = 'B', 'A'
                        else:
                            raise warnings.warn("Undetermined link order.")

                        # Stick together, first verifying linker
                        outputs['linker'] = values['linker_choice']
                        if outputs['linker'] == 'Custom':
                            if values['custom_linker']:
                                linkers['Custom'] = values['custom_linker']
                            else:
                                window['linked_log'].update("Cannot output linked sequence: custom linker chosen, "
                                                            "but not provided")

                        outputs['linker_seq'] = fxn.get_linker_seq(outputs['linker'], linkers)

                        outputs['linked'] = outputs['TR' + tr1 + '_stitched'] + \
                                            outputs['linker_seq'] + \
                                            outputs['TR' + tr2 + '_stitched']

                        outputs['linked_header'] = '_'.join([outputs['TR' + tr1 + '_out_str'],
                                                             outputs['linker'],
                                                             outputs['TR' + tr2 + '_out_str']])

                        outputs['linked_fasta'] = fxn.fastafy(outputs['linked_header'], outputs['linked'])

                        window['linked_out'].update(outputs['linked_fasta'])

                    else:
                        raise warnings.warn("Valid TRA and TRB chains required for linking.")

                except Exception as message:
                    warning_msgs['linked'] += str(message)

            warning_msgs['linked_out'] += ''.join([str(link_log[x].message) for x in range(len(link_log))
                                                    if 'DeprecationWarning' not in str(link_log[x].category)])

            if warning_msgs['linked_out']:
                window['linked_log'].update(warning_msgs['linked_out'])

        # Re-enable stitchr button once completed
        window['Run Stitchr'].update(disabled=False)

    elif event == 'Export output':

        # Only need to bother saving something if there's a stitched TCR to save
        if 'outputs' in dir():
            if len(outputs) > 0:
                out_str = ''
                for chain in ['TRA', 'TRB']:
                    if outputs[chain + '_fasta']:
                        out_str += outputs[chain + '_fasta']

                if values['chk_linker']:
                    out_str += outputs['linked_fasta']

                out_file = values['Export output']
                if out_file:
                    with open(out_file, 'w') as out_file:
                        out_file.write(out_str)


    elif event == 'linker_choice':

        window['chk_linker'].update(value=True)

        if values['linker_choice'] == 'Custom':
            window['custom_linker'].update(visible=True)
        else:
            window['custom_linker'].update(visible=False)

    elif event in ('Exit', None):
        break

window.close()

# TODO output a file of warnings?