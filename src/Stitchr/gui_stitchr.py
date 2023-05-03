# -*- coding: utf-8 -*-

"""
gui_stitchr.py

A graphical user interface for stitchr, powered by PySimpleGUI

"""


import PySimpleGUI as sg
import os
from . import stitchrfunctions as fxn
from . import stitchr as st
from . import thimble as th
import collections as coll
import warnings


__version__ = '1.3.2'
__author__ = 'Jamie Heather'
__email__ = 'jheather@mgh.harvard.edu'


def read_fasta_box(fasta_text):
    """
    :param fasta_text: Contents of text box containing FASTA format text
    """
    header, seq = None, []
    for line in fasta_text:
        line = line.rstrip()
        if line.startswith(">"):
            if header:
                yield header.replace('>', ''), ''.join(seq)
            header, seq = line, []
        else:
            seq.append(line)
    if header:
        yield header, ''.join(seq)


def switch_receptors(text, current_receptor):
    """
    :param text: a str containing a reference to a receptor which needs to be replaced
    :param current_receptor: str of 'TRA/TRB' or 'TRG/TRD', which indicates the required direction of change
    :return: that string with all references to a given locus replaced with the right one
    """
    switched_receptor = text
    changes = [
        ('Alpha', 'Gamma'),
        ('Beta', 'Delta'),
        ('TRA', 'TRG'),
        ('TRB', 'TRD'),
        ('AB', 'GD'),
        ('BA', 'DG')
    ]
    for change in changes:
        if current_receptor == 'TRA/TRB':
            switched_receptor = switched_receptor.replace(change[0], change[1])
        elif current_receptor == 'TRG/TRD':
            switched_receptor = switched_receptor.replace(change[1], change[0])
        else:
            raise ValueError('Unknown receptor configuration detected: ' + current_receptor)
    return switched_receptor


def change_receptors(receptor_str):
    """
    Updates all the actual fields when swapping receptor chains (i.e. a/b to g/d TCRs)
    :param receptor_str: string detailing the CURRENT receptor, i.e. either 'TRA/TRB' or 'TRG/TRD'
    :return: the new receptor used (couresty of switch_receptor)
    """
    # Allow the users to flick between stitching (and thus viewing parameters for) A/B and G/D TCR sequences
    new_receptor = switch_receptors(receptor_str, receptor_str)

    # First change the link order dropdown (which needs to be handled differently)
    window['link_order_choice'].update(values=link_orders[new_receptor], value=link_orders[new_receptor][1])

    # Then the rest of the fields
    fields_to_change = [x for x in window.AllKeysDict if x.endswith('_text')]
    fields_to_change += ['change_receptor']

    for window_field in fields_to_change:

        if window[window_field].Type == 'text':
            current_text = window[window_field].DisplayText
        elif window[window_field].Type == 'button':
            current_text = window[window_field].get_text()

        if window_field != 'change_receptor':
            new_text = switch_receptors(current_text, receptor_str)
        else:
            new_text = switch_receptors(current_text, new_receptor)

        window[window_field].update(new_text)

    # Last thing: change the current receptor being used
    return new_receptor


def upload_tcr_details(path_to_file, receptor_type, stated_species):
    """
    Read suitably formatted TCRs in template format to the GUI
    :param path_to_file: str detailing the full path to the input file
    :param receptor_type: str detailing the receptor in use, i.e. either 'TRA/TRB' or 'TRG/TRD'
    :param stated_species: str of the species on record before the upload, to remain if one not inferred
    :return: strs of receptor type and species again, in case they switch during the course of reading the TCR file in
    """

    # This section uses the Thimble format input file to populate the GUI fields

    # First check file details
    if not os.path.isfile(path_to_file):
        sg.Popup('Please use \'Find TCR input file\' button to try again.', title='TCR file not found')
    else:

        # Try to estimate species from input filename
        species_inference = fxn.infer_species(path_to_file)
        if species_inference:
            inferred_species = species_inference
            window['species_choice'].update(inferred_species)
        else:
            inferred_species = stated_species
            sg.Popup("Cannot infer species name from file name:\nplease set manually.")

        with open(path_to_file, 'r') as in_file:
            line_count = 0

            for line in in_file:
                bits = line.replace('\n', '').replace('\r', '').split('\t')

                # Use header line to check it's the right file format
                if line_count == 0:

                    if bits != th.in_headers[receptor_type]:
                        # Might be other loci
                        if bits == th.in_headers[switch_receptors(receptor_type, receptor_type)]:
                            receptor_type = change_receptors(receptor_type)

                        else:
                            sg.Popup("Input TCR file doesn't have expected columns.\n"
                                     "Please refer to template and try again.",
                                     title="TCR file error")
                            break

                # Then use the data line to replace any matching entries on the webform
                elif line_count == 1:
                    for x in range(len(th.in_headers[receptor_type])):

                        # Have to add TCR name field individually, as it's only provided once per line
                        if x == 0:
                            window['TR1_name'].update(bits[x])
                            window['TR2_name'].update(bits[x])

                        # Then only update those fields that are shared
                        elif 'Link' not in th.in_headers[receptor_type][x]:
                            window[th.locus_to_trx(th.in_headers[receptor_type][x])].update(bits[x])

                        # Update link order...
                        elif th.in_headers[receptor_type][x] == 'Link_order':
                            window['chk_linker'].update(value=True)

                            if bits[x]:
                                if bits[x] == 'AB':
                                    position = 0
                                elif bits[x] == 'BA':
                                    position = 1
                                elif bits[x] == 'GD':
                                    position = 0
                                elif bits[x] == 'DG':
                                    position = 1
                                else:
                                    raise warnings.warn("Invalid link order: " + bits[x])

                                window['link_order_choice'].update(values=link_orders[receptor_type],
                                                                   value=link_orders[receptor_type][position])

                        # ... and linker sequence
                        elif th.in_headers[receptor_type][x] == 'Linker':
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

        return receptor_type, inferred_species


def tidy_values(reference_chain, set_values):
    """
    Make sure the input values for a given chain are plausible
    :param reference_chain: string denoting the current chain
    :param set_values: the dict of values passed to the stitching function
    """
    for field in ['V', 'J', '_CDR3', '_leader', 'C']:
        if set_values[reference_chain + field]:
            set_values[reference_chain + field] = set_values[reference_chain + field].upper()

    return set_values


def main():

    # Define needed/starting variables
    extra_gene_text = ">TCRgenename*01\nATG\n"
    box_width = 70
    sz = (int(box_width * 0.9), 1)
    half_sz = (int(box_width * 0.44), 1)  # For half a column
    third_sz = (int(box_width / 3 - 1.3), 1)  # For one third of a column
    quart_sz = (int(box_width / 4 - 1.1), 1)  # For one quarter of a column
    out_box_font = ('Courier New', 10)
    fnt = 'Arial'

    global linkers
    linkers = fxn.get_linker_dict()
    global link_choices
    link_choices = list(linkers.keys()) + ['Custom']

    species_list = fxn.find_species_covered()

    receptor = 'TRA/TRB'  # Start off with a/b TCRs are the format in use
    receptor_list = ['TRA/TRB', 'TRG/TRD']

    global link_orders
    link_orders = {'TRA/TRB': ['AB', 'BA'],
                   'TRG/TRD': ['GD', 'DG']}

    examples_path = fxn.gui_examples_dir
    preferred = ''
    preferred_button_default = 'Preferred allele file'

    # General interface column
    col1 = [

        [sg.Button('Example data', size=quart_sz), sg.Button('Reset form', size=quart_sz)],

        [sg.FileBrowse(key="uploaded_tcr", size=quart_sz, button_text='Find TCR input file'),
         sg.Button("Upload TCR details", size=quart_sz)],

        [sg.Text('Species', size=half_sz, font=(fnt, 12))],
        [sg.Combo(species_list, key='species_choice', default_value='HUMAN', size=half_sz, enable_events=True)],

        [sg.Button(key='change_receptor', size=half_sz, enable_events=True, button_text='Change to TRG/TRD')],

        [sg.Text('Additional genes', size=half_sz, font=(fnt, 12))],
        [sg.MLine(default_text=extra_gene_text, size=(int(box_width / 2 - 1.3), 3), key='additional_genes')],

        [sg.Input(key='find_preferred_alleles', enable_events=True, visible=False)],
        [sg.FileBrowse(target='find_preferred_alleles', size=half_sz,
                       key='preferred_allele_button', button_text=preferred_button_default)],

        [sg.Checkbox('Link chains', key='chk_linker', enable_events=True, size=quart_sz, font=(fnt, 12)),
         sg.Combo(link_choices, key='linker_choice', default_value=link_choices[0], size=quart_sz, enable_events=True)],

        [sg.InputText('', key='custom_linker', visible=False, size=third_sz)],

        [sg.Text('Link order', size=quart_sz, font=(fnt, 12), justification='c'),
         sg.Combo(link_orders[receptor], key='link_order_choice', default_value=link_orders[receptor][1],
                  size=(8, 1), enable_events=True)],

        [sg.Checkbox('Seamless stitching', key='chk_seamless', enable_events=True, font=(fnt, 12))],

        [sg.Button('Run Stitchr', size=(int(box_width / 4), 2), font=(fnt, 20))],

        [sg.InputText(key='Export output', do_not_clear=False, enable_events=True, visible=False,
                      size=quart_sz),
         sg.FileSaveAs(file_types=(("FASTA files", "*.fasta"),), default_extension='.fasta',
                       size=quart_sz, button_text='Export output'), sg.Button('Exit', size=quart_sz)],

        [sg.Text('Linked out', key='linked_out_text')],
        [sg.MLine(default_text='', size=(int(box_width / 2), 10), key='linked_out', font=out_box_font)],

        [sg.Text('Linked log', key='linked_log_text')],
        [sg.MLine(default_text='', size=(int(box_width / 2), 5), key='linked_log', font=out_box_font)],
    ]

    # Alpha/gamma column
    col2 = [

        [sg.Text('Alpha chain TCR', font=(fnt, 16), key='TR1_title_text')],

        [sg.Text('TRAV', key='TR1_V_text')], [sg.InputText('', key='TR1V', size=sz)],

        [sg.Text('TRAJ', key='TR1_J_text')], [sg.InputText('', key='TR1J', size=sz)],

        [sg.Text('TRA CDR3 junction', key='TR1_CDR3_text')], [sg.InputText('', key='TR1_CDR3', size=sz)],

        [sg.Text('TRA name', key='TR1_name_text')], [sg.InputText('', key='TR1_name', size=sz)],

        [sg.Text('TRA leader', size=half_sz, key='TR1_l_title_text'), sg.Text('TRAC', key='TR1_c_title_text')],

        [sg.InputText('', key='TR1_leader', size=half_sz),
         sg.InputText('', key='TR1C', size=half_sz)],

        [sg.Text('5\' sequence', size=half_sz), sg.Text('3\' sequence')],
        [sg.InputText('', key='TR1_5_prime_seq', size=half_sz),
         sg.InputText('', key='TR1_3_prime_seq', size=half_sz)],

        [sg.Text('TRA out', key='TR1_out_title_text')],
        [sg.MLine(default_text='', size=(box_width - 9, 20), key='TR1_out', font=out_box_font)],

        [sg.Text('TRA log', key='TR1_log_title_text')],
        [sg.MLine(default_text='', size=(box_width - 9, 5), key='TR1_log', font=out_box_font)]

    ]

    # Beta/delta column
    col3 = [

        [sg.Text('Beta chain TCR', font=(fnt, 16), key='TR2_title_text')],

        [sg.Text('TRBV', key='TR2_V_text')], [sg.InputText('', key='TR2V', size=sz)],

        [sg.Text('TRBJ', key='TR2_J_text')], [sg.InputText('', key='TR2J', size=sz)],

        [sg.Text('TRB CDR3 junction', key='TR2_CDR3_text')], [sg.InputText('', key='TR2_CDR3', size=sz)],

        [sg.Text('TRB name', key='TR2_name_text')], [sg.InputText('', key='TR2_name', size=sz)],

        [sg.Text('TRB leader', size=half_sz, key='TR2_l_title_text'), sg.Text('TRBC', key='TR2_c_title_text')],

        [sg.InputText('', key='TR2_leader', size=half_sz),
         sg.InputText('', key='TR2C', size=half_sz)],

        [sg.Text('5\' sequence', size=half_sz), sg.Text('3\' sequence')],

        [sg.InputText('', key='TR2_5_prime_seq', size=half_sz),
         sg.InputText('', key='TR2_3_prime_seq', size=half_sz)],

        [sg.Text('TRB out', key='TR2_out_title_text')],
        [sg.MLine(default_text='', size=(box_width - 9, 20), key='TR2_out', font=out_box_font)],

        [sg.Text('TRB log', key='TR2_log_title_text')],
        [sg.MLine(default_text='', size=(box_width - 9, 5), key='TR2_log', font=out_box_font)]

    ]

    layout = [[sg.Column(col1, element_justification='c'),
               sg.Column(col2, element_justification='l'),
               sg.Column(col3, element_justification='l')]]

    global window
    window = sg.Window("stitchr", layout, finalize=True)

    window.bind('<Escape>', 'Exit')

    fields_to_reset = [
        'TR1V', 'TR1J', 'TR1_CDR3', 'TR1_name', 'TR1_leader', 'TR1C', 'TR1_5_prime_seq', 'TR1_3_prime_seq', 'TR1_out',
        'TR2V', 'TR2J', 'TR2_CDR3', 'TR2_name', 'TR2_leader', 'TR2C', 'TR2_5_prime_seq', 'TR2_3_prime_seq', 'TR2_out']

    convert_chains = {'TRA/TRB': {'TR1': 'TRA', 'TR2': 'TRB'},
                      'TRG/TRD': {'TR1': 'TRG', 'TR2': 'TRD'}}

    while True:
        event, values = window.read()

        # Prevent TypeError when closing via GUI 'X' button
        if event == sg.WINDOW_CLOSED:
            break

        # Determine species
        species = values['species_choice']

        if event == 'Example data':

            example_files = [x for x in os.listdir(examples_path) if x.endswith('tsv') and x[0] not in ['.', '~', '_']]
            example_matches = [x for x in example_files if species in x.upper() and receptor.replace('/', '-') in x]

            if len(example_matches) == 0:
                sg.Popup('No example ' + receptor + ' TCR example files available for species ' + species + '.')

            else:
                if len(example_matches) > 1:
                    sg.Popup('More than one ' + receptor + 'TCR example files available for species ' + species + ':\n'
                             'using the first alphabetically.')
                    example_matches.sort()

                receptor, species = upload_tcr_details(os.path.join(examples_path, example_matches[0]),
                                                       receptor, species)

        elif event == 'change_receptor':

            receptor = change_receptors(receptor)

        elif event == 'Reset form':

            window['species_choice'].update('HUMAN')
            species = 'HUMAN'
            for field in fields_to_reset:
                window[field].update('')

            for field in ['linked_out', 'linked_out_text',
                          'linked_log', 'linked_log_text',
                          'TR1_log', 'TR1_log_text',
                          'TR2_log', 'TR2_log_text']:
                if field.endswith('log') or field.endswith('out'):
                    window[field].update('')

            outputs = coll.defaultdict()
            window['link_order_choice'].update(values=link_orders[receptor], value=link_orders[receptor][1])
            window['additional_genes'].update(extra_gene_text)

            # Reset preferred alleles
            values['find_preferred_alleles'] = ''
            window['preferred_allele_button'].update(preferred_button_default)

        elif event == 'Upload TCR details':

            receptor, species = upload_tcr_details(values['uploaded_tcr'], receptor, species)

        elif event == 'find_preferred_alleles':

            preferred_file = values['find_preferred_alleles']
            window['preferred_allele_button'].update(os.path.basename(preferred_file))

        elif event == 'Run Stitchr':
            warning_msgs = coll.defaultdict(str)

            window['linked_out'].update('')
            window['linked_log'].update('')

            # Disable stitchr button while code is running
            window['Run Stitchr'].update(disabled=True)

            # Loop through both chains, determine which are asked for, and read data in
            codons = fxn.get_optimal_codons('', species)
            outputs = coll.defaultdict()

            # If additional genes provided, read in and run rudimentary checks
            if values['additional_genes'] != extra_gene_text + '\n':
                outputs['additional_fastas_raw'] = [x for x in read_fasta_box(
                    values['additional_genes'].split('\n') + ['>\n'])][:-1]

                # Check no redundant gene names
                if len(list(set([x[0] for x in outputs['additional_fastas_raw']]))) != \
                        len(outputs['additional_fastas_raw']):
                    window['additional_genes'].update("Multiple FASTAs detected with the same identifier name.\n"
                                                      "Additional genes ignored; correct and retry")

                # Check they all have allele numbers to match the expected gene name format
                outputs['additional_fastas'] = []
                for extra_gene in outputs['additional_fastas_raw']:

                    extra_gene = ([x.upper() for x in extra_gene])

                    if '*' in extra_gene[0]:
                        outputs['additional_fastas'].append(extra_gene)
                    else:
                        outputs['additional_fastas'].append((extra_gene[0] + '*01', extra_gene[1]))

                    # Also throw in an alert if non-DNA characters used
                    if not fxn.dna_check(extra_gene[1]):
                        if fxn.dna_check(extra_gene[1]) != extra_gene_text:
                            warnings.warn("Warning: user-provided gene " + extra_gene[0] +
                                          " contains non-DNA sequences.")

            # Check if seamless stitching selected
            if values['chk_seamless']:
                seamless = True
            else:
                seamless = False

            # Then stitch each individual chain...
            for ref_chain in ['TR1', 'TR2']:
                chain = convert_chains[receptor][ref_chain]

                window[ref_chain + '_out'].update('')
                window[ref_chain + '_log'].update('')

                with warnings.catch_warnings(record=True) as chain_log:
                    warnings.simplefilter("always")

                    if values[ref_chain + 'V'] and values[ref_chain + 'J'] and values[ref_chain + '_CDR3']:

                        values = tidy_values(ref_chain, values)

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
                                            raise warnings.warn("User provided gene/allele combination " +
                                                                extra_gene[0] + " already exists in TCR germline data. "
                                                                                "Please change and try.")
                                        else:
                                            tcr_dat[gene_type][gene][allele] = extra_gene[1].upper()
                                            functionality[gene][allele] = '?'

                            if values['find_preferred_alleles']:
                                preferred = fxn.get_preferred_alleles(values['find_preferred_alleles'],
                                                                      list(fxn.regions.values()), tcr_dat,
                                                                      partial, chain)
                            else:
                                preferred = ''

                            tcr_bits = {'v': values[ref_chain + 'V'], 'j': values[ref_chain + 'J'],
                                        'cdr3': values[ref_chain + '_CDR3'],
                                        'skip_c_checks': False, 'species': species, 'seamless': seamless,
                                        'name': values[ref_chain + '_name'].replace(' ', '_'),
                                        'l': values[ref_chain + '_leader'], 'c': values[ref_chain + 'C'],
                                        '5_prime_seq': values[ref_chain + '_5_prime_seq'],
                                        '3_prime_seq': values[ref_chain + '_3_prime_seq']}

                            # Can't do C checks if user providing genes, as it may be a C
                            if values['additional_genes'] != extra_gene_text + '\n':
                                tcr_bits['skip_c_checks'] = True

                            tcr_bits = fxn.autofill_input(tcr_bits, chain)

                            # Run the stitching
                            outputs[ref_chain + '_out_list'], \
                            outputs[ref_chain + '_stitched'], \
                            outputs[ref_chain + '_offset'] = st.stitch(tcr_bits, tcr_dat, functionality,
                                                                       partial, codons, 3, preferred)

                            outputs[ref_chain + '_out_str'] = '|'.join(outputs[ref_chain + '_out_list'])
                            outputs[ref_chain + '_fasta'] = fxn.fastafy('nt|' + outputs[ref_chain + '_out_str'],
                                                                        outputs[ref_chain + '_stitched'])

                            window[ref_chain + '_out'].update(outputs[ref_chain + '_fasta'])

                        except Exception as message:
                            warning_msgs[ref_chain + '_out'] = str(message)

                    elif values[ref_chain + 'V'] or values[ref_chain + 'J'] or values[ref_chain + '_CDR3']:
                        warnings.warn('V gene, J gene, and CDR3 sequence are all required to stitch a TCR chain.')

                warning_msgs[ref_chain + '_out'] += ''.join([str(chain_log[x].message) for x in range(len(chain_log))
                                                        if 'DeprecationWarning' not in str(chain_log[x].category)])

                window[ref_chain + '_log'].update(warning_msgs[ref_chain + '_out'])

            # ... and if asked for, link together
            if values['chk_linker']:

                with warnings.catch_warnings(record=True) as link_log:
                    warnings.simplefilter("always")

                    # Only link if both chains present
                    try:
                        if 'TR1_out_str' in outputs and 'TR2_out_str' in outputs:

                            # Determine order
                            if values['link_order_choice'] == 'BA' or values['link_order_choice'] == 'DG':
                                tr1, tr2 = '2', '1'
                            elif values['link_order_choice'] == 'AB' or values['link_order_choice'] == 'GD':
                                tr1, tr2 = '1', '2'
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
                            raise warnings.warn("Two valid chains required for linking.")

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
                    for chain in ['TR1', 'TR2']:
                        if chain + '_fasta' in outputs:
                            out_str += outputs[chain + '_fasta']

                    if values['chk_linker'] and ('TR1_out_str' in outputs and 'TR2_out_str' in outputs):
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
