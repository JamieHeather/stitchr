
import sys
import os
import pytest
import subprocess
from Stitchr import stitchrfunctions as fxn


def empty_file(file_name):
    with open(file_name, 'w') as out_file:
        out_file.write('')
    return file_name


def test_downloaded_data_handling():
    # Need to have data downloaded to test various functions; arbitrarily use canine and mouse
    dl_cmd = 'stitchrdl -s dog'
    subprocess.call(dl_cmd, shell=True)
    dl_cmd = 'stitchrdl -s mouse'
    subprocess.call(dl_cmd, shell=True)

    assert 'DOG' in fxn.find_species_covered()
    assert 'MOUSE' in fxn.find_species_covered()

    tmp_file = empty_file('temp_dog.txt')
    assert fxn.infer_species(tmp_file) == 'DOG'
    os.remove(tmp_file)

    tmp_file = empty_file('temp_dog_mouse.txt')
    assert fxn.infer_species(tmp_file) == ''
    os.remove(tmp_file)


def test_file_handling():
    in_header, in_seq = ['test_header', 'ACGT']
    tmp_file = 'temp.fasta'
    with open(tmp_file, 'w') as out_file:
        out_file.write(fxn.fastafy(in_header, in_seq))

    with fxn.opener(tmp_file) as in_file:
        for header, seq, null in fxn.read_fa(in_file):
            assert in_header == header and in_seq == seq, "Failed FASTA handling."

    os.remove(tmp_file)


def test_get_chain():
    assert fxn.get_chain('TRAV', 'TRAJ') == 'TRA'
    assert fxn.get_chain('TRDV', 'TRAJ') == 'TRA'
    assert fxn.get_chain('TRBV', 'TRBJ') == 'TRB'
    assert fxn.get_chain('TRDV', 'TRDJ') == 'TRD'
    assert fxn.get_chain('TRGV', 'TRGJ') == 'TRG'
    assert fxn.get_chain('IGHV', 'IGHJ') == 'IGH'
    assert fxn.get_chain('IGLV', 'IGLJ') == 'IGL'
    assert fxn.get_chain('IGKV', 'IGKJ') == 'IGK'
    assert not fxn.get_chain('TRAV', 'TRAJ') == 'TRB'
    with pytest.raises(ValueError):
        fxn.get_chain('TRAV', 'TRBJ')


def test_get_imgt_data():
    assert fxn.get_imgt_data('TRA', list(fxn.regions.values()), 'DOG')

    with pytest.raises(ValueError):
        fxn.get_imgt_data('TEST', list(fxn.regions.values()), 'DOG')

    with pytest.raises(IOError):
        fxn.get_imgt_data('TRA', list(fxn.regions.values()), 'UNICORN')

    imgt_dat = fxn.get_imgt_data('TRA', list(fxn.regions.values()), 'MOUSE')
    assert [isinstance(x, dict) for x in imgt_dat]
    assert list(imgt_dat[0].keys()) == ['LEADER', 'VARIABLE', 'JOINING', 'CONSTANT']
    assert list(imgt_dat[1].keys())[0].startswith('TR')

    preferred = fxn.get_preferred_alleles('templates/preferred-alleles/mouse_c57-bl6_example.tsv',
                                          list(fxn.regions.values()), imgt_dat[0], imgt_dat[2], 'TRA')
    assert preferred['LEADER']['TRAV15D-3'] == '02'



aa = ['*', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
def test_translation_handling():
    codon_file = os.path.join(fxn.data_dir, 'kazusa', 'DOG.txt')
    assert fxn.get_codon_frequencies(codon_file)
    codon_freqs = fxn.get_codon_frequencies(codon_file)
    assert [x for x in aa if x in codon_freqs.keys()] == aa
    assert [x for x in aa if x in list(set(fxn.codons.values()))] == aa
    assert fxn.get_optimal_codons(codon_file, 'DOG')
    codons = fxn.get_optimal_codons(codon_file, 'DOG')

    assert (fxn.translate_nt('AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCAC'
                             'CCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGG'
                             'GGGTGTAGTCGTGGTTTAATACTAGTATTCATCCTCGTCTTGATGCTGGTGTTTATTCTTGTTT') ==
            'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF')
    assert fxn.translate_nt('naaaaa') == 'xK'
    with pytest.raises(OSError):
        fxn.translate_nt('123')

    assert (fxn.rev_translate(''.join(aa).replace('*', ''), codons) ==
            'GCCTGCGACGAGTTCGGCCACATCAAGCTGATGAACCCCCAGAGGAGCACCGTGTGGTAC')


def test_misc_str_handling():
    assert fxn.strip_functionality('()[]test') == 'test'
    assert fxn.tidy_n_term('aaacccggg') == ('aaacccggg', 'KPG')
    assert fxn.tidy_n_term('aaacccgggtt') == ('aaacccggg', 'KPG')
    assert fxn.find_stop('0123456789') == 10
    assert fxn.find_stop('01234*6789') == 5
    assert fxn.dna_check('acgt')
    assert not fxn.dna_check('dna')
    assert fxn.check_suffix_prefix('thisdoesnt', 'overlap') == ''
    assert fxn.check_suffix_prefix('thisdoes', 'doesthough') == 'does'


def test_linker_handling():
    assert fxn.get_linker_dict()
    linkers = fxn.get_linker_dict()
    assert linkers['P2A'] == 'GGCAGCGGCGCCACCAACTTCAGCCTGCTGAAGCAGGCCGGCGACGTGGAGGAGAACCCCGGCCCC'
    linker_seq = fxn.get_linker_seq('P2A', linkers)
    assert linker_seq == 'GGCAGCGGCGCCACCAACTTCAGCCTGCTGAAGCAGGCCGGCGACGTGGAGGAGAACCCCGGCCCC'


# TODO STILL
    # actuall call an installed stitchrdl and test (can't currently test it due to intern
    # def get_additional_genes(imgt_data, imgt_functionality):
    # def tidy_c_term(c_term_nt, skip, c_region_motifs, c_gene):
    # def determine_v_interface(cdr3aa, n_term_nuc, n_term_amino):
    # def find_cdr3_c_term(cdr3_chunk, j_seq, strict):
    # def determine_j_interface(cdr3_cterm_aa, c_term_nuc, c_term_amino, gl_nt_j_len, j_warning_threshold):
    # def get_j_motifs(species):
    # def get_c_motifs(species):
    # def find_v_overlap(v_germline, nt_cdr3):
    # def find_j_overlap(nt_cdr3, j_germline):

    # def tweak_thimble_input(stitch_dict): # Todo move to a thimble test file?
