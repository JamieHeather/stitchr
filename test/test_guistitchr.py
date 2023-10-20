import filecmp

class TestOutFiles:
    """Ensures exported FASTAs are same as baseline versions.

    To test, run the GUI with the inputs below, exporting to output FASTAs.
    Ensure Link chains=P2A and Link order=BA or GD.

    test_alphabeta
    With TRA/TRB, Link chains=P2A, Link order=BA:
        1. 'Example data'
        2. 'Run Stitchr'
        3. 'Export output': stitchr/test/new_ab.fasta

    test_gammadelta
    With TRG/TRD, Link chains=P2A, Link order=GD:
        1. 'Example data'
        2. 'Run Stitchr'
        3. 'Export output': stitchr/test/new_gd.fasta

    test_template
        1. 'Find TCR input file': stitchr/templates/gui_input_example_human.tsv
        2. 'Upload TCR details'
        3. 'Run Stitchr'
        4. 'Export output': stitchr/test/new_template.fasta
    """

    def test_alphabeta(self):
        # filecmp.cmp() should return True
        assert filecmp.cmp('baseline_ab.fasta', 'new_ab.fasta')

    def test_gammadelta(self):
        # filecmp.cmp() should return True
        assert filecmp.cmp('baseline_gd.fasta', 'new_gd.fasta')

    def test_template(self):
        # filecmp.cmp() should return True
        assert filecmp.cmp('baseline_template.fasta', 'new_template.fasta')