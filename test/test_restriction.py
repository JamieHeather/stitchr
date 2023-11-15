from Stitchr import stitchrfunctions as fxn



class TestRestrict:
    """
    Ensures that the output of the restriction functions is the same as input.  No preparation needed.
    """
    def process(str, enzymes):
        sites = fxn.check_restricts(str, enzymes)
        sequence = fxn.wobble(str, sites, enzymes)
        sequence = fxn.translate_nt(sequence)
        return(sequence)
    
    def test_noSites(self):
        """
        Checks to see if two sequences are identical w/ and w/out editing with our process and no restriction sites
        """
        str = "ATGAAATCCTTGAGAGTTTTACTAGTGATCCTGTGGCTTCAGTTGAGCTGGGTTTGGAGCCAACAGAAGGAGGTGGAGCAGAATTCTGGACCCCTCAGTGTTCCAGAGGGAGCCATTGCCTCTCTCAACTGCACTTACAGTGACCGAGGTTCCCAGTCCTTCTTCTGGTACAGACAATATTCTGGGAAAAGCCCTGAGTTGATAATGTTCATATACTCCAATGGTGACAAAGAAGATGGAAGGTTTACAGCACAGCTCAATAAAGCCAGCCAGTATGTTTCTCTGCTCATCAGAGACTCCCAGCCCAGTGATTCAGCCACCTACCTCTGTGCCGTGAACTTCGGCGGAGGAAAGCTTATCTTCGGACAGGGAACGGAGTTATCTGTGAAACCCAATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTATATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCATGTGCAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTGTGATGTCAAGCTGGTCGAGAAAAGCTTTGAAACAGATACGAACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCCTCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCCAGC"         
        enzymes = ["BamHI", "SalI"]
        altered = TestRestrict.process(str, enzymes)
        unaltered = fxn.translate_nt(str)
        assert altered == unaltered

    def test_oneSite(self):
        """
        Checks to see if two sequences are identical w/ and w/out editing with our process and one restriction site
        """
        str = "ATGAAATCCTTGAGAGTTTTACTAGTGATCCTGTGGCTTCAGTTGAGCTGGGTTTGGAGCCAACAGAAGGAGGTGGAGCAGAATTCTGGACCCCTCAGTGTTCCAGAGGGAGCCATTGCCTCTCTCAACTGCACTTACAGTGACCGAGGTTCCCAGTCCTTCTTCTGGTACAGACAATATTCTGGGAAAAGCCCTGAGTTGATAATGTTCATATACTCCAATGGTGACAAAGAAGATGGAAGGTTTACAGCACAGCTCAATAAAGCCAGCCAGTATGTTTCTCTGCTCATCAGAGACTCCCAGCCCAGTGATTCAGCCACCTACCTCTGTGCCGTGAACTTCGGCGGAGGAAAGCTTATCTTCGGACAGGGAACGGAGTTATCTGTGAAACCCAATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACTCTAAATCCAGTGACAAGTCTGTCTGCCTATTCACCGATTTTGATTCTCAAACAAATGTGTCACAAAGTAAGGATTCTGATGTGTATATCACAGACAAAACTGTGCTAGACATGAGGTCTATGGACTTCAAGAGCAACAGTGCTGTGGCCTGGAGCAACAAATCTGACTTTGCATGTGCAAACGCCTTCAACAACAGCATTATTCCAGAAGACACCTTCTTCCCCAGCCCAGAAAGTTCCTGTGATGTCAAGCTGGTCGAGAAAAGCTTTGAAACAGATACGAACCTAAACTTTCAAAACCTGTCAGTGATTGGGTTCCGAATCCTCCTCCTGAAAGTGGCCGGGTTTAATCTGCTCATGACGCTGCGGCTGTGGTCCAGC"         
        enzymes = ["EcoRI"]
        altered = TestRestrict.process(str, enzymes)
        unaltered = fxn.translate_nt(str)
        assert altered == unaltered
