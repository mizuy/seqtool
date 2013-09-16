from seqtool.nucleotide.alignment import make_alignment
from nose.tools import *
from nose.tools import set_trace

seq1 = 'AGACGGAGTTTGTGAGTGGTTTTTGGTCGGAGGGACGGGGTGGGTTGAGT'
seq0 = 'GTTTTTGGTYGGGGAAGGAYGGGGTGGGTGAGTYGTGYGTTTTTTYGGGYG'

seq1_mid_loc = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 27.5, 28, 29, 30, 31, 32]
seq0_mid_loc = [18+x for x in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9.333333333333332, 9.666666666666668, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]]

text_all = ['AGACGGAGTTTGTGAGTGGTTTTTGGTC--GGAGGGACGGGGTGGGTTGAGT                  ',
            '                  ||||||||||  ||| |||||||||||| |||||                  ',
            '                  GTTTTTGGTYGGGGAAGGAYGGGGTGGG-TGAGTYGTGYGTTTTTTYGGGYG']
text_shrinked = ['..GTGAGTGGTTTTTGGTC--GGAGGGACGGGGTGGGTTGAGT        ',
                 '         ||||||||||  ||| |||||||||||| |||||        ',
                 '         GTTTTTGGTYGGGGAAGGAYGGGGTGGG-TGAGTYGTGYG..']

def test_alignment():
    a = make_alignment(seq0, seq1)
    
    eq_(a.score, 49)
    
    eq_([x for x in a.text_all().splitlines()], text_all)
    eq_([x for x in a.text_shrinked(9,8).splitlines()], text_shrinked)

    eq_(list(a.get_loc([x for x in range(len(seq0))])),seq1_mid_loc)
    eq_(list(a.reversed().get_loc([x for x in range(len(seq1))])),seq0_mid_loc)


