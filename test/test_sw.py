import seqtool.nucleotide.sw as sw
from nose.tools import *
from nose.tools import set_trace

seq1 = 'AGACGGAGTTTGTGAGTGGTTTTTGGTCGGAGGGACGGGGTGGGTTGAGT'
seq0 = 'GTTTTTGGTYGGGGAAGGAYGGGGTGGGTGAGTYGTGYGTTTTTTYGGGYG'

seq1_mid_loc = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,         12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 27.5, 28, 29, 30, 31, 32]
seq1_loc     = [0 for x in range(18)] + seq1_mid_loc
seq0_mid_loc = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,       28, 29, 30, 31, 32]
seq0_loc     = seq0_mid_loc + [33+x for x in range(18)]

text_local = ['            AGACGGAGTTTGTGAGTGGTTTTTGGTC--GGAGGGACGGGGTGGGTTGAGT',
              '                              ||||||||||  ||| |||||||||||| |||||',
              '                              GTTTTTGGTYGGGGAAGGAYGGGGTGGG-TGAGTYGTGYGTTTTTTYGGGYG']


def test_alignment():
    a = sw.Alignment(seq0, seq1)
    
    eq_(a.score, 49)
    eq_(a.correspondance_map()['G']['G'], 16)
    eq_(a.correspondance_str()['G'], 'GGGGGAGGGGGGGGGGG')
    
    eq_([x.rstrip() for x in a.text_local().splitlines()], text_local)

    eq_(a.get_mid_loc(seq1_loc,index=1),seq0_mid_loc)
    eq_(a.get_mid_loc(seq0_loc,index=0),seq1_mid_loc)


