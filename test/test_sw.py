import seqtool.nucleotide.sw as sw
from nose.tools import *
from nose.tools import set_trace

seq1 = 'AGACGGAGTTTGTGAGTGGTTTTTGGTCGGGAGGGACGGGGTGGGTTGAGT'
seq0 = 'GTTTTTGGTYGGGGAAGGAYGGGGTGGGTGAGTYGTGYGTTTTTTYGGGYG'

seq1_loc = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,     11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 27.5, 28, 29, 30, 31, 32]
seq0_loc = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,       28, 29, 30, 31, 32]

text_local = ['            AGACGGAGTTTGTGAGTGGTTTTTGGTC-GGGAGGGACGGGGTGGGTTGAGT',
              '                              |||||||||| |||| |||||||||||| |||||',
              '                              GTTTTTGGTYGGGGAAGGAYGGGGTGGG-TGAGTYGTGYGTTTTTTYGGGYG']


def test_alignment():
    a = sw.Alignment(seq0, seq1)
    
    eq_(a.score, 54)
    eq_(a.correspondance_map()['G']['G'], 17)
    eq_(a.correspondance_str()['G'], 'GGGGGGAGGGGGGGGGGG')
    
    eq_([x.rstrip() for x in a.text_local().splitlines()], text_local)

    eq_(a.get_loc(seq1_loc,index=1),seq0_loc)
    eq_(a.get_loc(seq0_loc,index=0),seq1_loc)


