import seqtool.nucleotide.cpg
from nose.tools import *

def test_cpg():
    eq_(seqtool.nucleotide.cpg.count_cpg('ATGCGC'), 1)
    
