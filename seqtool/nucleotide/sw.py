from .sw_c import smith_waterman

def bar_nucleotide(i,j):
    if i == j:
        return '|'
    elif i=='N' or j=='N':
        return ':'

    elif i == 'Y':
        if j=='C' or j=='T':
            return '|'
    elif i == 'R':
        if j=='G' or j=='A':
            return '|'

    elif j == 'Y':
        if i=='C' or i=='T':
            return '|'
    elif j == 'R':
        if i=='G' or i=='A':
            return '|'

    return ' '

def match_bar(s0, s1):
    assert(len(s0)==len(s1))
    return ''.join(bar_nucleotide(s0[i],s1[i]) for i in xrange(len(s0)))

class Alignment(object):
    def __init__(self, seq_s0, seq_s1):
        self.aseq0, self.aseq1, self.score = smith_waterman(seq_s0,seq_s1)

    def text_all(self):
        p,q = self.aseq0.adjust, self.aseq1.adjust
        r = max(p,q)
        return '{}\n{}'.format(' '*(r-p) + self.aseq0.combined(), 
                               ' '*(r-q) + self.aseq1.combined() )

    def text_local(self, upstream=30, downstream=30):
        p,q = self.aseq0.adjust, self.aseq1.adjust
        return '{}\n{}\n{}'.format(self.aseq0.local(upstream,downstream), 
                               ' '*upstream+match_bar(self.aseq0.mid, self.aseq1.mid),
                               self.aseq1.local(upstream,downstream) )

    def length(self):
        return len(self.aseq0.mid)

    def score_density(self):
        return 1. * self.score / self.length()

def print_sw(s0, s1):
    a = Alignment(s0, s1)

    print 'alignment of :'
    print a.aseq0.seq
    print a.aseq1.seq
    print 'is:'
    print a.text()

if __name__=='__main__':
    print_sw('ACACACTA','AGCACACA')
    print_sw('GCCCTAGCG', 'GCGCAATG')
