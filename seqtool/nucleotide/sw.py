from .sw_c import smith_waterman

def match_nucleotide(a,b):
    return a==b

def match_bar(s0, s1):
    assert(len(s0)==len(s1))
    return ''.join('|' if match_nucleotide(s0[i],s1[i]) else ' ' for i in xrange(len(s0)))

class Alignment(object):
    def __init__(self, seq_s0, seq_s1):
        self.aseq0, self.aseq1 = smith_waterman(seq_s0,seq_s1)

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
