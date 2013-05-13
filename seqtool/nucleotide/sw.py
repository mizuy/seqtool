from .sw_c import smith_waterman
from collections import defaultdict

__all__ = ['Alignment']

class AlingnedSeq(object):
    def __init__(self, seq, first, last, mid_gap):
        self.seq = seq
        self.location = (first,last)
        self.first = seq[:first]
        self.mid = seq[first:last]
        self.last = seq[:last]

        self.mid_gap = mid_gap

        self.adjust = len(self.first)

    def combined(self):
        return self.first + self.mid_gap + self.last

    def local(self, upstream, downstream):
        return self.first[-upstream:].rjust(upstream) + self.mid_gap + self.last[:downstream].ljust(downstream)

M_MATCH = 0
M_MISMATCH = 1
M_NOTBAD = 2
M_GAP = 3

def does_match(i,j):
    if i == j:
        return M_MATCH
    elif i=='-' or j=='-':
        return M_GAP
    elif i=='N' or j=='N':
        return M_NOTBAD

    elif i == 'Y':
        if j=='C' or j=='T':
            return M_MATCH
    elif i == 'R':
        if j=='G' or j=='A':
            return M_MATCH

    elif j == 'Y':
        if i=='C' or i=='T':
            return M_MATCH
    elif j == 'R':
        if i=='G' or i=='A':
            return M_MATCH
    return M_MISMATCH

def match_bar(s0, s1):
    assert(len(s0)==len(s1))
    return ''.join('| : '[does_match(s0[i],s1[i])] for i in xrange(len(s0)))

class Alignment(object):
    def __init__(self, seq_s0, seq_s1):
        s, first, last, mid, mv = smith_waterman(str(seq_s0),str(seq_s1))
        self.aseq0 = AlingnedSeq(s[0], first[0], last[0], mid[0])
        self.aseq1 = AlingnedSeq(s[1], first[1], last[1], mid[1])
        self.score = mv

    def text_all(self):
        p,q = self.aseq0.adjust, self.aseq1.adjust
        r = max(p,q)
        return '{}\n{}'.format(' '*(r-p) + self.aseq0.combined(), 
                               ' '*(r-q) + self.aseq1.combined() )

    def text_local(self, upstream=30, downstream=30):
        p,q = self.aseq0.adjust, self.aseq1.adjust
        return '{}\n{}\n{}'.format(self.aseq1.local(upstream,downstream), 
                               ' '*upstream+self.match_bar(),
                               self.aseq0.local(upstream,downstream) )

    def match_bar(self):
        return match_bar(self.aseq0.mid_gap, self.aseq1.mid_gap)

    def correspondance_map(self):
        m = defaultdict(lambda: defaultdict(int))

        length = len(self.aseq0.mid_gap)
        for i in range(length):
            p = self.aseq0.mid_gap[i]
            q = self.aseq1.mid_gap[i]
            m[q][p] += 1

        return m

    def correspondance_str(self):
        m = defaultdict(str)

        length = len(self.aseq0.mid_gap)
        for i in range(length):
            p = self.aseq0.mid_gap[i]
            q = self.aseq1.mid_gap[i]
            m[q] += p

        return m

    def compare_bar(self, index):
        """
        view   : ATTTG-CA
        match  : |   |  :
        base   : A---AACN

        compare: Ag-CA
        gap    :  <   
        base'  : AAACN
        """
        if index == 0:
            base = self.aseq1.mid_gap
            based = self.aseq1.mid
            view = self.aseq0.mid_gap
        else:
            base = self.aseq0.mid_gap
            based = self.aseq0.mid
            view = self.aseq1.mid_gap

        compare = []
        gap = []
        next_gap = False

        length = len(base)
        for i in range(length):
            bi = base[i]
            vi = view[i]
            m = does_match(bi,vi)
            if m==M_MATCH or m==M_NOTBAD:
                compare.append(vi)
                gap.append('<' if next_gap else '_')
                next_gap = False
            elif m==M_GAP:
                if bi=='-':
                    next_gap = True
                else:
                    compare.append('-')
                    gap.append('<' if next_gap else '_')
                    next_gap = False
            elif m==M_MISMATCH:
                compare.append(vi.lower())
                gap.append('<' if next_gap else '_')
                next_gap = False
            else:
                raise Exception('unreachable')

        assert(len(compare)==len(gap)==len(based))

        return ''.join(compare), ''.join(gap)

    def length(self):
        return len(self.aseq0.mid)

    def score_density(self):
        return 1. * self.score / self.length()

    def score_text(self):
        return '{} = {:.2f} * {}'.format(self.score, self.score_density(), self.length())

if __name__=='__main__':
    print Alignment('ACACACTA','AGCACACA').text_local()
    print Alignment('GCCCTAGCG', 'GCGCAATG').text_local()
