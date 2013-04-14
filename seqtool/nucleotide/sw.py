from .sw_c import smith_waterman

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

def bar_nucleotide(i,j):
    m = does_match(i,j)
    if m==M_MATCH:
        return '|'
    elif m==M_MISMATCH:
        return ' '
    elif m==M_NOTBAD:
        return ':'
    elif m==M_GAP:
        return ' '
    return ' '

def match_bar(s0, s1):
    assert(len(s0)==len(s1))
    return ''.join(bar_nucleotide(s0[i],s1[i]) for i in xrange(len(s0)))

class Alignment(object):
    def __init__(self, seq_s0, seq_s1):
        self.aseq0, self.aseq1, self.score = smith_waterman(str(seq_s0),str(seq_s1))

    def text_all(self):
        p,q = self.aseq0.adjust, self.aseq1.adjust
        r = max(p,q)
        return '{}\n{}'.format(' '*(r-p) + self.aseq0.combined(), 
                               ' '*(r-q) + self.aseq1.combined() )

    def text_local(self, upstream=30, downstream=30):
        p,q = self.aseq0.adjust, self.aseq1.adjust
        return '{}\n{}\n{}'.format(self.aseq0.local(upstream,downstream), 
                               ' '*upstream+self.match_bar(),
                               self.aseq1.local(upstream,downstream) )

    def match_bar(self):
        return match_bar(self.aseq0.mid_gap, self.aseq1.mid_gap)

    def compare_bar(self, index):
        """
        view   : ATTTG-CA
        mb     : |   |:
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

        compare = ''
        gap = ''
        next_gap = True

        length = len(base)
        for i in range(length):
            bi = base[i]
            vi = view[i]
            m = does_match(bi,vi)
            if m==M_MATCH or m==M_NOTBAD:
                compare += vi
                gap += ('<' if next_gap else '-')
                next_gap = False
            elif m==M_GAP:
                if bi=='-':
                    next_gap = True
                else:
                    compare += '-'
                    gap += ('<' if next_gap else '-')
                    next_gap = False
            elif m==M_MISMATCH:
                compare += vi.lower()
                gap += ('<' if next_gap else '-')
                next_gap = False
            else:
                raise Exception('unreachable')

        assert(len(compare)==len(gap)==len(based))

        return compare, gap

    def length(self):
        return len(self.aseq0.mid)

    def score_density(self):
        return 1. * self.score / self.length()

    def score_text(self):
        return '{} = {:.2f} * {}'.format(self.score, self.score_density(), self.length())


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
