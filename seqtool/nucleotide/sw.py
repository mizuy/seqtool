from .sw_c import smith_waterman
from collections import defaultdict

__all__ = ['Alignment']


class GappedSequence:
    def __init__(self, gapseq):
        assert(isinstance(gapseq, str))
        self.gapseq = gapseq
        self.length = len(gapseq)
        
        self.gap_num = gapseq.count('-')
        self.seq = gapseq.replace('-','')
        self.index = list(self.calc_gap_index(gapseq))

        self.seq2mid = {}
        for i,(s,g,gl,gapleft) in enumerate(self.index):
            if s is not None:
                self.seq2mid[s] = i

    def __str__(self):
        return self.gapseq
    def __len__(self):
        return len(self.gapseq)

    def get_seqindex(self,i):
        return self.index[0][i]
        
    @classmethod
    def calc_gap_index(cls, midgap):
        """
            midgap      A--A---A-A
            midindex    0123456789
            seqindex    0  1   2 3
            gapindex     01 012 1
            gaplen       22 333 1
            gapleft      00 111 2
            (gapright)   11 222 3
        
        
        >>> calc_gap_index("A--A---A-A")
        [(0,0), (1,2), (2,2), (0,0), (1,3), (2,3), (3,3), (0,0), (1,1), (0,0)]
        """
        l = len(midgap)
        gap_num = midgap.count('-')

        seqindex = [None for x in range(l)]
        gapindex = [None for x in range(l)]
        gaplen = [None for x in range(l)]
        gapleft = [None for x in range(l)]

        seqindex_i = 0
        gapindex_i = 0
        
        gap_flag = False
        gap_from = 0

        for i in range(l):
            if midgap[i]=='-':
                if gap_flag:
                    pass
                else:
                    # gap start
                    gap_flag = True
                    gapindex_i = 0
                    gap_from = i
                    
                gapindex[i] = gapindex_i
                gapindex_i += 1
            else:
                if gap_flag:
                    # gap end
                    ll = i - gap_from
                    gaplen[gap_from:i] = [ll for i in range(ll)]
                    gap_flag = False
                else:
                    pass

                seqindex[i] = seqindex_i
                seqindex_i += 1

        seqindex_i = 0
        for i in range(l):
            if midgap[i]=='-':
                pass
            else:
                seqindex_i += 1
            gapleft[i] = seqindex_i

        return zip(seqindex, gapindex, gaplen, gapleft)

class AlingnedSeq(object):
    def __init__(self, seq, first, last, mid_gap):
        self.seq = seq
        self.location = (first,last)
        self.first = seq[:first]
        self.mid = seq[first:last]
        self.last = seq[last:]

        self.mid_gap = GappedSequence(mid_gap)

        self.adjust = len(self.first)

    def combined(self):
        return self.first + str(self.mid_gap) + self.last

    def local(self, upstream, downstream):
        return self.first[-upstream:].rjust(upstream) + str(self.mid_gap) + self.last[:downstream].ljust(downstream)

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
    return ''.join('| : '[does_match(s0[i],s1[i])] for i in range(len(s0)))

def dividing_point(left, right, ratio):
    return left + (right-left)*ratio
    
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
        return match_bar(str(self.aseq0.mid_gap), str(self.aseq1.mid_gap))

    def correspondance_map(self):
        m = defaultdict(lambda: defaultdict(int))

        length = len(str(self.aseq0.mid_gap))
        for i in range(length):
            p = str(self.aseq0.mid_gap)[i]
            q = str(self.aseq1.mid_gap)[i]
            m[q][p] += 1

        return m

    def correspondance_str(self):
        m = defaultdict(str)

        length = len(str(self.aseq0.mid_gap))
        for i in range(length):
            p = str(self.aseq0.mid_gap)[i]
            q = str(self.aseq1.mid_gap)[i]
            m[q] += p

        return m

    def get_mid_loc(self, refe_location):
        p,q = self.aseq0.location
        loc = refe_location[p:q]
        refe = self.aseq0.mid_gap
        targ = self.aseq1.mid_gap
        for i,base in enumerate(targ.gapseq):
            seqindex,gapindex,gaplen,gapleft = refe.index[i]
            if seqindex is not None:
                yield loc[seqindex]
            else:
                # gap
                ratio = (gapindex + 1) / (gaplen + 1)
                yield dividing_point(loc[gapleft], loc[gapleft+1], ratio)
                
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
            base = str(self.aseq1.mid_gap)
            based = self.aseq1.mid
            view = str(self.aseq0.mid_gap)
        else:
            base = str(self.aseq0.mid_gap)
            based = self.aseq0.mid
            view = str(self.aseq1.mid_gap)

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
    print(Alignment('ACACACTA','AGCACACA').text_local())
    print(Alignment('GCCCTAGCG', 'GCGCAATG').text_local())
