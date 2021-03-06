from .sw_c import smith_waterman
from collections import defaultdict

#from enum import Enum

__all__ = ['make_alignment']


from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def reverse_complement(seq_gap):
    return str(Seq(seq_gap, generic_dna).reverse_complement())

SCORE_THRESHOLD = 1.5

class GapBase:
    def __init__(self, index_gap, index_nongap=None, gapindex=None, gaplen=None, gapleft=None):
        self.index_gap = index_gap
        self.index_nongap = index_nongap
        self.gapindex = gapindex
        self.gaplen = gaplen
        self.gapleft = gapleft

    @property
    def gapright(self):
        return self.gapleft + 1 if self.isgap else None
    @property
    def isgap(self):
        return self.index_nongap is None
    @property
    def gaplerp(self):
        return (self.gapindex+1.)/(self.gaplen+1.) if self.isgap else None

class GapSeq:
    def __init__(self, seq_gap):
        """
        >>> g0 = GapSeq("A--TAT-GCC")
        >>> g0.gap_num
        3
        >>> len(g0)
        10
        >>> g0.seq_gap
        'A--TAT-GCC'
        >>> g0.seq_nongap
        'ATATGCC'
        >>> [x.index_gap for x in g0._index]
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        >>> [x.index_nongap for x in g0._index]
        [0, None, None, 1, 2, 3, None, 4, 5, 6]
        >>> [x.gapindex for x in g0._index]
        [None, 0, 1, None, None, None, 0, None, None, None]
        >>> [x.gaplen for x in g0._index]
        [None, 2, 2, None, None, None, 1, None, None, None]
        >>> [x.gapleft for x in g0._index]
        [None, 0, 0, None, None, None, 3, None, None, None]
        >>> g1 = GapSeq("ATGT-TGG-C")
        """
        assert(isinstance(seq_gap, str))

        self.seq_gap = seq_gap
        self.seq_nongap = seq_gap.replace('-','')
        
        self._index = self.calc_gap_index(seq_gap)

        self.gap_num = seq_gap.count('-')

        self._seq2mid = {}
        for gb in self._index:
            if not gb.isgap:
                self._seq2mid[gb.index_nongap] = gb.index_gap

    def __str__(self):
        return self.seq_gap

    def __repr__(self):
        return 'GapSeq("{}")'.format(self.seq_gap)

    def __len__(self):
        return len(self.seq_gap)

    def __getitem__(self, index_gap):
        return self.get_gap(index_gap)

    def get_gap(self, index_gap):
        return self._index[index_gap]

    def get_nongap(self, index_nongap):
        index_gap = self._seq2mid[index_nongap]
        return self[index_gap]
        
    @classmethod
    def calc_gap_index(cls, seq_gap):
        """
            seq             ATGC
            seq_gap          A--T---G-C
            index_gap       0123456789
            index_nongap    0  1   2 3
            gapindex         01 012 1
            gaplen           22 333 1
            gapleft          00 111 2
        """
        l = len(seq_gap)
        ret = [GapBase(x) for x in range(l)]

        index_nongap_i = 0
        gapindex_i = 0
        
        gap_flag = False
        gap_from = 0

        def gap_start(i):
            nonlocal gap_flag, gapindex_i, gap_from
            gap_flag = True
            gapindex_i = 0
            gap_from = i
            
        def gap_end(i):
            nonlocal gap_flag, gapindex_i, gap_from
            gaplen = i - gap_from
            for gb in ret[gap_from:i]:
                gb.gaplen = gaplen
            gap_flag = False
        
        for i,seq_gapi in enumerate(seq_gap):
            if seq_gapi == '-':
                if not gap_flag:
                    gap_start(i)
                    
                ret[i].gapindex = gapindex_i
                gapindex_i += 1
                ret[i].gapleft = index_nongap_i - 1
            else:
                if gap_flag:
                    gap_end(i)

                ret[i].index_nongap = index_nongap_i
                index_nongap_i += 1
        if gap_flag:
            gap_end(l)

        assert(not gap_flag)

        return ret

class AlignedSeq(object):
    def __init__(self, seq, first, last, mid_gap):
        """
        ATG--CAT-CATGC
          <------>
        >>> a = AlignedSeq("ATGCATCATGC",2,7,'G--CAT-C')
        >>> a.first
        'AT'
        >>> a.mid_gap.seq_nongap
        'GCATC'
        >>> a.mid_gap.seq_gap
        'G--CAT-C'
        >>> a.last
        'ATGC'
        >>> a.location
        (2, 7)
        >>> a.combined_gap
        'ATG--CAT-CATGC'
        >>> a.combined_nongap
        'ATGCATCATGC'
        >>> a.local(3,3)
        ' ATG--CAT-CA..'
        >>> a.local(10,10)
        '        ATG--CAT-CATGC      '
        """
        self.seq = seq
        self.location = (first,last)
        self.first = seq[:first]
        #self.mid = seq[first:last]
        self.last = seq[last:]

        self.mid_gap = GapSeq(mid_gap)

        self.combined_gap = self.first + self.mid_gap.seq_gap + self.last
        self.combined_nongap = self.seq

        self.adjust = len(self.first)

    def __repr__(self):
        return 'AlignedSeq("{}")'.format(self.combined_gap)

    def reverse_complement(self):
        seq = reverse_complement(self.seq)
        mid_gap = reverse_complement(self.mid_gap.seq_gap)
        l = len(self.seq)
        p,q = self.location
        
        return AlignedSeq(seq, l-q, l-p, mid_gap)

    def local(self, upstream, downstream):
        return self.first_adjust(upstream) + self.mid_gap.seq_gap + self.last_adjust(downstream)

    def first_adjust(self, length):
        fl = len(self.first)
        if fl <= length:
            return self.first.rjust(length)
        else:
            cc = min(2, length)
            return '.'*cc + self.first[fl-(length-cc):]
            
    def last_adjust(self, length):
        ll = len(self.last)
        if ll <= length:
            return self.last.ljust(length)
        else:
            cc = min(2, length)
            return self.last[:(length-cc)] + '.'*cc

class BaseMatch:
    MATCH = 0
    MISMATCH = 1
    NOTBAD = 2
    GAP = 3

    @classmethod
    def does_match(cls, i,j):
        if i == j:
            return cls.MATCH
        elif i=='-' or j=='-':
            return cls.GAP
        elif i=='N' or j=='N':
            return cls.NOTBAD

        elif i == 'Y':
            if j=='C' or j=='T':
                return cls.MATCH
        elif i == 'R':
            if j=='G' or j=='A':
                return cls.MATCH

        elif j == 'Y':
            if i=='C' or i=='T':
                return cls.MATCH
        elif j == 'R':
            if i=='G' or i=='A':
                return cls.MATCH
        return cls.MISMATCH

    @classmethod
    def match_bar(cls, s0, s1):
        """
        ATNCATGC
        ||:||| |
        ATGCATCC

        >>> match_bar('ATNCATGC','ATGCATCC')
        '||:||| |'
        """
        assert(len(s0)==len(s1))
        return ''.join('| : '[cls.does_match(s0[i],s1[i])] for i in range(len(s0)))

def dividing_point(left, right, ratio):
    """
    >>> dividing_point(1, 11, 0.5)
    6.0
    >>> dividing_point(1, 11, 0.2)
    3.0
    """
    return left + (right-left)*ratio

def make_alignment(seq_s0,seq_s1):
    s, first, last, mid, mv = smith_waterman(str(seq_s0),str(seq_s1))
    aseq0 = AlignedSeq(s[0], first[0], last[0], mid[0])
    aseq1 = AlignedSeq(s[1], first[1], last[1], mid[1])
    return Alignment(aseq0, aseq1, mv)

class CorrespondanceMap:
    def __init__(self, s0, s1):
        self.s0 = s0
        self.s1 = s1
        assert(len(s0)==len(s1))
        self.m_map = defaultdict(lambda: defaultdict(int))
        self.m_str = defaultdict(str)

        for p,q in zip(s0,s1):
            self.m_map[p][q] += 1
            self.m_str[p] += q

    def get_bsa_result(self):
        y = self.m_str['Y']
        g = self.m_str['R']
        return y or g
            
    def text_map(self):
        ret = ''
        for k,v in self.m_map.items():
            ret += str(k)+':{'+', '.join('{}:{}'.format(kk,vv) for kk,vv in list(v.items()))+'}'
        return ret
        
    def text_str(self):
        ret = ''
        for k,v in self.m_str.items():
            ret += str(k)+': '+v+'\n'
        return ret
    
    
class Alignment(object):
    def __init__(self, aseq0, aseq1, mv):
        assert(len(aseq0.mid_gap)==len(aseq1.mid_gap))
        self.aseq0 = aseq0
        self.aseq1 = aseq1
        self.score = mv

    def __repr__(self):
        return 'Alignment(score={}, seq0={}, seq1={})'.format(self.score, self.aseq0, self.aseq1)
        
    def reversed(self):
        return Alignment(self.aseq1, self.aseq0, self.score)

    def reverse_complement(self):
        return Alignment(self.aseq0.reverse_complement(), self.aseq1.reverse_complement(), self.score)

    def text_all(self):
        u,d = self.get_longest_first_last_length()
        return '{}\n{}\n{}'.format(self.aseq1.local(u,d),
                                   ' '*u + self.match_bar() + ' '*d,
                                   self.aseq0.local(u,d) )

    def text_shrinked(self, upstream=30, downstream=30):
        u,d = self.get_longest_first_last_length()
        u = min(upstream, u)
        d = min(downstream, d)
        
        return '{}\n{}\n{}'.format(self.aseq1.local(u,d),
                                   ' '*u+self.match_bar()+' '*d,
                                   self.aseq0.local(u,d) )

    def match_bar(self):
        return BaseMatch.match_bar(self.aseq0.mid_gap.seq_gap, self.aseq1.mid_gap.seq_gap)

    def correspondance_map(self):
        u,d = len(self.aseq0.first),len(self.aseq0.last)
        s0 = self.aseq0.local(u,d)
        s1 = self.aseq1.local(u,d)
        assert(len(self.aseq0.mid_gap)==len(self.aseq1.mid_gap))
        return CorrespondanceMap(s0,s1)

    def get_common_first_last_length(self):
        u = min(len(self.aseq0.first),len(self.aseq1.first))
        d = min(len(self.aseq0.last),len(self.aseq1.last))
        return u,d

    def get_longest_first_last_length(self):
        u = max(len(self.aseq0.first),len(self.aseq1.first))
        d = max(len(self.aseq0.last),len(self.aseq1.last))
        return u,d
        
    def get_loc(self, seq0_location, upstream=0, downstream=0):
        u,d = upstream,downstream
        p,q = self.aseq0.location

        yield from seq0_location[p-u:p]
        
        loc = seq0_location[p:q]
        refe = self.aseq0.mid_gap
        targ = self.aseq1.mid_gap

        for i,base in enumerate(targ.seq_gap):
            gapbase = refe[i]
            if gapbase.isgap:
                yield dividing_point(loc[gapbase.gapleft], loc[gapbase.gapright], gapbase.gaplerp)
            else:
                yield loc[gapbase.index_nongap]

        yield from seq0_location[q:q+d]

    def get_match_result(self):
        u,d = len(self.aseq0.first),len(self.aseq0.last)
        compare,gap = self.compare_bar()
        return '-'*u + compare + '-'*d
        
    def compare_bar(self):
        """
        seq1   : ATTTG-CA
        match  : |   |  :
        seq0   : A---AACN

        seq1_compare: Ag-CA
        gap         :  <
        seq0        : AAACN
        """

        base = self.aseq0.mid_gap.seq_gap
        view = self.aseq1.mid_gap.seq_gap

        compare = []
        gap = []

        next_gap = False
        def append_gap():
            nonlocal gap, next_gap
            gap.append('<' if next_gap else '_')
            next_gap = False

        for bi,vi in zip(base,view):
            m = BaseMatch.does_match(bi,vi)
            if m==BaseMatch.MATCH or m==BaseMatch.NOTBAD:
                compare.append(vi)
                append_gap()
            elif m==BaseMatch.GAP and bi=='-':
                next_gap = True
            elif m==BaseMatch.GAP and vi=='-':
                compare.append('-')
                append_gap()
            elif m==BaseMatch.MISMATCH:
                compare.append(vi.lower())
                append_gap()
            else:
                raise Exception('unreachable')

        assert(len(compare)==len(gap))

        return ''.join(compare), ''.join(gap)

    def length(self):
        return len(self.aseq0.mid_gap)

    def score_density(self):
        return 1. * self.score / self.length()

    def criteria(self):
        return self.score_density() > SCORE_THRESHOLD

    def score_text(self):
        return '{} = {:.2f} * {}'.format(self.score, self.score_density(), self.length())
                         
if __name__=='__main__':
    print(alignment('ACACACTA','AGCACACA').text_local())
    print(alignment('GCCCTAGCG', 'GCGCAATG').text_local())
