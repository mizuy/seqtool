from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from collections import defaultdict
import re
import nucleotide,htmlwriter

class RestrictionEnzyme(object):
    def __init__(self, pattern,sense_cut,antisense_cut):
        pass
    def digest(self, dsdna):
        pass

class Enzymes(object):
    XHOI =   RestrictionEnzyme(pattern='C+TCGA+G',sense_cut=1,antisense_cut=-1)
    BAMHI =  RestrictionEnzyme(pattern='G^GATC^C',sense_cut=1,antisense_cut=-1)
    ECHORI = RestrictionEnzyme(pattern='GAATTC',sense_cut=1, antisense_cut=-1)

class RestrictionEnzymeMap(object):
    def __init__(self, enzymes):
        pass
    def list_all(self, dsdna):
        pass
    def list_unique(self, dsdna):
        pass

class dsDNA(object):
    def __init__(self, sequence, left_cut=0, right_cut=0):
        self._seq = Seq.Seq(str(sequence))
        self._cseq = self._seq.complement()
        common_len = len(sequence)-abs(left_cut)+abs(right_cut)
        if common_len <= 0:
            raise ValueError('the dsDNA fragment has no common region between sense and antisense fragments: %s'%common_len)

        self.left = left_cut
        self.right = right_cut
        
    def __len__(self):
        return len(self._seq)

    def _repr_cut(self, seq, left, right):
        assert left>=0 and right>=0
        return ' '*left + str(seq[left:len(seq)-right]) + ' '*right

    def to_string(self):
        sense = self._repr_cut(self._seq,max(0,self.left),max(0,self.right))
        antisense = self._repr_cut(self._cseq, max(0,-self.left), max(0,-self.right))
        return (sense,antisense)
    def print_strand(self):
        (s,a) = self.to_string()
        print "5'-%s-3'"%s
        print "3'-%s-5'"%a

    def __repr__(self):
        def sign(v):
            if v>0:
                return +1
            elif v==0:
                return 0
            else:
                return -1
        ch = ['+','','|']
        lc = ch[sign(self.left)+1]
        rc = ch[sign(self.right)+1]
        l = abs(self.left)
        r = len(self)-abs(self.right)
        assert l<r
        j = self._seq[:l] + lc + self._seq[l:r] + rc + self._seq[r:]
        return 'dsDNA("%s")'%j
    def ligation(self, other):
        pass

class ssDNA(object):
    def __init__(self, sequence):
        pass

def ligation(fragments):
    pass

class DecorationStack(object):
    def __init__(self):
        self.s = []
        pass
    def push(self, sclass, color, bgcolor):
        self.s.append((sclass,color,bgcolor))
    def pop(self):
        self.s = self.s[:-1]
    def decoration(self):
        sclass = ' '.join([sclass for sclass,color,bgcolor in self.s])
        color = '#000000'
        bgcolor = '#ffffff'
        for sclass,f,b in self.s:
            if f:
                color = f
            if b:
                bgcolor = b
        return sclass, 'color:%s; background-color:%s;'%(color,bgcolor)

class AnnotatedSequence(object):
    def __init__(self, seq):
        self.seq = 'N'+seq
        self.anno_start = defaultdict(lambda:[])
        self.anno_end = defaultdict(lambda:[])
        self.anno_point = defaultdict(lambda:[])

    def add_region(self, index_from, index_to, name, color):
        '''
        ...ATGC(name)[background_colored_sequence]ATGC...
        '''
        self.anno_start[index_from].append((name,('region',None,color)))
        self.anno_end[index_to].append((name,('region',None,color)))
    def add_point(self, index, name):
        '''
        ...TGAC(name)TGAC...
        '''
        self.anno_point[index].append((name))
    def add_pattern(self, pattern, name, color):
        '''
        ...ATGC(name)[underlined_colored_sequence]ATGC...
        '''
        
        pp,pc = self._search(pattern)
        for f,t in pp:
            n = '>'+name if name else None
            self.anno_start[f].append((n,('pattern',color, None)))
            self.anno_end[t].append((n,('pattern',color,None)))
        for f,t in pc:
            n = '<'+name if name else None
            self.anno_start[f].append((n,('pattern',color,None)))
            self.anno_end[t].append((n,('pattern',color,None)))

    def _search(self, pattern):
        seqs = str(self.seq)
        primer = Seq.Seq(pattern)
        cprimer = primer.reverse_complement()
        reg = re.compile('(%s)|(%s)'%(nucleotide._gen_re_seq(primer),nucleotide._gen_re_seq(cprimer)))
        pp = []
        pc = []
        start = 0
        while True:
            m = reg.search(seqs, start)
            if not m:
                break
            start = m.start()+1
            if m.group(1):
                pp.append((m.start(),m.end()))
            if m.group(2):
                pc.append((m.start(),m.end()))
        return pp,pc

    def print_html(self, w):
        stack = DecorationStack()
        w.push('pre',style='word-wrap:break-word;')
        for i in range(1,len(self.seq)):
            for name,(sclass,color,bgcolor) in self.anno_start[i]:
                stack.push(sclass,color,bgcolor)
                sclass,style=stack.decoration()
                w.push('span',**{'class':sclass, 'style':style})
                if name:
                    w.text('(%s)'%name)

            for name in self.anno_point[i]:
                w.text('(%s)'%name)

            w.text(self.seq[i])

            for name,(sclass,color,bgcolor) in self.anno_end[i]:
                stack.pop()
                w.pop()
        w.pop()
