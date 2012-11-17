from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

from collections import defaultdict

from . import xmlwriter

def drop_small_region(items, length):
    ret = []
    for p,q in items:
        if q-p >= length:
            ret.append((p,q))
    return ret

def region_growing(items, gap=0):
    if len(items)<2:
        return items
    p,q = items[0]
    rest = items[1:]
    for i in range(len(rest)):
        r,s = rest[i]
        if r<q+gap:
            q = max(q,s)
            continue
        else:
            return [(p,q)]+region_growing(rest[i:],gap)

# NOTE: this cpg island implementation is ad-hoc.
class cpgisland_searcher(object):
    def __init__(self, length, window):
        self.length = length
        self.window = window
        self.in_region = False
        self.islands = []
        self.start = 0

    def p(self, i, in_region):
        if self.in_region and (not in_region):
            """region exit"""
            self.islands.append((self.start, i-1))
            self.in_region = False
        elif (not self.in_region) and in_region:
            """region enter"""
            self.start = i
            self.in_region = True

    def finish(self):
        self.p(self.length, False)
        h = self.window/2
        
        connected = region_growing([(p-h, q+h) for p,q in self.islands])
        # Two individual CpG islands were connected if they were separated by less than 100 bp
        gap_connected = region_growing(connected, 100)
        return drop_small_region(gap_connected, 500)

def seq_cpg_analysis(seq, window):
    """
    calculate gc percent for entire sequence.
    """
    seqstr = str(seq).upper()
    l = len(seq)
    h = int(window/2)
    gc_per = []
    obs = []

    sr = cpgisland_searcher(l,window)
    for i in xrange(0,l):
        p = max(0,i-h)
        q = min(i+h,l)
        n = q-p

        c = seqstr.count('C',p,q)
        g = seqstr.count('G',p,q)
        cg = seqstr.count('CG',p,q)

        gcp = 1.*(c+g)/n
        oe = 1.*n*cg/(c*g) if (c*g)!=0 else 0

        gc_per.append(gcp)
        obs.append(oe)
    
        island = (gcp > 0.55) and (oe > 0.65)
        sr.p(i, island)

    cpg_islands = sr.finish()

    return gc_per, obs, cpg_islands

def seq_gc_percent(seq, window):
    """
    calculate gc percent for entire sequence.
    """
    seqstr = str(seq).upper()
    l = len(seq)
    h = int(window/2)
    ret = []
    for i in xrange(0,l):
        p = max(0,i-h)
        q = min(i+h,l)
        n = q-p

        c = seqstr.count('C',p,q)
        g = seqstr.count('G',p,q)
        gcp = 1.*(c+g)/n

        ret.append(gcp)
    return ret

def seq_cpg_obs_per_exp(seq, window, step):
    """
    calculate obs/exp cpg for entire sequence.
    """
    seqstr = str(seq).upper()
    l = len(seq)
    h = int(window/2)
    ret = []
    for i in xrange(0,l,step):
        p = max(0,i-h)
        q = min(i+h,l)
        n = q-p

        c = seqstr.count('C',p,q)
        g = seqstr.count('G',p,q)
        cg = seqstr.count('CG',p,q)
        oe = 1.*n*cg/(c*g) if (c*g)!=0 else 0

        ret.append(oe)
    return ret

def cpg_obs_per_exp(seq):
    """
    CpG islands in vertebrate genomes {Gardiner-Garden, 1987, p02206}
    'Obs/Exp CpG' = N * 'Number of CpG' / 'Number of C' * 'Number of G'
    where, N is the total number of nucleotide in the sequence being analyzed.

    2 >= Obs/Exp >= 0
    """
    n = len(seq)
    c = 0
    g = 0
    cpg = 0
    for i,b in xenumerate(seq):
        b = b.upper()
        if b=='C': c+=1
        if b=='G': g+=1
        if i+1<n and b=='C' and seq[i+1]=='G':
            cpg += 1

    if c*g == 0:
        return 0
    return 1.*n*cpg / (c*g)


#def cpg_sites(seq): return [i for i in range(0,len(seq)-1) if seq[i]=='C' and seq[i+1]=='G']
def cpg_sites(seq):
    seqstr = str(seq)
    ret = []
    j = 0
    while 1:
        j = seqstr.find('CG', j)
        if j<0:
            break
        ret.append(j)
        j += 1
    return ret

def base_color(n):
    return {'A':'#00FF00',
            'T':'#FF0000',
            'G':'#000000',
            'C':'#0000FF'}[n.upper()]

def tm_gc(seq):
    gc = seq.count('G')+seq.count('C')
    at = seq.count('A')+seq.count('T')
    return '4x%s+2x%s=%s'%(gc,at,4*gc+2*at)

def is_cpg(seq,i):
    if i+1>=len(seq):
        return False
    if seq[i]=='C' and seq[i+1]=='G':
        return True
    return False

def is_repeat(seq,i,repeatno=8):
    v = seq[i]
    tail = 0
    for t in seq[i+1:]:
        if t!=v:
            break
        tail += 1
    head = 0
    for t in seq[:i][::-1]:
        if t!=v:
            break
        head += 1
    return tail+head+1 >= repeatno

def count_cpg(seq):
    count = 0
    for i,n in enumerate(seq[:-1]):
        if n=='C' and seq[i+1]=='G':
            count += 1
    return count

def bisulfite(seq, methyl):
    key = '_bisulfite_' + ('met' if methyl else 'unmet')
    if not hasattr(seq, key):
        muta = seq.tomutable()
        old = 'X'
        for i,c in enumerate(muta[:-2]):
            if c=='C':
                if not(muta[i+1]=='G' and methyl):
                    muta[i] = 'T'

        ret = muta.toseq()

        setattr(seq, key, ret)

    return getattr(seq, key)

class ColorMap(object):
    def __init__(self):
        self.colors = defaultdict(lambda: (0,0,0))

    def add_color(self, i, red, green, blue):
        r,g,b = self.colors[i]
        self.colors[i] = (r+red, g+green, b+blue)

    def get_color(self, i):
        return '#%02x%02x%02x'%self.colors[i]

def pprint_sequence_html(w, seq, get_color):
    b = xmlwriter.builder(w)
    seqstr = str(seq)
    with b.pre:
        for i in range(0, len(seqstr), 100):
            w.write('<span style="black">%4d: </span>'%i)
            oldcol = None
            for ii in range(i, min(i+100, len(seqstr)), 10):
                for iii in range(ii, min(ii+10, len(seqstr)), 1):
                    color = get_color(iii)
                    if oldcol!=color:
                        if oldcol:
                            w.write('</span>')
                        w.write('<span style="color:%s">'%color)
                        oldcol = color
                    w.write(seqstr[iii])
                w.write(' ')
            w.write('\n')
            if oldcol:
                w.write('</span>')
