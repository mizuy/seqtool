from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

from collections import defaultdict

from . import xmlwriter

def cpg_obs_per_exp(seq):
    """
    CpG islands in vertebrate genomes {Gardiner-Garden, 1987, p02206}
    'Obs/Exp CpG' = N * 'Number of CpG' / 'Number of C' * 'Number of G'
    where, N is the total number of nucleotide in the sequence being analyzed.
    """
    n = len(seq)
    c = 0
    g = 0
    cpg = 0
    for i,b in enumerate(seq):
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
