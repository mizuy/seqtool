from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC

from collections import defaultdict

from ..util import xmlwriter

def to_seq(s):
    if isinstance(s, Seq.Seq):
        return s

    t = str(s).upper()
    if all(n in IUPAC.unambiguous_dna.letters for n in t):
        return Seq.Seq(t,IUPAC.unambiguous_dna)
    elif all(n in IUPAC.ambiguous_dna.letters for n in t):
        return Seq.Seq(t,IUPAC.ambiguous_dna)

    raise ValueError("invalid DNA sequence, base must be one of %s: %s"%(IUPAC.ambiguous_dna.letters,s))

def base_color(n):
    return {'A':'#00FF00',
            'T':'#FF0000',
            'G':'#000000',
            'C':'#0000FF'}[n.upper()]

def tm_gc(seq):
    gc = seq.count('G')+seq.count('C')
    at = seq.count('A')+seq.count('T')
    return '4x%s+2x%s=%s'%(gc,at,4*gc+2*at)

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
