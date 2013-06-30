from . import to_unambiguous_seq


def drop_small_region(items, length):
    ret = []
    for p, q in items:
        if q-p >= length:
            ret.append((p, q))
    return ret

def region_growing(items, gap=0):
    if len(items) < 2:
        return items
    p, q = items[0]
    rest = items[1:]
    for i in range(len(rest)):
        r, s = rest[i]
        if r < q + gap:
            q = max(q, s)
            continue
        else:
            return [(p,q)]+region_growing(rest[i:],gap)
    return [(p,q)]

# http://genomewiki.ucsc.edu/index.php/CpG_Islands
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
    for i in range(0,l):
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

def cpg_obs_per_exp(seq):
    """
    CpG islands in vertebrate genomes {Gardiner-Garden, 1987, p02206}
    'Obs/Exp CpG' = N * 'Number of CpG' / 'Number of C' * 'Number of G'
    where, N is the total number of nucleotide in the sequence being analyzed.

    2 >= Obs/Exp >= 0

    >>> cpg_obs_per_exp('GC')
    0.0
    >>> cpg_obs_per_exp('CG')
    2.0
    >>> cpg_obs_per_exp('CGCGCGCG')
    2.0
    >>> cpg_obs_per_exp('CGGCCGGCCGGC')
    1.0
    """
    n = len(seq)

    seqstr = str(seq)
    n = len(seqstr)
    c = seqstr.count('C')
    g = seqstr.count('G')
    cpg = seqstr.count('CG')
    oe = 1.*n*cpg/(c*g) if (c*g)!=0 else 0

    return oe

#def cpg_sites(seq): return [i for i in range(0,len(seq)-1) if seq[i]=='C' and seq[i+1]=='G']
def cpg_sites(seq, range_=(None,None)):
    """
    return all cpg location of seq

    >>> cpg_sites('ATGCCGCGATCG')
    [4,6,10]
    >>> cpg_sites('ATGCCGCGATCG',(2,6))
    [4,6]
    """
    seqstr = str(seq)
    length = len(seq)

    p,q = range_
    p = p or 0
    if q:
        q = min(q+1, length)
    else:
        q = length

    seqstr = seqstr

    j = p
    while 1:
        j = seqstr.find('CG', j)
        if j<0 or q<=j:
            break
        yield j
        j += 1

def is_cpg(seq,i):
    """
    >>> is_cpg('ATGCGC', 3)
    True
    >>> is_cpg('ATGCGC', 4)
    False
    """
    if i+1>=len(seq):
        return False
    if seq[i:i+2]=='CG':
        return True
    return False

def count_cpg(seq, range_=(None,None)):
    """
    >>> count_cpg('ATGCCGCGATCG')
    3
    >>> count_cpg('ATGCCGCGATCG',(2,6))
    2
    >>> count_cpg('ATGCCGCGATCG',(2,6))
    30
    """
    p,q = range_
    p = p or 0
    q = q or -1

    return str(seq).count('CG',p,q)


def _bisulfite_conversion(seq):
    """
    >>> _bisulfite_conversion('ATGCCGATGC')
    Seq.Seq('ATGTYGATGT')
    """
    seqstr = str(seq)
    l = len(seqstr)
    muta = seq.tomutable()

    j = 0
    while j+1<l:
        j = seqstr.find('C', j)

        if not (0<=j and j+1<l):
            break
        if seqstr[j+1]=='G':
            muta[j] = 'Y'
        else:
            muta[j] = 'T'

        j += 1

    return muta.toseq()

def bisulfite_conversion(seq, sense=True):
    return asymmetric_conversion(seq, _bisulfite_conversion, sense=sense)

def asymmetric_conversion(seq, conv, sense):
    if sense:
        return conv(seq)
    else:
        return conv(seq.reverse_complement()).reverse_complement()

def to_unambiguous(bsseq, methyl=True):
    trans = str.maketrans('YR','CG') if methyl else str.maketrans('YR','TA')
    return to_unambiguous_seq(str(bsseq).translate(trans))

def bisulfite_conversion_unambiguous(seq, sense, methyl):
    return to_unambiguous(bisulfite_conversion(seq, sense=sense), methyl=methyl)

def bisulfite(seq, methyl, sense=True):
    """
    >>> import Bio.Seq as Seq
    >>> bisulfite(Seq.Seq('ATGCGC'), True)
    Seq('ATGCGT', Alphabet())
    >>> bisulfite(Seq.Seq('ATGCGC'), False)
    Seq('ATGTGT', Alphabet())
    """
    key = '_bisulfite_' + ('met' if methyl else 'unmet') + '_' + ('sense' if sense else 'asense')

    if not hasattr(seq, key):
        ret = bisulfite_conversion_unambiguous(seq, sense=sense, methyl=methyl)
        setattr(seq, key, ret)

    return getattr(seq, key)

def gc_ratio(seq):
    c = seq.count('C')
    c += seq.count('G')
    c += seq.count('Y')
    c += seq.count('R')
    return 100. * c/len(seq)

def _c2t_conversion(seq):
    muta = seq.tomutable()
    for i,c in enumerate(muta):
        if c=='C':
            muta[i] = 'Y'
    return muta.toseq()

def c2t_conversion(seq, sense=True):
    return asymmetric_conversion(seq, lambda x: _c2t_conversion(x), sense)


class BisulfiteTemplate:
    def __init__(self, origin_seq):
        self.origin = origin_seq
        self.sense = bisulfite_conversion(True)
        self.asense = bisulfite_conversion(False)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
