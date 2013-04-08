from __future__ import absolute_import

from Bio.SeqUtils import GC

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

def is_cpg(seq,i):
    if i+1>=len(seq):
        return False
    if seq[i]=='C' and seq[i+1]=='G':
        return True
    return False

def count_cpg(seq):
    count = 0
    for i,n in enumerate(seq[:-1]):
        if n=='C' and seq[i+1]=='G':
            count += 1
    return count

def bisulfite_conversion(seq, methyl):
    muta = seq.tomutable()
    for i,c in enumerate(muta[:-1]):
        if c=='C':
            if not(muta[i+1]=='G' and methyl):
                muta[i] = 'T'
    if muta[-1]=='C':
        muta[-1]='T'
    return muta.toseq()

def bisulfite(seq, methyl, sense=True):
    key = '_bisulfite_' + ('met' if methyl else 'unmet') + '_' + ('sense' if sense else 'asense')

    if not hasattr(seq, key):
        if sense:
            ret = bisulfite_conversion(seq, methyl)
        else:
            ret = bisulfite_conversion(seq.reverse_complement(), methyl)
        setattr(seq, key, ret)

    return getattr(seq, key)

def gc_ratio(seq):
    return GC(seq)
