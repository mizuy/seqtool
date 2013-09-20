__all__ = ['seq_cpg_analysis']

def drop_small_region(items, length):
    """
    >>> drop_small_region([(10,20),(100,300),(1000,2000)], 0)
    [(10, 20), (100, 300), (1000, 2000)]
    >>> drop_small_region([(10,20),(100,300),(1000,2000)], 200)
    [(100, 300), (1000, 2000)]
    >>> drop_small_region([(10,20),(100,300),(1000,2000)], 201)
    [(1000, 2000)]
    >>> drop_small_region([(10,20),(100,300),(1000,2000)], 2000)
    []
    """
    ret = []
    for p, q in items:
        if q-p >= length:
            ret.append((p, q))
    return ret

def region_growing(items, gap=0):
    """
    >>> region_growing([(10,20),(100,300),(1000,2000)])
    [(10, 20), (100, 300), (1000, 2000)]
    >>> region_growing([(10,20),(15,30),(30,50)])
    [(10, 50)]
    """
    if len(items) < 2:
        return items
    p, q = items[0]
    rest = items[1:]
    for i in range(len(rest)):
        r, s = rest[i]
        assert(p <= r) # invalid input
        #         p             q
        #                r            s
        if r <= q + gap:
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
        self.met_criteria = False
        self.islands = []
        self.start = 0

    def p(self, i, met_criteria):
        """point each base"""
        if self.met_criteria and (not met_criteria):
            """region exit"""
            self.islands.append((self.start, i-1))
            self.met_criteria = False
        elif (not self.met_criteria) and met_criteria:
            """region enter"""
            self.start = i
            self.met_criteria = True

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

        island_criteria = (gcp > 0.55) and (oe > 0.65)
        sr.p(i, island_criteria)

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

