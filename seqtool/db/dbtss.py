from __future__ import absolute_import

import cPickle as pickle

from collections import defaultdict, OrderedDict
from .locus import Locus
from ..util.namedlist import NamedList
import numpy
import os

directory = os.path.dirname(os.path.abspath(__file__))

database_dir = os.path.join(directory,'../../_db/dbtss/')
default_cache_file = os.path.join(directory,'../../_cache/bed.cache')

import bisect

def settings(text, sp=','):
    return [[x.strip() for x in c.split(sp)] for c in text.strip().split('\n')]

"""
bed/
# Definition of each column (tab limited):
-chromosome
-TSStag_start
-tag_end
-strand
-Number of tags
"""


t_dbtss = """
Brain1: adult-tissue/ambion_brain.bed
Brain2: adult-tissue/ambion_brain2.bed
C-Brain: adult-tissue/clontech_brain.bed
Heart: adult-tissue/ambion_heart.bed
C-Heart: adult-tissue/clontech_heart.bed
Lung: adult-tissue/ambion_lung.bed
Breast: adult-tissue/ambion_breast.bed
Kidney: adult-tissue/ambion_kidney.bed
C-Kidney: adult-tissue/clontech_kidney.bed
Liver: adult-tissue/ambion_liver.bed
Colon: adult-tissue/ambion_colon.bed
Lymph: adult-tissue/ambion_lymph.bed
Adipose: adult-tissue/ambion_adipose.bed
Muscle: adult-tissue/ambion_muscle.bed
Thyroid: adult-tissue/ambion_thyroid.bed
Adrenal: adult-tissue/ambion_adrenal.bed
Ovary: adult-tissue/ambion_ovary.bed
Prostate: adult-tissue/ambion_prostate.bed
Testis: adult-tissue/ambion_testis.bed
"""
l_dbtss = settings(t_dbtss, ':')
tissues = [k for k,v in l_dbtss]

"""
db = Dbtss()
db.load()
for tissue in db.get_tissues():
    for loc, value in tissue.search('chrX', '+', 3045, 3145):
        print loc, value
"""

class SortedLocVal(object):
    def __init__(self, locations, values):
        self.locs = locations
        self.vals = values

    def search(self, lower, higher):
        assert(lower <= higher)

        # faster than bisect
        l_lower = self.locs.searchsorted(lower, side='left')
        l_higher = self.locs.searchsorted(higher, side='right')

        #l_lower = self.locs.bisect.bisect_left(tlocs, lower)
        #l_higher   = bisect.bisect_right(tlocs, higher)

        for i in xrange(l_lower, l_higher):
            yield self.locs[i], self.vals[i]

class DoubleStrand(object):
    def __init__(self, sense, antisense):
        self.sense = sense
        self.antisense = antisense

    def get_strand(self, sense):
        return self.sense if sense else self.antisense

    def search(self, sense, lower, higher):
        assert(lower <= higher)
        strand = self.get_strand(sense)
        
        return strand.search(lower,higher)

    def search_pos(self, pos):
        return self.search(pos.sense, pos.lower, pos.higher)

    def search_pos_index(self, pos, start, stop):
        p = pos.index_5(start)
        q = pos.index_5(stop)
        if p > q:
            p,q = q,p
        return self.search(pos.sense, p, q)


class Dbtss(object):
    def __init__(self):
        self.tissues = None

        if os.path.exists(default_cache_file):
            try:
                #print "loading cache %s"%default_cache_file
                self.tissues = pickle.load(open(default_cache_file,'rb'))
                #print "done."
            except pickle.PickleError:
                print "pickle error: %s"%default_cache_file
                self.tissues = None

        if not self.tissues:
            self.load()
            pickle.dump(self.tissues, open(default_cache_file,'wb'), -1)

    def load(self):
        self.tissues = OrderedDict()
        for name,fname in l_dbtss:
            self.tissues[name] = {}

            bed_file = os.path.join(database_dir, fname)
            with open(bed_file,'r') as fileobj:
                locs = defaultdict(lambda : defaultdict(list))
                vals = defaultdict(lambda : defaultdict(list))

                print 'loading bed...: %s'%bed_file
                for l in fileobj:
                    ll = l.split()
                    chromosome = ll[0]
                    location = int(ll[1])
                    strand = 0 if ll[3].strip()=='+' else 1
                    value = int(ll[4])

                    locs[chromosome][strand].append(location)
                    vals[chromosome][strand].append(value)

                for ch in locs.keys():
                    sense =     SortedLocVal(numpy.array(locs[ch][0]),
                                             numpy.array(vals[ch][0]))
                    antisense = SortedLocVal(numpy.array(locs[ch][1]),
                                             numpy.array(vals[ch][1]))
                    self.tissues[name][ch] = DoubleStrand(sense, antisense)

    def get_strand(self, tissue, chrom):
        return self.tissues[tissue][chrom]

    def search(self, tissue, chrom, strand, lower, higher):
        assert(lower <= higher)
        return self.tissues[tissue][chrom].search(strand,lower,higher)

_dbtss = None
def get_dbtss():
    global _dbtss
    if not _dbtss:
        _dbtss = Dbtss()
    return _dbtss


class TissueLocus(object):
    def __init__(self, tissue, locus):
        self.name = tissue
        self.locus = locus
        self.pos = locus.pos
        self.db = get_dbtss().get_strand(tissue, locus.chrom)

        self.maxtag = 1
        for loc, val in self.db.search_pos(locus.pos):
            self.maxtag = max(self.maxtag, val)

    def count_range(self, start, stop):
        if not start or not stop:
            iterator = self.db.search_pos(self.locus.pos)
        else:
            iterator = self.db.search_pos_index(self.locus.pos, start, stop)
            
        c = 0
        for loc,val in iterator:
            c += val

        return c

    def items(self):
        iterator = self.db.search_pos(self.locus.pos)
        for k,v in iterator:
            yield self.pos.rel_5(k),v

class TissuesetLocus(object):
    def __init__(self, tissue_names, locus):
        self.tissue_names = tissue_names
        self.tl = NamedList()
        self.maxtag = 1

        for tissue in tissue_names:
            self.tl.append(TissueLocus(tissue, locus))

        for r in self.tl:
            self.maxtag = max(self.maxtag, r.maxtag)

    def __iter__(self):
        return iter(self.tl)

    def get_tissue(self, name):
        return self.tl[name]

    def count_tags_header(self):
        return self.tissue_names
        
    def count_tags(self, start, stop):
        return [str(r.count_range(start, stop)) for r in self.tl]
