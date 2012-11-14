from __future__ import absolute_import

import cPickle as pickle

from collections import defaultdict, OrderedDict
from .locus import Locus
import os

directory = os.path.dirname(os.path.abspath(__file__))
database_dir = os.path.join(directory,'../../_db/dbtss/')

default_cache_file = os.path.join(directory,'../../_cache/bed.cache')

import bisect

def settings(text, sp=','):
    return [[x.strip() for x in c.split(sp)] for c in text.strip().split('\n')]

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

"""
db = Dbtss()
db.load()
for tissue in db.get_tissues():
    for loc, value in tissue.search('chrX', '+', 3045, 3145):
        print loc, value
"""

class Dbtss(object):
    def __init__(self):
        self.tissues = None

        if os.path.exists(default_cache_file):
            try:
                print "loading %s"%default_cache_file
                self.tissues = pickle.load(open(default_cache_file,'rb'))
                print "done."
            except pickle.PickleError:
                print "pickle error: %s"%default_cache_file
                self.tissues = None

        if not self.tissues:
            self.tissues = OrderedDict()
            for name,fname in l_dbtss:
                tabfile = os.path.join(database_dir, fname)
                tab = TssTab(name, open(tabfile,'r'))
                self.tissues[name] = tab
            pickle.dump(self.tissues, open(default_cache_file,'wb'))

    def gets_tissue(self, name):
        return self.tissues[name]
    
    def __getitem__(self,name):
        return self.tissues[name]

class TssTab(object):
    def __init__(self, name, fileobj):
        self.name = name
        self.locs = defaultdict(list)
        self.vals = defaultdict(list)

        print 'loading tss...: %s'%name
        for l in fileobj:
            ll = l.split()
            chromosome = ll[0]
            location = int(ll[1])
            strand = ll[3].strip()=='+'
            value = int(ll[4])

            self.locs[(chromosome,strand)].append(location)
            self.vals[(chromosome,strand)].append(value)

        print 'done.'

    def search_locus(self, locus):
        return self.search(locus.chromosome, locus.strand, locus.lower, locus.higher)

    def search_locus_position(self, locus, start, stop):
        p = locus.index_5(start)
        q = locus.index_5(stop)
        if p > q:
            p,q = q,p
        return self.search(locus.chromosome, locus.strand, p, q)

    def search(self, chromosome, strand, lower, higher):
        assert(lower <= higher)
        tlocs = self.locs[(chromosome,strand)]
        tvals = self.vals[(chromosome,strand)]
        
        l_lower = bisect.bisect_left(tlocs, lower)
        l_higher   = bisect.bisect_right(tlocs, higher)

        for i in xrange(l_lower, l_higher):
            yield tlocs[i], tvals[i]

class RegionTss(object):
    def __init__(self, name, tissue, locus):
        self.name = name
        self.maxtag = 1
        self.locus = locus
        self.db = dbtss[tissue]
        for loc, val in self.db.search_locus(locus):
            self.maxtag = max(self.maxtag, val)

    def count_range(self, start, stop):
        if not start or not stop:
            iterator = self.db.search_locus(self.locus)
        else:
            iterator = self.db.search_locus_position(self.locus, start, stop)
            
        c = 0
        for loc,val in iterator:
            c += val

        return c

    def items(self):
        iterator = self.db.search_locus(self.locus)
        for k,v in iterator:
            yield self.locus.rel_5(k),v

            
dbtss = Dbtss()
