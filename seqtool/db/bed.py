import pickle as pickle

from collections import defaultdict, OrderedDict
from ..util.namedlist import NamedList
import numpy
import os, glob

__all__ = ['BedDB']

"""
db = BedDB()
db.load()
for tissue in db.get_tissues():
    for loc, value in tissue.search('chrX', '+', 3045, 3145):
        print loc, value
"""

"""
GenomeLocus = Genome * Locus
Genome = Doublestrand * chrom
Doublestrand = 2* SortedLocVal
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

        for i in range(l_lower, l_higher):
            yield self.locs[i], self.vals[i]

class DoubleStrand(object):
    def __init__(self, sense, antisense):
        assert(isinstance(sense, SortedLocVal))
        assert(isinstance(antisense, SortedLocVal))
        self.sense = sense
        self.antisense = antisense

    def get_strand(self, sense):
        return self.sense if sense else self.antisense

    def search(self, sense, lower, higher):
        assert(lower <= higher)
        return self.get_strand(sense).search(lower,higher)

    def search_pos(self, pos):
        return self.search(pos.sense, pos.lower, pos.higher)

    def search_pos_index(self, pos, start, stop):
        p = pos.index_5(start)
        q = pos.index_5(stop)
        if p > q:
            p,q = q,p
        return self.search(pos.sense, p, q)

class Genome(object):
    def __init__(self, name):
        self.name = name
        self.chroms = {}

    def read(self, filename):
        with open(filename,'r') as fileobj:
            self.readfp(fileobj)

    def readfp(self, fileobj):
        self.chroms = {}

        locs = defaultdict(lambda : defaultdict(list))
        vals = defaultdict(lambda : defaultdict(list))

        for l in fileobj:
            """
            bed/
            # Definition of each column (tab limited):
            -chromosome
            -TSStag_start
            -tag_end
            -strand
            -Number of tags
            """
            ll = l.split()
            chrom = ll[0]
            location = int(ll[1])
            strand = 0 if ll[3].strip()=='+' else 1
            value = int(ll[4])

            locs[chrom][strand].append(location)
            vals[chrom][strand].append(value)

        for ch in list(locs.keys()):
            sense =     SortedLocVal(numpy.array(locs[ch][0]),
                                     numpy.array(vals[ch][0]))
            antisense = SortedLocVal(numpy.array(locs[ch][1]),
                                     numpy.array(vals[ch][1]))
            self.chroms[ch] = DoubleStrand(sense, antisense)

    def get_ds(self, chrom):
        return self.chroms[chrom]

    def search(self, chrom, strand, lower, higher):
        assert(lower <= higher)
        return self.chroms[chrom].search(strand,lower,higher)

    def get_locus(self, locus):
        return GenomeLocus(locus.pos, self.get_ds(locus.chrom), self.name)

class GenomeLocus(object):
    def __init__(self, pos, ds, genome_name):
        self.pos = pos
        self.ds = ds
        self.name = genome_name

        self.maxtag = 1
        for loc, val in self.ds.search_pos(self.pos):
            self.maxtag = max(self.maxtag, val)

    def count_range(self, start, stop):
        if not start or not stop:
            iterator = self.ds.search_pos(self.pos)
        else:
            iterator = self.ds.search_pos_index(self.pos, start, stop)
        c = 0
        for loc,val in iterator:
            c += val
        return c

    def items(self):
        iterator = self.ds.search_pos(self.pos)
        for k,v in iterator:
            yield self.pos.rel_5(k), v

class GenomeSetLocus(object):
    def __init__(self, locus):
        self.locus = locus
        self.gl = NamedList()
        self.maxtag = 1

    def add_genome(self, name, genome):
        r = genome.get_locus(self.locus)
        self.gl[name] = r
        self.maxtag = max(self.maxtag, r.maxtag)

    def __iter__(self):
        return iter(self.gl)

    def get_genomelocus(self, name):
        return self.gl[name]

    def count_tags_header(self):
        return self.gl.keys()

    def count_tags(self, start, stop):
        return [str(r.count_range(start, stop)) for r in self.gl]

class BedDB(object):
    def __init__(self):
        self.genomeset = None

    def cache_file(self, cache_dir):
        return os.path.join(cache_dir, 'beddb.pickle.cache')

    def load_cache(self, cache_dir):
        cache_file = self.cache_file(cache_dir)

        # load from cache
        if os.path.exists(cache_file):
            try:
                self.genomeset = pickle.load(open(cache_file,'rb'))
            except pickle.PickleError:
                print("pickle error: %s"%cache_file)
                self.genomeset = None

    def load_file(self, bed_dir, cache_dir):
        cache_file = self.cache_file(cache_dir)

        # De-novo loading
        self.genomeset = OrderedDict()

        for fname in glob.glob(os.path.join(bed_dir,'*.bed')):
            print('loading: ',fname)
            name = os.path.splitext(os.path.basename(fname))[0]
            g = Genome(name)
            g.read(fname)
            self.genomeset[name] = g

        # save to cache
        pickle.dump(self.genomeset, open(cache_file,'wb'), -1)

    def load(self, database_dir, cache_dir):
        self.load_cache(cache_dir)

        if not self.genomeset:
            self.load_file(database_dir, cache_dir)

    def get_genomeset_locus(self, tissues, locus):
        gsl = GenomeSetLocus(locus)
        for genome_name in tissues:
            if genome_name in self.genomeset:
                gsl.add_genome(genome_name, self.genomeset[genome_name])
            else:
                print('unkown tissue: ', genome_name)
        return gsl
