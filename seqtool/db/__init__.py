from .locus import Locus
from .dbtss import Dbtss
from .sql import GeneTable


class NoSuchGene(Exception):
    def __init__(self, m):
        self.gene_text = str(m)
    def __str__(self):
        return 'No Such Gene "%s"'%self.gene_text


__all__ = ['initialize', 'dbtss', 'genetable', 'Locus', 'clear', 'load_dbtss', 'load_table' ,'NoSuchGene']

dbtss = None
genetable = None

def initialize(cache_dir, email):
    global dbtss
    global genetable
    dbtss = Dbtss()
    dbtss.load_cache(cache_dir)
    genetable = GeneTable(cache_dir, email)


def clear(cache_dir):
    genetable = GeneTable(cache_dir, None)
    genetable.clear()
    """
    dbtss = Dbtss()
    dbtss.load_file(args.bed_dir, args.cache_dir)
    dbtss.clear()
    """

def load_dbtss(cache_dir, bed_dir):
    dbtss = Dbtss()
    dbtss.load_file(bed_dir, cache_dir)

def load_table(cache_dir, chrom_tab_file, hgnc_tab_file, ucsc_tab_file):
    print("loading to...: ",cache_dir)
    genetable = GeneTable(cache_dir, None)
    genetable.load(chrom_tab_file, hgnc_tab_file, ucsc_tab_file)


