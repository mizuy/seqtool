import os
import cPickle as pickle

__all__ = ['GeneTable', 'NoSuchGene']

directory = os.path.dirname(os.path.abspath(__file__))
default_tab_file = os.path.join(directory,'../../_db/hgnc.tab')
default_cache_file = os.path.join(directory,'../../_cache/hgnc.cache')

class NoSuchGene(Exception):
    def __init__(self, m):
        self.gene_text = str(m)
    def __str__(self):
        return 'No Such Gene "%s"'%self.gene_text

class GeneTable(object):
    def __init__(self, tab_file=default_tab_file, cache_file=default_cache_file):
        self.tab_file = tab_file
        self.cache_file = cache_file
        if not self.load_cache(self.cache_file):
            self.load(tab_file)

    def load_cache(self, cache_file):
        if not os.path.exists(cache_file):
            return False
        try:
            s,g = pickle.load(open(cache_file,'rb'))
        except pickle.PickleError:
            return False
        self.symbol = s
        self.gene_id = g
        return True

    def load(self, tab_file):
        self.symbol = {}
        self.gene_id = {}
        with open(self.tab_file, 'r') as f:
            lines = f.readlines()
            title = lines[0].strip().split('\t')
            assert(title[1]=="Approved Symbol")
            assert(title[2]=="Entrez Gene ID (mapped data supplied by NCBI)")
            for line in lines[1:]:
                l = line.strip().split('\t')
                if len(l) < 4:
                    continue
                try:
                    symbol = l[1]
                    gene_id = int(l[2])
                except:
                    continue
                self.symbol[gene_id] = symbol
                self.gene_id[symbol] = gene_id

        pickle.dump((self.symbol, self.gene_id), open(self.cache_file,'wb'))

    def get_symbol(self, gene_id):
        try:
            return self.symbol[gene_id]
        except LookupError:
            raise NoSuchGene(gene_id)
    
    def get_gene_id(self, symbol):
        try:
            return self.gene_id[symbol]
        except LookupError:
            raise NoSuchGene(symbol)
