
import os

from ..util.parser import TreekvParser

from ..util.namedlist import DefaultNamedList

from .. import db
from . import dbtss_block

LENGTH_THRESHOLD = 800

__all__ = ['TssvFile']

from seqview import Seqviews

class TssvFile(Seqviews):
    def __init__(self, seqdb):
        self.tissueset = []
        self.genes = DefaultNamedList()
        super(TssvFile,self).__init__()

    def tss_count_csv(self):
        content = ''
        content += ', '.join(['tss \\ tissue']+[n for n in self.tissueset]) + '\n'

        for name, start, end, gsl in self.genes:
            content += ','.join([name] + gsl.count_tags(start,stop)) + '\n'

        return content

    def write_csv(self, outputfile):
        with open(outputfile, 'w') as f:
            f.write(self.tss_count_csv())

    def load_tssv(self, filename):
        parser = TreekvParser()

        genes = DefaultNamedList() # i need ordered default dict....

        with open(filename,'r') as fileobj:
            tree = parser.readfp(fileobj, filename)

            e = None
        
            kv = tree.get_one('tss/tissues')
            if kv:
                self.tissueset = kv.value_list()
            else:
                self.tissueset = db.dbtss.get_genomeset_list()

            genes = DefaultNamedList()
            kv = tree.get_one('genes')
            if kv:
                for kv in kv.items():
                    name = kv.key
                    lq = kv.value_list()
                    gene = lq[0]
                    if lq[1]=='-':
                        start,stop = None,None
                    else:
                        start,stop = [int(x) for x in lq[1].split('-')]

                    genes[gene].append((name,start,stop))

            self.genes = []
            for gene,value in genes.items():

                e = SeqviewEntity.create_gene(symbol, gene_id)
                block = dbtss_block.DbtssBlock(e.template, tissueset)
                e.add_block(block)

                counter = 1
                for name,start,stop in value:
                    if not start or not stop:
                        p = e.template.transcript_start_site
                        start,end = p-200, p+200
                        self.genes.append((gene+" Assumed TSS", start, end, block.gsl))
                    else:
                        if not name:
                            name = "%s TSS No. %s" % (gene, counter)
                            counter += 1
                        self.genes.append((name, start, end, block.gsl))

                self.append(e)
