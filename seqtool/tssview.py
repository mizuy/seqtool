from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from collections import defaultdict, OrderedDict
from cStringIO import StringIO

import os

from .memoize import memoize
from .nucleotide import bisulfite
from .pcr import Primer, PCR, primers_write_html
from .parser import parse_file
from .primers import load_primer_list_file
from .prompt import prompt
from .dbtss import TssFile
from .db import gene
from . import seqtrack
from . import xmlwriter
from .listdict import ListDict

from .parser import SettingFile
from .seqview import GenomicTemplate, GeneBankEntry, seqview_css, SubFileSystem
from .db.dbtss import Dbtss, RegionTss
from .db import gene as genedb

class TissueSet(object):
    def __init__(self, names):
        self.tissues = names

    def __iter__(self):
        return iter(self.tissues)

class GeneTssEntry(GeneBankEntry):
    def __init__(self, gene_symbol, tissueset):
        super(GeneTssEntry, self).__init__(gene_symbol)

        self.gene_symbol = gene_symbol
        gene_id, symbol = genedb.get_gene_from_text(gene_symbol)
        self.gene_id = gene_id
        self.locus = genedb.get_gene_locus(gene_id)

        self.tissues = tissueset
        self.rtss = []
        for t in tissueset:
            self.rtss.append(RegionTss(t, t, self.locus))

        self.tsss = []

        self.load_genbank(genedb.get_genomic_context_genbank(gene_id))

    
    def add_tss_site(self, tss_name, start=None, end=None):
        self.tsss.append((tss_name, start, end))

    def tss_count_csv(self):
        content = ''
        for name, start, stop in self.tsss:
            l = [name] + [str(r.count_range(start, stop)) for r in self.rtss]
            content += ','.join(l) + '\n'
        return content

    def render_genome(self, width):
        assert(self.template)
        length = len(self.template.seq)

        t = seqtrack.TrackGroup()
        t.add(seqtrack.Track(0,10))
        t.add(seqtrack.SequenceTrack(self.template.seq, self.template.features, -1* self.template.transcript_start_site))
        for rt in self.rtss:
            t.add(seqtrack.HbarTrack('', length))
            t.add(seqtrack.RegiontssTrack(rt, self.template.seq))
        t.add(seqtrack.HbarTrack('', length))

        return t.svg(width)


    def write_html(self, b, subfs):
        """
        subfs must have 2 methods
        def write(self, filename, content_text)
        def get_link_path(self, filename)
        """
        genome_n = 'genome.svg'

        # writing svgs
        subfs.write(genome_n, self.render_genome(5000))

        # link path for svg files
        genome_l = subfs.get_link_path(genome_n)

        # writing html
        b.h1(self.gene_symbol)
        with b.div(**{'class':'images'}):
            with b.a(href=genome_l):
                b.img(src=genome_l, width='1000px')

class TssvFile(object):
    def __init__(self, filename):
        with open(filename,'r') as f:
            self.parse(f)

    def __iter__(self):
        return iter(self.entries)

    def write_html(self, b, subfs):
        for gt in self.entries:
            subsubfs = subfs.get_subfs(gt.gene_symbol)
            gt.write_html(b, subsubfs)

    def csv(self):
        content = ''
        content += ', '.join(['tss \\ tissue']+[n for n in self.tissueset]) + '\n'
        for gt in self.entries:
            content += gt.tss_count_csv()
        return content

    def parse(self, fileobj):
        tissues = []
        genes = OrderedDict() # i need ordered default dict....

        s = SettingFile()
        s.parse(fileobj)

        def error(msg,l,lineno):
            print ':%s: %s: "%s"'%(lineno,msg,l)
        
        for block in s:
            if block.name=='tss':
                for line, lineno in block:
                    ls = [x.strip() for x in line.split(':')]
                    if len(ls)!=2:
                        error('unknown line', line, lineno)
                        continue
                    tissues.append(ls[0])
                    # TODO ls[1] for tss tab file.
            elif block.name=='genes':
                for line, lineno in block:
                    lp = [x.strip() for x in line.split(':')]
                    name = lp[0]
                    lq = [x.strip() for x in lp[1].split(',')]
                    gene = lq[0]
                    if lq[1]=='-':
                        start,stop = None,None
                    else:
                        start,stop = [int(x) for x in lq[1].split('-')]
                        
                    if gene not in genes:
                        genes[gene] = [(name,start,stop)]
                    else:
                        genes[gene].append((name,start,stop))

        self.tissueset = TissueSet(tissues)
        self.entries = []
        for gene in genes.keys():
            g = GeneTssEntry(gene, self.tissueset)
            for name,start,stop in genes[gene]:
                g.add_tss_site(name, start, stop)
            self.entries.append(g)

seqview_css = '''
    .images{}
    .image{border: solid 1px;}
    .template{margin-left: 1em;}
    .pcr{margin: 1em; padding: 1em; border: solid 1px;}
    .products{margin-left: 2em;}
    .copybox{margin-left:4em;}

    .primerpairtable{
      font-family: monospace
    }
'''

def main():
    import sys, os
    from argparse import ArgumentParser

    parser = ArgumentParser(prog='tssview', description='pretty HTML+SVG report of dbtss of multiple genes')
    parser.add_argument("tssvfile", help=".tssv file")
    parser.add_argument("-o", "--output", dest="output", help="output filename")
    
    args = parser.parse_args()

    inputfile = args.tssvfile
    
    if args.output:
        outputfile = args.output
    else:
        output_dir = os.path.abspath(os.path.dirname(inputfile))
        outputfile = os.path.join(output_dir, os.path.basename(inputfile).rpartition('.')[0] + '.html')


    output_dir = os.path.abspath(os.path.dirname(outputfile))
    output_basename = os.path.basename(outputfile).rpartition('.')[0]

    tssv = TssvFile(inputfile)

    output_csv = os.path.join(output_dir, output_basename + '.csv')
    open(output_csv, 'w').write(tssv.csv())

    subfs = SubFileSystem(output_dir, output_basename)

    with open(outputfile,'w') as output:
        html = xmlwriter.XmlWriter(output)
        b = xmlwriter.builder(html)
        with b.html:
            with b.head:
                with b.style(type='text/css'):
                    b.text(seqview_css)
        with b.body:
            tssv.write_html(b, subfs)

    subfs.finish()

if __name__=='__main__':
    main()
