from __future__ import absolute_import

import os
from collections import OrderedDict
from cStringIO import StringIO

from .memoize import memoize
from .nucleotide import bisulfite
from .pcr import Primer, PCR, primers_write_html
from .parser import parse_file
from .primers import load_primer_list_file
from .prompt import prompt
from .dbtss import TssFile
from . import db
from . import seqtrack
from . import xmlwriter
from .listdict import ListDict

from .parser import SettingFile
from .seqview import GenomicTemplate, GenebankEntry, GenebankTssEntry, seqview_css, SubFileSystem

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

class TssvFile(object):
    def __init__(self, fileobj):
        self.parse(fileobj)

    def __iter__(self):
        return iter(self.entries)

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

        self.tissueset = tissues
        self.entries = []
        for gene in genes.keys():
            gene_id, symbol = db.get_gene_from_text(gene)

            sv = GenebankTssEntry(gene)
            sv.load_gene(gene_id)
            sv.set_tissueset(self.tissueset)

            for name,start,stop in genes[gene]:
                if not start or not stop:
                    sv.add_default_tss(name)
                else:
                    sv.add_tss(start, stop, name)
            self.entries.append(sv)


    def write_csv(self, outputfile):
        with open(outputfile, 'w') as f:
            f.write(', '.join(['tss \\ tissue']+[n for n in self.tissueset]) + '\n')
            for gt in self.entries:
                f.write(gt.tss_count_csv())

    def write_html(self, outputp):
        subfs = SubFileSystem(outputp.dir, outputp.suffix)

        with open(outputp.path,'w') as output:
            html = xmlwriter.XmlWriter(output)
            b = xmlwriter.builder(html)
            with b.html:
                with b.head:
                    with b.style(type='text/css'):
                        b.text(seqview_css)
            with b.body:
                for gt in self.entries:
                    subsubfs = subfs.get_subfs(gt.name)
                    gt.write_html(b, subsubfs)

        subfs.finish()

    

