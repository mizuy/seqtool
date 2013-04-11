from __future__ import absolute_import

from Bio import Restriction
from Bio import Seq
from Bio.Alphabet import IUPAC

import os

from ..util import xmlwriter
from ..nucleotide.pcr import Primer
from ..parser import SettingFile

from ..util.subfs import SubFileSystem
from ..util.dirutils import Filepath
from .. import db

from .css import seqview_css

from . import template as temp
from . import block

from .baseseq_renderer import BaseseqRenderer
from .outline_renderer import OutlineRenderer

from collections import OrderedDict

LENGTH_THRESHOLD = 800

__all__ = ['Seqview', 'SeqvFile']

class SeqviewEntity(object):
    def __init__(self, name, template):
        self.name = name
        self.template = template

        self.primers = block.PrimersHolder()
        self.restrictions = []

        self.pcrs = block.PcrsBlock(self.template, self.primers)
        self.bs_pcrs = block.BsPcrsBlock(self.template, self.primers)
        self.rt_pcrs = block.RtPcrsBlock(self.template, self.primers)

        self.dbtss = block.DbtssBlock(self.template)
        self.bsa = block.BsaBlock(self.bs_pcrs)

        self.blocks = [self.dbtss,
                     self.bsa,
                     self.pcrs,
                     self.bs_pcrs,
                     self.rt_pcrs]

        self.show_bisulfite = False

    def set_show_bisulfite(self, b):
        self.show_bisulfite = b

    @classmethod
    def create_genbank(cls, name, content):
        template = temp.GenbankTemplate(content, None)
        return cls(name, template)

    @classmethod
    def create_sequence(cls, name, content):
        sequence = Seq.Seq(content.upper(), IUPAC.unambiguous_dna)
        template = temp.SequenceTemplate(sequence)
        return cls(name, template)

    @classmethod
    def create_gene(cls, name, gene_id):
        locus = db.get_gene_locus(gene_id).expand(1000,1000)
        template = temp.GenbankTemplate(db.get_locus_genbank(locus), locus)
        return cls(name, template)

    def svg_genome(self):
        # todo change track_genome to svg_genome(self, t)

        scale = 1.
        length = len(self.template.seq)
        if length > LENGTH_THRESHOLD:
            scale = 1.*length/LENGTH_THRESHOLD

        t = OutlineRenderer(scale)

        t.add_padding(10)
        start = -1* self.template.transcript_start_site
        t.add_sequence_track(self.template.seq, self.template.features, start)

        for block in self.blocks:
            t.add_hline(length, 10)
            block.svg_genome(t)

        return t.svg()

    def svg_transcript(self):
        t = OutlineRenderer(1)

        for tr in self.template.transcripts:
            t.add_transcript_track(tr.name, tr.seq, tr.feature)
            self.rt_pcrs.svg_transcript(t, tr)
        return t.svg()

    def svg_baseseq(self):
        aseq = BaseseqRenderer(self.template.seq, self.show_bisulfite)

        for p in self.primers:
            aseq.add_primer(p)

        aseq.add_restriction_batch(Restriction.RestrictionBatch(self.restrictions))

        return aseq.track().svg()

    def has_transcripts(self):
        return not not self.template.transcripts

    def write_html(self, b, subfs):
        svgs = [('genome.svg', 'Outline', self.svg_genome(), True),
                ('baseseq.svg', 'Base Sequence', self.svg_baseseq(), False)]
        if self.has_transcripts():
            svgs.insert(1, ('transcript.svg', 'Transcript', self.svg_transcript(), True))

        b.h2(self.name)
        for filename, name, svg, show in svgs:
            subfs.write(filename, svg)
            link = subfs.get_link_path(filename)

            b.h3(name)
            with b.a(href=link):
                if show:
                    #b.write_raw(svg.svg_node())
                    b.img(src=link, width='1000px')
                else:
                    with b.a(href=link):
                        b.text(filename)

        #with b.div(**{'class':'primers'}):
        #    b.h2('Primers')
        #    primers_write_html(b.get_writer(), self.primers)

        b.h3('Analysis')
        for block in self.blocks:
            block.write_html(b, subfs)

class Seqview(object):
    def __init__(self):
        self.entries = []

    def append(self, entity):
        self.entries.append(entity)

    def top(self):
        if len(self.entries)>0:
            return self.entries[-1]
        else:
            return None

    def write_html(self, outputp):
        subfs = SubFileSystem(outputp.dir, outputp.prefix)

        with open(outputp.path,'w') as output:
            html = xmlwriter.XmlWriter(output)
            b = xmlwriter.builder(html)
            with b.html:
                with b.head:
                    with b.style(type='text/css'):
                        b.text(seqview_css)
            with b.body:
                count = 0
                for gt in self.entries:
                    count += 1
                    name = gt.name or '%s'%count
                    subsubfs = subfs.get_subfs(name)
                    gt.write_html(b, subsubfs)

        subfs.finish()

    def write_csv(self, outputfile):
        with open(outputfile, 'w') as f:
            f.write(', '.join(['tss \\ tissue']+[n for n in self.tissueset]) + '\n')
            for gt in self.entries:
                f.write(gt.dbtss.tss_count_csv())

    def load_genbank(self, name, content):
        self.append(SeqviewEntity.create_genbank(name, content))

    def load_sequence(self, name, content):
        self.append(SeqviewEntity.create_sequence(name, content))

    def load_gene(self, name, gene_id):
        self.append(SeqviewEntity.create_gene(name, gene_id))


class SeqvFile(Seqview):
    def __init__(self):
        super(SeqvFile,self).__init__()

    def load_genbankentry(self, genbankentry):
        self.append(genbankentry)

    def load_seqvfileentry(self, filename):
        inputp = Filepath(filename)
        relative_path = lambda x: inputp.relative(x)

        with open(filename,'r') as fileobj:
            bn = os.path.basename(filename)
            for category, name, value, em in self.parse(fileobj):
                e = self.top()
                if not name and not value:
                    pass
                else:
                    if category == 'general':
                        if name=='genbank':
                            with open(relative_path(value), 'r') as f:
                                self.load_genbank(bn, f.read())
                        elif name=='sequence':
                            self.load_sequence(bn, value)
                        elif name=='gene':
                            try:
                                gene_id, gene_symbol = db.get_gene_from_text(value)
                            except db.NoSuchGene,e:
                                em('gene entry: No such Gene %s'%value)
                                continue
                            self.load_gene(bn, gene_id)
                        elif name=='primers':
                            e.primers.load(relative_path(value))
                        elif name=='tss':
                            e.set_tissueset([x.strip() for x in value.split(',')])
                        elif name=='bsa':
                            e.bsa.set_celllines([x.strip() for x in value.split(',')])

                        elif name =='restriction':
                            for v in value.split(','):
                                if v not in Restriction.AllEnzymes:
                                    em('No such Restriction Enzyme: {0}'.format(v))
                                    continue
                                e.restrictions.append(v)

                        elif name == 'show_bisulfite':
                            if value in ['True','1','true','t','T']:
                                e.set_show_bisulfite(True)
                            elif value in ['False','0','false','f','F']:
                                e.set_show_bisulfite(False)
                            else:
                                em('Unkown Boolean value: {0}. must be True or False'.format(value))

                    elif category == 'primer':
                        e.primers.add(Primer(name, value))
                    elif category == 'motif':
                        pass
                    elif category in ['pcr','rt_pcr','bs_pcr']:
                        name = name.split(',')[0].strip()
                        ls = value.split(',')
                        if len(ls)!=2:
                            em('you must specify 2 primer names separated by "," for each pcr: %s'%name)
                            continue
                        fw = ls[0].strip()
                        rv = ls[1].strip()

                        if category == 'pcr':
                            e.pcrs.add(name, fw, rv)
                        elif category == 'rt_pcr':
                            e.rt_pcrs.add(name, fw, rv)
                        elif category == 'bs_pcr':
                            e.bs_pcrs.add(name, fw, rv)
                        else:
                            raise ValueError
                    elif category == 'bsa':
                        n = [n.strip() for n in name.split(',')]
                        if not len(n)>=2:
                            em('each bsa must have at least 2 key; cell line name and pcr name')
                            continue
                        pcrname = n[0].strip()
                        cellline = n[1].strip().upper()
                        #annotations = n[2:]
                        if not pcrname or not cellline:
                            em('empty pcr or cellline name: %s, %s'%(pcrname,cellline))
                            continue

                        e.bsa.add(cellline, pcrname, value.strip().upper())
                    else:
                        em('unkown category: %s'%category)

    def parse(self, fileobj):
        s = SettingFile()
        s.parse(fileobj)

        def error(msg,l,lineno):
            print ':%s: %s: "%s"'%(lineno,msg,l)
        
        for block in s:
            category = block.name
            yield category, None, None, lambda x:error(x, block.line, block.lineno)
            for line,lineno in block:
                ls = line.split(':')
                if len(ls)!=2:
                    error('unknown line', line, lineno)
                    continue
                name = ls[0].strip()
                value = ls[1].strip()
                yield category, name, value, lambda x:error(x, line, lineno)

class TssvFile(Seqview):
    def __init__(self, fileobj):
        super(TssvFile,self).__init__()
        self.parse(fileobj)

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
        for gene in genes.keys():
            gene_id, symbol = db.get_gene_from_text(gene)

            self.load_gene(symbol, gene_id)
            self.top().dbtss.set_tissueset(self.tissueset)

            for name,start,stop in genes[gene]:
                if not start or not stop:
                    self.top().dbtss.add_default_tss(name)
                else:
                    self.top().dbtss.add_tss(start, stop, name)
