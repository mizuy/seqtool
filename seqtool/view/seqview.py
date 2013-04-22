from __future__ import absolute_import

from Bio import Restriction
from Bio import Seq
from Bio.Alphabet import IUPAC

import os

from ..util import xmlwriter
from ..nucleotide.primer import Primer,Primers
from ..util.parser import TreekvParser

from ..util.subfs import SubFileSystem
from ..util.dirutils import Filepath
from .. import db

from .css import seqview_css

from . import template as temp
from . import block, bsa_block, dbtss_block

from .baseseq_renderer import BaseseqRenderer
from .outline_renderer import OutlineRenderer

from collections import OrderedDict

LENGTH_THRESHOLD = 800

__all__ = ['Seqviews', 'TssvFile']

class SeqviewEntity(object):
    def __init__(self, name, template):
        self.name = name
        self.template = template

        self.primers = Primers()
        self.restrictions = []

        self.pcrs = block.PcrsBlock(self.template, self.primers)
        self.bs_pcrs = block.BsPcrsBlock(self.template, self.primers)
        self.rt_pcrs = block.RtPcrsBlock(self.template, self.primers)

        self.dbtss = dbtss_block.DbtssBlock(self.template)
        self.bsa = bsa_block.BsaBlock(self.bs_pcrs)

        self.blocks = [self.dbtss,
                     self.bsa,
                     self.pcrs,
                     self.bs_pcrs,
                     self.rt_pcrs]

        self.show_bisulfite = False

    def set_show_bisulfite(self, b):
        self.show_bisulfite = b

    @classmethod
    def load_seqv(cls, filename):
        parser = TreekvParser()

        inputp = Filepath(filename)
        relative_path = lambda x: inputp.relative(x)

        post_processing = []

        def ignore_child(kv):
            if kv.has_items():
                print 'Ignoring unkown child items'
                print kv.lineinfo.error_msg()

        with open(filename,'r') as fileobj:
            bn = os.path.basename(filename)

            tree = parser.readfp(fileobj, filename)

            e = None

            kv = tree.get_one('general/genbank')
            if kv:
                with open(relative_path(kv.value), 'r') as f:
                    e = cls.create_genbank(bn, f.read())
            else:
                kv = tree.get_one('general/sequence')
                if kv:
                    e = cls.create_sequence(bn, kv.value)
                else:
                    kv = tree.get_one('general/gene')
                    if kv:
                        try:
                            gene_id, gene_symbol = db.get_gene_from_text(kv.value)
                            e = cls.create_gene(gene_symbol, gene_id)
                        except db.NoSuchGene,e:
                            print 'gene entry: No such Gene %s'%kv.value
            if not e:
                print 'unkown sequence. exit.'
                return None

            kv = tree.get_one('general/primers')
            if kv:
                e.primers.load_file(relative_path(kv.value))

            kv = tree.get_one('general/tss')
            if kv:
                e.dbtss.set_tissueset(kv.value_list())

            kv = tree.get_one('general/restriction')
            if kv:
                for v in kv.value_list():
                    if v not in Restriction.AllEnzymes:
                        print 'No such Restriction Enzyme: {0}'.format(v)
                        continue
                    e.restrictions.append(v)

            kv = tree.get_one('general/show_bisulfite')
            if kv:
                if kv.value in ['True','1','true','t','T']:
                    e.set_show_bisulfite(True)
                elif kv.value in ['False','0','false','f','F']:
                    e.set_show_bisulfite(False)
                else:
                    print 'Unkown Boolean value: {0}. must be True or False'.format(kv.value)

            kv = tree.get_one('general/bsa')
            kv = tree.get_one('general/bsa/result')
            if kv:
                bsa_file = relative_path(kv.value)
                def lazy():
                    e.bsa.read(relative_path(bsa_file))
                post_processing.append(lazy)
            kv = tree.get_one('general/bsa/celllines')
            if kv:
                e.bsa.set_celllines(kv.value_list())

            for kv in tree.get_one('general').get_unused():
                print kv.lineinfo.error_msg('Ignored.')

            kv = tree.get_one('primers')
            if kv:
                for kv in kv.items():
                    e.primers.add(Primer(kv.key, kv.value))

            def pcr_line(kv):
                name = kv.key.split(',')[0].strip()
                ls = kv.value.split(',')
                if len(ls)!=2:
                    print 'you must specify 2 primer names separated by "," for each pcr: %s'%name
                    return None, None, None
                fw = ls[0].strip()
                rv = ls[1].strip()
                return name, fw, rv

            pp = tree.get_one('pcr')
            if pp:
                for kv in pp.items():
                    name, fw, rv = pcr_line(kv)
                    if not name:
                        continue
                    e.pcrs.add(name, fw, rv)

            pp = tree.get_one('rt_pcr')
            if pp:
                for kv in pp.items():
                    name, fw, rv = pcr_line(kv)
                    if not name:
                        continue
                    e.rt_pcrs.add(name, fw, rv)

            pp = tree.get_one('bs_pcr')
            if pp:
                for kv in pp.items():
                    name, fw, rv = pcr_line(kv)
                    if not name:
                        continue
                    e.bs_pcrs.add(name, fw, rv)

        for p in post_processing:
            p()
        return e

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

class Seqviews(object):
    def __init__(self):
        self.entries = []

    def append(self, entity):
        self.entries.append(entity)

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

    def load_seqv(self, filename):
        e = SeqviewEntity.load_seqv(filename)
        self.append(e)
        return e

    def load_gene(self,name, gene_id):
        e = SeqviewEntity.create_gene(name, gene_id)
        self.append(e)
        return e

class TssvFile(Seqviews):
    def __init__(self):
        self.tissueset = []
        super(TssvFile,self).__init__()

    def write_csv(self, outputfile):
        with open(outputfile, 'w') as f:
            f.write(', '.join(['tss \\ tissue']+[n for n in self.tissueset]) + '\n')
            for gt in self.entries:
                f.write(gt.dbtss.tss_count_csv())

    def load_tssv(self, filename):
        parser = TreekvParser()

        genes = OrderedDict() # i need ordered default dict....

        with open(filename,'r') as fileobj:
            tree = parser.readfp(fileobj, filename)

            e = None
        
            kv = tree.get_one('tss/tissues')
            if kv:
                for t in kv.value_list():
                    self.tissueset.append(t)

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
                        
                    if gene not in genes:
                        genes[gene] = [(name,start,stop)]
                    else:
                        genes[gene].append((name,start,stop))

            for gene in genes.keys():
                gene_id, symbol = db.get_gene_from_text(gene)

                e = SeqviewEntity.create_gene(symbol, gene_id)
                e.dbtss.set_tissueset(self.tissueset)

                for name,start,stop in genes[gene]:
                    if not start or not stop:
                        e.dbtss.add_default_tss(name)
                    else:
                        e.dbtss.add_tss(start, stop, name)

                self.append(e)
