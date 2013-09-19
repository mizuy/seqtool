

from Bio import Restriction
from Bio import Seq
from Bio.Alphabet import IUPAC

import os

from ..util import report
from ..nucleotide.primer import Primer,Primers
from ..util.parser import TreekvParser

from ..util.dirutils import Filepath
from .. import db

from . import template as temp
from . import block, bsa_block, dbtss_block

from .baseseq_renderer import BaseseqRenderer
from .outline_renderer import OutlineRenderer

LENGTH_THRESHOLD = 20000

__all__ = ['Seqview']

class Seqview(object):
    def __init__(self, name, template):
        self.name = name
        self.template = template

        self.primers = Primers()
        self.restrictions = []

        self.pcrs = block.PcrsBlock(self.template, self.primers)
        self.bs_pcrs = block.BsPcrsBlock(self.template, self.primers)
        self.rt_pcrs = block.RtPcrsBlock(self.template, self.primers)

        self.bsa = bsa_block.BsaBlock(self.bs_pcrs)

        self.blocks = [self.bsa,
                       self.pcrs,
                       self.bs_pcrs,
                       self.rt_pcrs]

        self.show_bisulfite = False

    def add_block(self, block):
        self.blocks.append(block)

    def add_dbtss(self, tissues):
        block = dbtss_block.DbtssBlock(self.template, tissues)
        self.add_block(block)

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
                print('Ignoring unkown child items')
                print(kv.lineinfo.error_msg())

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
                            gene_id, gene_symbol = db.genetable.get_gene_from_text(kv.value)
                            e = cls.create_gene(gene_symbol, gene_id)
                        except db.NoSuchGene as e:
                            print('gene entry: No such Gene %s'%kv.value)
            if not e:
                print('unkown sequence. exit.')
                return None

            kv = tree.get_one('general/primers')
            if kv:
                for v in kv.value_list():
                    e.primers.load_file(relative_path(v))

            kv = tree.get_one('general/tss')
            if kv:
                e.add_dbtss(kv.value_list())

            kv = tree.get_one('general/restriction')
            if kv:
                for v in kv.value_list():
                    # AllEnzymes is not RestrictionBatch, but set....
                    if v not in Restriction.CommOnly.elements():
                        print('No such Restriction Enzyme: {0}'.format(v))
                        continue
                    e.restrictions.append(v)

            kv = tree.get_one('general/show_bisulfite')
            if kv:
                if kv.value in ['True','1','true','t','T']:
                    e.set_show_bisulfite(True)
                elif kv.value in ['False','0','false','f','F']:
                    e.set_show_bisulfite(False)
                else:
                    print('Unkown Boolean value: {0}. must be True or False'.format(kv.value))

            kv = tree.get_one('general/bsa')
            if kv:
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
                print(kv.lineinfo.error_msg('Ignored.'))

            kv = tree.get_one('primers')
            if kv:
                for kv in list(kv.items()):
                    e.primers.append(Primer(kv.key, kv.value))

            def pcr_line(kv):
                name = kv.key.split(',')[0].strip()
                ls = kv.value.split(',')
                if len(ls)!=2:
                    print('you must specify 2 primer names separated by "," for each pcr: %s'%name)
                    return None, None, None
                fw = ls[0].strip()
                rv = ls[1].strip()
                return name, fw, rv

            pp = tree.get_one('pcr')
            if pp:
                for kv in list(pp.items()):
                    name, fw, rv = pcr_line(kv)
                    if not name:
                        continue
                    e.pcrs.add(name, fw, rv)

            pp = tree.get_one('rt_pcr')
            if pp:
                for kv in list(pp.items()):
                    name, fw, rv = pcr_line(kv)
                    if not name:
                        continue
                    e.rt_pcrs.add(name, fw, rv)

            pp = tree.get_one('bs_pcr')
            if pp:
                for kv in list(pp.items()):
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
        locus = db.genetable.get_gene_locus(gene_id).expand(1000,1000)
        genbank = db.genetable.get_locus_genbank(locus)
        template = temp.GenbankTemplate(genbank, locus)
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

    def baseseq_renderer(self):
        aseq = BaseseqRenderer(self.template.seq, self.show_bisulfite)

        for p in self.primers:
            aseq.add_primer(p)

        aseq.add_restriction_batch(Restriction.RestrictionBatch(self.restrictions))
        return aseq

    def fasta(self):
        ret = '>{} : {}\n'.format(self.template.name, self.template.description)
        ss = str(self.template.seq)
        l = len(ss)

        for i in range(0, l, 70):
            ret += ss[i:i+100]
            ret += '\n'

        return ret

    def has_transcripts(self):
        return not not self.template.transcripts

    def html_content(self, b, toc, subfs):
        bsr =  self.baseseq_renderer()
        svgs = [('genome.svg', 'Outline', self.svg_genome(), True),
                ('baseseq.svg', 'Base Sequence', bsr.track(width = 120).svg(), False),
                ('seq.seq', 'sequence file', str(self.template.seq), False)]
        if self.has_transcripts():
            svgs.insert(1, ('transcript.svg', 'Transcript', self.svg_transcript(), True))

        def product_callback(product):
            aseq = BaseseqRenderer(product.template, self.show_bisulfite)

            for p in self.primers:
                aseq.add_primer(p)

            aseq.add_restriction_batch(Restriction.RestrictionBatch(self.restrictions))

            p = max(0, product.start_0 - 20)
            q = min(product.end_0 + 20, len(product.template))
            s = aseq.track_partial(p, q, width = None).svg()

            filename = 'pcr{}.svg'.format(product_callback.count)
            product_callback.count += 1

            subfs.write(filename, s)

            b.a('Primer alignment view', href = subfs.get_link_path(filename))
        product_callback.count = 0
    
        with toc.section(self.name):
            for filename, name, svg, show in svgs:
                subfs.write(filename, svg)
                link = subfs.get_link_path(filename)

                with toc.section(name):
                    with b.div:
                        b.h2(name)
                        with b.a(href=link):
                            if show:
                                #b.write_raw(svg.svg_node())
                                b.img(src=link, width='1000px')
                            else:
                                b.a(filename, href=link)

            for block in self.blocks:
                if 'html_content' in dir(block):
                    block.html_content(b, toc, subfs, product_callback)

    def write_html(self, outputp):
        report.write_html(outputp, self.name, self.html_content)
