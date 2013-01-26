from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from collections import defaultdict
from cStringIO import StringIO

import os

from .memoize import memoize
from .nucleotide import bisulfite
from .pcr import Primer, PCR, primers_write_html
from .parser import parse_file, SettingFile
from .primers import load_primer_list_file
from .prompt import prompt
from . import seqsvg
from . import xmlwriter

from .subfs import SubFileSystem
from .dirutils import Filepath
from . import db
from .db.dbtss import TissuesetLocus

from .bisulfite_sequencing import BisulfiteSequencingResult
from .namedlist import NamedList, DefaultNamedList

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

def all_primers(pcrs):
    ret = set()
    for p in pcrs:
        if p.products:
            ret.add(p.fw)
            ret.add(p.rv)
    return ret

class BaseTemplate(object):
    def __init__(self):
        pass
    
    @property
    def name(self):
        return ""

    @property
    def description(self):
        return ""

    @property
    def seq(self):
        raise NotImplementedError();

    @property
    @memoize
    def bs_met(self):
        return bisulfite(self.seq, True)

    @property
    @memoize
    def bs_unmet(self):
        return bisulfite(self.seq, False)

    @property
    def transcripts(self):
        return []

    @property
    @memoize
    def transcript_start_site(self):
        return 0
        
    @property
    def features(self):
        return []

class SequenceTemplate(BaseTemplate):
    def __init__(self, sequence):
        self.seq_ = sequence

    @property
    def seq(self):
        return self.seq_;

class GenomicTemplate(BaseTemplate):
    def __init__(self, genbank_content):
        #with prompt('loading genbank...', '...done.'):
        self.genbank_ = SeqIO.read(StringIO(genbank_content), "genbank")

    @property
    def genbank(self):
        return self.genbank_

    @property
    def name(self):
        return self.name_

    @property
    def description(self):
        return self.genbank_.description

    @property
    def seq(self):
        return self.genbank_.seq

    @property
    def transcripts(self):
        for f in self.features:
            if f.type=='mRNA':
                yield f, f.extract(self.genbank_.seq)
    
    @property
    @memoize
    def transcript_start_site(self):
        v = [int(t.location.start) for t in self.features if t.type=='mRNA']
        return min(v) if v else 0
        
    @property
    def features(self):
        return self.genbank_.features

class GenebankEntry(object):
    def __init__(self, name=None):
        self.name = name
        self.template = None
        self.loaded_ = False
        self.locus = None

    def load_sequence(self, sequence):
        self.locus = None
        self.template = SequenceTemplate(sequence)

    def load_genbank(self, content, locus=None):
        self.locus = None
        self.template = GenomicTemplate(content)

    def load_gene(self, gene_id, upstream=1000, downstream=1000):
        self.locus = db.get_gene_locus(gene_id).expand(upstream,downstream)
        self.template = GenomicTemplate(db.get_locus_genbank(self.locus))

    def track_genome(self):
        assert(self.template)

        scale = 1
        length = len(self.template.seq)
        if length > 1000:
            scale = length/1000

        t = seqsvg.SeqviewTrack(scale)
        t.add_padding(10)
        start = -1* self.template.transcript_start_site
        t.add_sequence_track(self.template.seq, self.template.features, start)

        return t

    def has_transcripts(self):
        return not not self.template.transcripts

    def write_html(self, b, subfs):
        """
        subfs must have 2 methods
        def write(self, filename, content_text)
        def get_link_path(self, filename)
        """
        genome_n = 'genome.svg'

        # writing svgs
        subfs.write(genome_n, self.track_genome().svg())

        # link path for svg files
        genome_l = subfs.get_link_path(genome_n)

        # writing html
        b.h1(self.name)
        with b.div(**{'class':'images'}):
            b.h2(self.template.description)
            with b.a(href=genome_l):
                b.img(src=genome_l, width='1000px')

class GenebankTssEntry(GenebankEntry):
    def __init__(self, name=None):
        self.tss = None
        self.tsl = None

        self.tsss = []
        self.tss_name_counter = 1

        super(GenebankTssEntry, self).__init__(name)
    
    def set_tissueset(self, tissues):
        if not self.locus:
            print "No Locus Defined."
            return False
        self.tsl = TissuesetLocus(tissues, self.locus)
        return True

    def add_default_tss(self, tss_name="Assumed TSS"):
        p = self.template.transcript_start_site
        start,end = p-200, p+200
        self.tsss.append((tss_name, start, end))

    def add_tss(self, start, end, name=None):
        assert(not not start and not not end)
        if not name:
            name = "%s TSS No. %s" % (self.name, self.tss_name_counter)
            self.tss_name_counter += 1
        self.tsss.append((name, start, end))
    
    def tss_count_csv(self):
        content = ''
        for name, start, stop in self.tsss:
            content += ','.join([name] + self.tsl.count_tags(start,stop)) + '\n'
        return content

    def track_genome(self):
        assert(self.template)
        length = len(self.template.seq)

        t = super(GenebankTssEntry, self).track_genome()

        if self.tsl:
            for r in self.tsl:
                t.add_hline(length)
                t.add_dbtss_track(r, self.tsl.maxtag, self.template.seq)
            t.add_hline(length)

        return t

    def write_tss_count_csv(self, subfs):
        content = ''
        content += ', '.join(['range \\ tissue']+[t.name for t in self.tss])
        content += '\n'
        for name,start,end in self.tss_count:
            content += ', '.join([name]+[str(t.count_range(start,end)) for t in self.tss])
            content += '\n'
        subfs.write('tss.csv', content)

    def write_html(self, b, subfs):
        genome_n = 'genome.svg'

        # writing svgs
        subfs.write(genome_n, self.track_genome().svg())

        # link path for svg files
        genome_l = subfs.get_link_path(genome_n)

        # writing html
        b.h1(self.name)
        with b.div(**{'class':'images'}):
            with b.a(href=genome_l):
                b.img(src=genome_l, width='1000px')


class SeqvFileEntry(GenebankTssEntry):
    def __init__(self, name=None):
        self.primers = NamedList()
        self.motifs = NamedList()
        self.pcrs = NamedList()
        self.bs_pcrs = NamedList()
        self.rt_pcrs = NamedList()
        self.bsas = DefaultNamedList(BisulfiteSequencingResult)
        self.bsa_list = None

        super(SeqvFileEntry, self).__init__(name)

    def load_primers(self, filename):
        #with prompt('loading primers: '+filename) as pr:
        with open(filename,'r') as f:
            for i in load_primer_list_file(f):
                self.primers.append(i)
                #pr.progress()

    def set_bsa_list(self, bsa_list):
        self.bsa_list = bsa_list

    def add_primer(self, primer):
        self.primers.append(primer)

    def add_motif(self, name, seq):
        seq.name = name
        self.motifs.append(seq)

    def get_primer(self, name, default_name):
        try:
            return self.primers[name]
        except KeyError:
            if all(n.upper() in 'ATGC' for n in name):
                p = Primer(default_name, Seq.Seq(name,IUPAC.unambiguous_dna))
                self.add_primer(p)
                return p
            elif all(n.upper() in ('ATGC'+IUPAC.ambiguous_dna.letters) for n in name):
                p = Primer(default_name, Seq.Seq(name,IUPAC.ambiguous_dna))
                self.add_primer(p)
                return p
            else:
                raise KeyError("no such primer: %s"%name)
                
    def get_primers(self, fw_primer, rv_primer, pcr_name):
        fw = self.get_primer(fw_primer, 'PCR-FW(%s)'%pcr_name)
        rv = self.get_primer(rv_primer, 'PCR-RV(%s)'%pcr_name)
        return fw, rv

    def get_pcr(self, name, fw_primer, rv_primer):
        assert(self.template)
        fw, rv = self.get_primers(fw_primer, rv_primer, name)
        return PCR(name, self.template.seq, fw, rv)

    def add_pcr(self, name, fw_primer, rv_primer):
        self.pcrs.append(self.get_pcr(name,fw_primer,rv_primer))

    def add_bs_pcr(self, name, fw_primer, rv_primer):
        self.bs_pcrs.append(self.get_pcr(name,fw_primer,rv_primer))

    def add_rt_pcr(self, name, fw_primer, rv_primer):
        self.rt_pcrs.append(self.get_pcr(name,fw_primer,rv_primer))

    def add_bsa(self, cellline, pcr_name, result):
        assert(self.template)
        assert isinstance(cellline,str)
        assert isinstance(pcr_name,str)
        assert isinstance(result,str)
        try:
            pcr = self.bs_pcrs[pcr_name]
        except KeyError:
            raise ValueError('no such pcr: %s'%pcr_name)

        self.bsas[cellline].add_bsa_result(pcr, result)

    def track_genome(self):
        assert(self.template)
        length = len(self.template.seq)

        t = super(SeqvFileEntry, self).track_genome()

        if self.bsa_list:
            for name in self.bsa_list:
                bsa_map, start, end = self.bsas[name].get_map()
                t.add_bsa_track(name, bsa_map, start, end)
        else:
            for name, bsa in self.bsas.items():
                bsa_map, start, end = bsa.get_map()
                t.add_bsa_track(name, bsa_map, start, end)

        #for m in self.motifs:
            #t.add_hline(length)
            #t.add_primer_track(m.name, Primer(m.name,m), self.template.seq))

        t.add_hline(length, 10)
        t.add_pcrs_track('Genome PCR', self.pcrs)

        t.add_hline(length, 10)
        t.add_bs_pcrs_track('Bisulfite PCR', self.bs_pcrs)

        t.add_hline(length, 10)
        t.add_pcrs_track('RT-PCR', self.rt_pcrs)

        return t

    def track_transcript(self):
        assert(self.template)

        t = seqsvg.SeqviewTrack(1)

        for feature, seq in self.template.transcripts:
            length = len(seq)
            name = feature.qualifiers['product'][0]

            t.add_transcript_track(name, seq, feature)

            pcrs = [PCR(pcr.name, seq, pcr.fw, pcr.rv) for pcr in self.rt_pcrs]
            t.add_pcrs_track('RT-PCR', pcrs)
        
        return t

    def write_html(self, b, subfs):
        """
        subfs must have 2 methods
        def write(self, filename, content_text)
        def get_link_path(self, filename)
        """
        genome_n = 'genome.svg'
        transcript_n = 'transcript.svg'

        # writing svgs
        subfs.write(genome_n, self.track_genome().svg())
        if self.has_transcripts():
            subfs.write(transcript_n, self.track_transcript().svg())

        # link path for svg files
        genome_l = subfs.get_link_path(genome_n)
        if self.has_transcripts():
            transcript_l = subfs.get_link_path(transcript_n)

        # writing html
        b.h1(self.template.description)
        with b.div(**{'class':'images'}):
            b.h2('images')
            with b.div:
                b.h3('genome overview')
                with b.a(href=genome_l):
                    #b.write_raw(self.track_genome().svg_node())
                    b.img(src=genome_l, width='1000px')
                if self.has_transcripts():
                    b.h3('transcript overview')
                    with b.a(href=transcript_l):
                        b.img(src=transcript_l,width='1000px')

        with b.div(**{'class':'primers'}):
            b.h2('Primers')
            primers_write_html(b.get_writer(), self.primers)

        def write_products(products):
            if len(products) > 0:
                for c in products:
                    c.write_html(b.get_writer())
            else:
                b.p('no products')

        with b.div(**{'class':'pcrs'}):
            b.h2('PCRs')
            for pcr in self.pcrs:
                with b.div(**{'class':'pcr'}):
                    b.h3(pcr.name + " (Genome PCR)")
                    pcr.primers.write_html(b.get_writer())

                    b.h4('products')
                    with b.div(**{'class':'products'}):
                        write_products(pcr.products)

            for pcr in self.bs_pcrs:
                with b.div(**{'class':'pcr'}):
                    b.h3(pcr.name + " (Bisulfite PCR)")
                    pcr.primers.write_html(b.get_writer())

                    b.h4('products')
                    with b.div(**{'class':'products'}):
                        b.h5('template = Bisulfite-Treated (Methyl)')
                        write_products(pcr.bs_met_products)
                        b.h5('template = Bisulfite-Treated (Unmethyl)')
                        write_products(pcr.bs_unmet_products)
                        b.h5('template = Genome')
                        write_products(pcr.products)

            if self.has_transcripts():
                for pcr in self.rt_pcrs:
                    with b.div(**{'class':'pcr'}):
                        b.h3(pcr.name + " (RT PCR)")
                        pcr.primers.write_html(b.get_writer())

                        b.h4('products')
                        with b.div(**{'class':'products'}):
                            for feature, seq in self.template.transcripts:
                                length = len(seq)
                                name = feature.qualifiers['product'][0]
                                b.h5('template = transcripts: %s'%name)
                                write_products(PCR(pcr.name, seq, pcr.fw, pcr.rv).products)

                            b.h5('template = Genome')
                            write_products(pcr.products)


class SeqvFile(object):
    def __init__(self):
        self.entries = []

    def load_genbankentry(self, genbankentry):
        self.entries.append(genbankentry)

    def load_seqvfileentry(self, filename):
        inputp = Filepath(filename)
        relative_path = lambda x: inputp.relative(x)

        tss_tissues = []

        with open(filename,'r') as fileobj:
            for category, name, value, em in self.parse(fileobj):
                if not name and not value:
                    if category == 'general':
                        self.entries.append(SeqvFileEntry())
                else:
                    if not self.entries:
                        raise ValueError('no entries')

                    e = self.entries[-1]
                    if category == 'general':
                        if name=='genbank':
                            with open(relative_path(value), 'r') as f:
                                e.load_genbank(f.read())
                        elif name=='sequence':
                            e.load_sequence(Seq.Seq(value.upper(), IUPAC.unambiguous_dna))
                        elif name=='gene':
                            try:
                                gene_id, gene_symbol = db.get_gene_from_text(value)
                            except db.NoSuchGene,e:
                                em('gene entry: No such Gene %s'%value)
                                continue
                            e.load_gene(gene_id)
                        elif name=='primers':
                            e.load_primers(relative_path(value))
                        elif name=='tss':
                            e.set_tissueset([x.strip() for x in value.split(',')])
                        elif name=='bsa':
                            e.set_bsa_list([x.strip() for x in value.split(',')])
                    elif category == 'primer':
                        e.add_primer(Primer(name, Seq.Seq(value.upper(),IUPAC.ambiguous_dna)))
                    elif category == 'motif':
                        e.add_motif(name,Seq.Seq(value,IUPAC.ambiguous_dna))
                    elif category in ['pcr','rt_pcr','bs_pcr']:
                        name = name.split(',')[0].strip()
                        ls = value.split(',')
                        if len(ls)!=2:
                            em('you must specify 2 primer names separated by "," for each pcr: %s'%name)
                            continue
                        fw = ls[0].strip()
                        rv = ls[1].strip()

                        if category == 'pcr':
                            e.add_pcr(name, fw, rv)
                        elif category == 'rt_pcr':
                            e.add_rt_pcr(name, fw, rv)
                        elif category == 'bs_pcr':
                            e.add_bs_pcr(name, fw, rv)
                        else:
                            raise ValueError
                    elif category == 'bsa':
                        n = [n.strip() for n in name.split(',')]
                        if not len(n)>=2:
                            em('each bsa must have at least 2 key; cell line name and pcr name')
                            continue
                        pcrname = n[0].strip()
                        cellline = n[1].strip().upper()
                        annotations = n[2:]
                        if not pcrname or not cellline:
                            em('empty pcr or cellline name: %s, %s'%(pcrname,cellline))
                            continue

                        e.add_bsa(cellline, pcrname, value.strip().upper())
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


    def write_html(self, outputp):
        with open(outputp.path,'w') as f:
            html = xmlwriter.XmlWriter(f)
            b = xmlwriter.builder(html)
            with b.html:
                with b.head:
                    with b.style(type='text/css'):
                        b.text(seqview_css)
            with b.body:
                for e in self.entries:
                    subfs = SubFileSystem(outputp.dir, outputp.prefix)

                    e.write_html(b, subfs)

                    subfs.finish()

