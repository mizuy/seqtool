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
from .db import entrez
from . import seqtrack
from . import xmlwriter
from .listdict import ListDict
from .subfs import SubFileSystem
from .dirutils import Filepath
from .db.dbtss import RegionTissueset

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

class BisulfiteSequenceEntry(object):
    def __init__(self, name, pcr, result):
        self.name = name
        self.pcr = pcr
        self.result = result

        p = pcr.bs_met_products
        if not len(p)==1:
            raise ValueError('number of pcr products of %s must be 1 but %s'%(pcr_name,len(products)))
        self.product = p[0]
        if not all(i in 'MUP?' for i in result):
            raise ValueError('bsp result must contain only M,U,P or ?')

        self.cpg = self.product.detectable_cpg()
        if len(result)!=self.cpg:
            raise ValueError('%s has %s detectable CpG sites, but result gives %s'%(pcr_name,self.cpg,len(result)))

        self.cpg_sites = self.product.cpg_sites()

        self.bsp_map = [(n,result[i]) for i,n in enumerate(self.cpg_sites)]

class BisulfiteSequenceEntries(object):
    def __init__(self):
        self.entries = []
    def append(self, entry):
        self.entries.append(entry)

    def combine(self):
        start_i = min(e.product.start_i for e in self.entries)
        end_i = max(e.product.end_i for e in self.entries)

        results = defaultdict(str)
        for e in self.entries:
            for i,n in enumerate(e.cpg_sites):
                if e.result[i]!='?':
                    results[n] += e.result[i]

        bsp_map = []
        for index, result in results.items():
            if all(r in 'M' for r in result):
                c = 'M'
            elif all(r in 'U' for r in result):
                c = 'U'
            else:
                c = 'P'
            bsp_map.append( (index, c) )
            
        return bsp_map, start_i, end_i

class BisulfiteSequence(object):
    def __init__(self):
        self.keys = []
        self.entries = defaultdict(BisulfiteSequenceEntries)
    def __iter__(self):
        for k in self.keys:
            yield k, self.entries[k]

    def __getitem__(self, name):
        if name not in self.keys:
            raise KeyError
        return self.entries[name]

    def add(self, name, entry):
        if name not in self.keys:
            self.keys.append(name)
        self.entries[name].append(entry)

class GenomicTemplate(object):
    def __init__(self, genbank_content):
        with prompt('loading genbank...', '...done.'):
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
    @memoize
    def bs_met(self):
        return bisulfite(self.genbank_.seq, True)

    @property
    @memoize
    def bs_unmet(self):
        return bisulfite(self.genbank_.seq, False)

    @property
    def transcripts(self):
        for f in self.features:
            if f.type=='mRNA':
                yield f, f.extract(self.genbank_.seq)
    
    @property
    @memoize
    def transcript_start_site(self):
        v = [int(t.location.start) for t in self.features if t.type=='mRNA']
        if v:
            i = min(v)
        else:
            i = 0
        return i
        
    @property
    def features(self):
        return self.genbank_.features

class GeneBankEntry(object):
    def __init__(self, name=None):
        self.name_ = name
        self.template = None
        self.loaded_ = False

    def load_genbank(self, content):
        self.template = GenomicTemplate(content)

    def render_genome(self, width):
        assert(self.template)
        length = len(self.template.seq)

        t = seqtrack.TrackGroup()
        t.add(seqtrack.Track(0,10))
        t.add(seqtrack.SequenceTrack(self.template.seq, self.template.features, -1* self.template.transcript_start_site))

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
        b.h1(self.template.description)
        with b.div(**{'class':'images'}):
            b.h2('images')
            with b.div:
                b.h3('genome overview')
                with b.a(href=genome_l):
                    b.img(src=genome_l, width='1000px')

class SeqvFileEntry(GeneBankEntry):
    def __init__(self, name=None):
        self.primers = ListDict()
        self.motifs = ListDict()

        # RegionTissueSet * tsss
        self.rts = None
        self.tsss = []
        self.tss_name_counter = 1

        self.pcrs = ListDict()
        self.bs_pcrs = ListDict()
        self.rt_pcrs = ListDict()

        self.bsps = BisulfiteSequence()

        super(SeqvFileEntry, self).__init__(name)

    def add_default_tss(self, tss_name="Assumed TSS"):
        p = self.template.transcript_start_site
        start,end = p-200, p+200
        self.tsss.append((tss_name, start, end))

    def add_tss_tissues(self, tissues, locus):
        self.rts = RegionTissueset(tissues, locus)

    def add_tss(self, start, end, name=None):
        assert(not not start and not not end)
        if not name:
            name = "%s TSS No. %s" % (self.name, self.tss_name_counter)
            self.tss_name_counter += 1
        self.tsss.append((name, start, end))
    
    def tss_count_csv(self):
        content = ''
        for name, start, stop in self.tsss:
            content += ','.join([name] + self.rts.count_tags(start,stop)) + '\n'
        return content

    def load_primers(self, filename):
        with prompt('loading primers: '+filename) as pr:
            with open(filename,'r') as f:
                for i in load_primer_list_file(f):
                    self.primers.append(i)
                    pr.progress()


    def add_primer(self, primer):
        self.primers.append(primer)

    def add_motif(self, name, seq):
        seq.name = name
        self.motifs.append(seq)

    def get_primer(self, name, pcr_name):
        try:
            return self.primers[name]
        except KeyError:
            if all(n.upper() in 'ATGC' for n in name):
                p = Primer(pcr_name, Seq.Seq(name,IUPAC.unambiguous_dna))
                self.add_primer(p)
                return p
            else:
                raise KeyError('no such forward primer: %s'%name)

    def get_primers(self, fw_primer, rv_primer, name):
        fw = self.get_primer(fw_primer, 'PCR-FW(%s)'%name)
        rv = self.get_primer(rv_primer, 'PCR-RV(%s)'%name)
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

    def add_bsp(self, cellline, pcr_name, result):
        assert(self.template)
        assert isinstance(cellline,str)
        assert isinstance(pcr_name,str)
        assert isinstance(result,str)
        try:
            pcr = self.bs_pcrs[pcr_name]
        except KeyError:
            raise ValueError('no such pcr: %s'%pcr_name)

        self.bsps.add(cellline, BisulfiteSequenceEntry(cellline, pcr, result))

    def render_genome(self, width):
        assert(self.template)
        length = len(self.template.seq)

        t = seqtrack.TrackGroup()
        t.add(seqtrack.Track(0,10))
        t.add(seqtrack.SequenceTrack(self.template.seq, self.template.features, -1* self.template.transcript_start_site))

        if self.rts:
            for r in self.rts:
                t.add(seqtrack.HbarTrack('', length))
                t.add(seqtrack.DbtssTrack(r, self.rts.maxtag, self.template.seq))
            t.add(seqtrack.HbarTrack('', length))

        for name, bsp in self.bsps:
            bsp_map, start, end = bsp.combine()
            t.add(seqtrack.BisulfiteSequenceTrack(name, bsp_map, start, end))

        for m in self.motifs:
            t.add(seqtrack.HbarTrack('', length))
            t.add(seqtrack.PrimerTrack(m.name, Primer(m.name,m), self.template.seq))

        t.add(seqtrack.HbarTrack('', length))
        t.add(seqtrack.PcrsTrack('genome', self.pcrs))

        t.add(seqtrack.HbarTrack('', length))
        t.add(seqtrack.BsPcrsTrack('bisulfite pcr', self.bs_pcrs))

        t.add(seqtrack.HbarTrack('', length))
        t.add(seqtrack.PcrsTrack('rt pcr', self.rt_pcrs))

        return t.svg(width)

    def render_transcript(self, width):
        assert(self.template)

        t = seqtrack.TrackGroup()

        for feature, seq in self.template.transcripts:
            length = len(seq)
            name = feature.qualifiers['product'][0]

            t.add(seqtrack.TranscriptTrack(name, seq, feature))

            pcrs = [PCR(pcr.name, seq, pcr.fw, pcr.rv) for pcr in self.rt_pcrs]
            t.add(seqtrack.PcrsTrack('RT-PCR', pcrs))
        
        return t.svg(width)

    def render_bsp(self, width=None, scale=2, init=50):
        assert(self.template)
        return ''
        '''
        TODO:

        if not self.record:
            raise ValueError('load genbank first')
        if not self.has_bisulfite_sequence_result():
            raise ValueError('no bisulfite sequence result')

        if not self.cpg_location:
            self.cpg_location = self._calc_cpg_location()

        count = len(self.cpg_location)
        def cx(x):
            return init+x*scale
        w = SeqCanvasOverview(0, cx(count))

        color = {'M':'#F00','U':'#00F','P':'#AA0'}
        
        y = 1
        for cl in self.celllines:
            w.put_line(init, cx(count), y)
            w.put_char(Point(0, y-0.5), cl)
            for bs in self.bisulfite_sequence[cl]:
                for x,i in enumerate(bs.product.cpg_sites()):
                    loc = cx(self.cpg_location[i])
                    result = bs.result[x]
                    if not result=='?':
                        w.put_yline(loc,y-.5,y,color[result])
            y += 1

        w.save(filename, format='PNG',height=y)
        return True
        '''

    def write_tss_count_csv(self, subfs):
        content = ''
        content += ', '.join(['range \\ tissue']+[t.name for t in self.tss])
        content += '\n'
        for name,start,end in self.tss_count:
            content += ', '.join([name]+[str(t.count_range(start,end)) for t in self.tss])
            content += '\n'
        subfs.write('tss.csv', content)

    def write_html(self, b, subfs):
        """
        subfs must have 2 methods
        def write(self, filename, content_text)
        def get_link_path(self, filename)
        """
        genome_n = 'genome.svg'
        transcript_n = 'transcript.svg'
        bsp_n = 'bsp.svg'

        # writing svgs
        subfs.write(genome_n, self.render_genome(5000))
        subfs.write(transcript_n, self.render_transcript(1200))
        subfs.write(bsp_n, self.render_bsp())

        # link path for svg files
        genome_l = subfs.get_link_path(genome_n)
        transcript_l = subfs.get_link_path(transcript_n)
        bsp_l = subfs.get_link_path(bsp_n)

        # writing html
        b.h1(self.template.description)
        with b.div(**{'class':'images'}):
            b.h2('images')
            with b.div:
                b.h3('genome overview')
                with b.a(href=genome_l):
                    b.img(src=genome_l, width='1000px')
                b.h3('transcript overview')
                with b.a(href=transcript_l):
                    b.img(src=transcript_l,width='1000px')
                b.h3('bsp overview')
                b.text('Not Implemented Yet')
                #with b.a(href=bsp_l):
                #b.img(src=bsp_l,width='1000px')

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
                        elif name=='primers':
                            e.load_primers(relative_path(value))

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
                    elif category == 'bsp':
                        n = [n.strip() for n in name.split(',')]
                        if not len(n)>=2:
                            em('each bsp must have at least 2 key; cell line name and pcr name')
                            continue
                        pcrname = n[0].strip()
                        cellline = n[1].strip().upper()
                        annotations = n[2:]
                        if not pcrname or not cellline:
                            em('empty pcr or cellline name: %s, %s'%(pcrname,cellline))
                            continue

                        e.add_bsp(cellline, pcrname, value.strip().upper())
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
                    subfs = SubFileSystem(outputp.dir, outputp.suffix)

                    e.write_html(b, subfs)

                    subfs.finish()

