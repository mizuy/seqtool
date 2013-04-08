from __future__ import absolute_import

from ..nucleotide.pcr import Primer, PCR, load_primer_list_file

from ..util.subfs import SubFileSystem
from ..db.dbtss import TissuesetLocus

from .bisulfite_sequencing import BisulfiteSequencingResult
from ..util.namedlist import NamedList, DefaultNamedList

class Primers(object):
    def __init__(self):
        self.primers = NamedList()

    def __len__(self):
        return len(self.primers)

    def __iter__(self):
        return iter(self.primers)

    def load(self, filename):
        with open(filename,'r') as f:
            for i in load_primer_list_file(f):
                self.primers.append(i)

    def add(self, primer):
        self.primers.add(primer)

    def get(self, name, default_name):
        try:
            return self.primers[name]
        except KeyError:
            try:
                return Primer(default_name, name)
            except:
                raise KeyError("no such primer: %s"%name)
                
    def get_pcr_pair(self, fw_primer, rv_primer, pcr_name):
        fw = self.get(fw_primer, 'PCR-FW(%s)'%pcr_name)
        rv = self.get(rv_primer, 'PCR-RV(%s)'%pcr_name)
        return fw, rv

class Pcrs(object):
    def __init__(self, template, primers):
        self.pcrs = NamedList()
        self.primers = primers
        self.template = template

    def __len__(self):
        return len(self.pcrs)

    def __iter__(self):
        return iter(self.pcrs)

    def _get(self, name, fw_primer, rv_primer):
        fw, rv = self.primers.get_pcr_pair(fw_primer, rv_primer, name)
        return PCR(name, self.template.seq, fw, rv)

    def add(self, name, fw_primer, rv_primer):
        pcr = self._get(name,fw_primer,rv_primer)
        self.pcrs.append(pcr)
        return pcr

    def get(self, name):
        return self.pcrs.get(name)
    
    def used_primers(pcrs):
        ret = set()
        for p in self.pcrs:
            if p.products:
                ret.add(p.fw)
                ret.add(p.rv)
        return ret


class BaseBlock(object):
    def svg_genome(self, t):
        pass

    def svg_transcript(self, t, transcript):
        pass

    def write_html(self, b, subfs):
        pass


class PcrsBlock(BaseBlock):
    def __init__(self, template, primers):
        self.kind = 'Genome PCR'
        self.pcrs = Pcrs(template, primers)
        self.template = template

    def __iter__(self):
        return iter(self.pcrs)

    def add(self, name, fw, rv):
        self.pcrs.add(name, fw, rv)

    def svg_genome(self, t):
        t.add_pcrs_track(self.kind, self.pcrs.pcrs)

    def svg_transcripts(self, t, transcript):
        pass

    def write_products(self, b, products):
        if len(products) > 0:
            for c in products:
                c.write_html(b.get_writer())
        else:
            b.p('no products')

    def write_html(self, b, subfs):
        for pcr in self.pcrs:
            with b.div(**{'class':'pcr'}):
                b.h3(pcr.name + " (%s)"%self.kind)
                pcr.primers.write_html(b.get_writer())

                b.h4('products')
                with b.div(**{'class':'products'}):
                    self.write_products(b, pcr.products)

class BsPcrsBlock(PcrsBlock):
    """
    Bisulfite Sequence sample usually contains Genomic DNA as well.
    so, targets  are bisulite-methyl-DNA, unmethyl-DNA and Genomic DNA.
    """
    def __init__(self, template, primers):
        super(BsPcrsBlock,self).__init__(template, primers)
        self.kind = 'Bisulfite PCR'

    def svg_genome(self, t):
        t.add_bs_pcrs_track(self.kind, self.pcrs.pcrs)

    def write_html(self, b, subfs):
        for pcr in self.pcrs:
            with b.div(**{'class':'pcr'}):
                b.h3(pcr.name + " (%s)"%self.kind)
                pcr.primers.write_html(b.get_writer())

                b.h4('products')
                with b.div(**{'class':'products'}):
                    b.h5('template = Bisulfite-Treated (Methyl)')
                    self.write_products(b, pcr.bs_met_products)
                    b.h5('template = Bisulfite-Treated (Unmethyl)')
                    self.write_products(b, pcr.bs_unmet_products)
                    b.h5('template = Genome')
                    self.write_products(b, pcr.products)

class RtPcrsBlock(PcrsBlock):
    """
    RT-PCR sample usually contains Genomic DNA as well.
    so, targets of RT-PCR are all transcripts and Genomic DNA.
    """
    def __init__(self, template, primers):
        super(RtPcrsBlock,self).__init__(template, primers)
        self.kind = 'RT-PCR'

    def svg_genome(self, t):
        t.add_pcrs_track(self.kind, self.pcrs.pcrs)

    def svg_transcript(self, t, transcript):
        pcrs = [PCR(pcr.name, transcript.seq, pcr.fw, pcr.rv) for pcr in self.pcrs.pcrs]
        t.add_pcrs_track('RT-PCR', pcrs)

    def write_html(self, b, subfs):
        for pcr in self.pcrs:
            with b.div(**{'class':'pcr'}):
                b.h3(pcr.name + " (RT PCR)")
                pcr.primers.write_html(b.get_writer())

                b.h4('products')
                with b.div(**{'class':'products'}):
                    for transcript in self.template.transcripts:
                        b.h5('template = transcripts: %s'%transcript.name)
                        self.write_products(b, PCR(pcr.name, transcript.seq, pcr.fw, pcr.rv).products)

                    b.h5('template = Genome')
                    self.write_products(b, pcr.products)


class DbtssBlock(BaseBlock):
    def __init__(self, template):
        self.template = template
        self.tss = None
        self.tsl = None

        self.tsss = []
        self.tss_name_counter = 1
    
    def set_tissueset(self, tissues):
        if not self.template.locus:
            print "No Locus Defined."
            return False
        self.tsl = TissuesetLocus(tissues, self.template.locus)
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

    def svg_genome(self, t):
        assert(self.template)
        length = len(self.template)

        if self.tsl:
            for r in self.tsl:
                t.add_hline(length)
                t.add_dbtss_track(r, self.tsl.maxtag, self.template.seq)
            t.add_hline(length)

    def write_tss_count_csv(self, subfs):
        content = ''
        content += ', '.join(['range \\ tissue']+[t.name for t in self.tss])
        content += '\n'
        for name,start,end in self.tss_count:
            content += ', '.join([name]+[str(t.count_range(start,end)) for t in self.tss])
            content += '\n'
        subfs.write('tss.csv', content)


class BsaBlock(BaseBlock):
    def __init__(self, bs_pcrs):
        self.bsas = DefaultNamedList(BisulfiteSequencingResult)
        self.celllines = None
        self.bs_pcrs = bs_pcrs

    def set_celllines(self, celllines):
        self.celllines = celllines

    def add(self, cellline, pcr_name, result):
        assert isinstance(cellline,str)
        assert isinstance(pcr_name,str)
        assert isinstance(result,str)
        try:
            pcr = self.bs_pcrs.pcrs.get(pcr_name)
        except KeyError:
            raise ValueError('no such pcr: %s'%pcr_name)

        self.bsas[cellline].add_bsa_result(pcr, result)

    def svg_genome(self, t):
        if self.celllines:
            for name in self.celllines:
                bsa_map, start, end = self.bsas[name].get_map()
                t.add_bsa_track(name, bsa_map, start, end)
        else:
            for name, bsa in self.bsas.items():
                bsa_map, start, end = bsa.get_map()
                t.add_bsa_track(name, bsa_map, start, end)
