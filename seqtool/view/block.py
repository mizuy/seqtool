

from ..nucleotide.pcr import PCR
from ..nucleotide.primer import Primer
from ..util.namedlist import NamedList

class PcrsHolder(object):
    def __init__(self, template, primers):
        self.pcrs = NamedList()
        self.primers = primers
        self.template = template

    def __len__(self):
        return len(self.pcrs)

    def __iter__(self):
        return iter(self.pcrs)

    def _get(self, name, fw_primer, rv_primer):
        fw = self.primers.get_default(fw_primer, 'PCR-FW(%s)'%name)
        rv = self.primers.get_default(rv_primer, 'PCR-RV(%s)'%name)
        return PCR(name, self.template.seq, fw, rv)

    def add(self, name, fw_primer, rv_primer):
        pcr = self._get(name,fw_primer,rv_primer)
        self.pcrs.append(pcr)
        return pcr

    def get(self, name):
        return self.pcrs.get(name)
    
    def used_primers(self, pcrs):
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
        self.pcrs = PcrsHolder(template, primers)
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
                    b.h5('template = Bisulfite-Treated Sense Strand (CpG Methylated)')
                    self.write_products(b, pcr.bs_products(methyl=True,sense=True))
                    b.h5('template = Bisulfite-Treated Sense Strand (CpG Unmethylated)')
                    self.write_products(b, pcr.bs_products(methyl=False,sense=True))
                    b.h5('template = Bisulfite-Treated Antisense Strand (CpG Methylated)')
                    self.write_products(b, pcr.bs_products(methyl=True,sense=False))
                    b.h5('template = Bisulfite-Treated Antisense Strand (CpG Unmethylated)')
                    self.write_products(b, pcr.bs_products(methyl=False,sense=False))
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
