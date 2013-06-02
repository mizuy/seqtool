from ..nucleotide.pcr import PCR, BisulfitePCR, RtPCR
from ..util.namedlist import NamedList
from ..util import report

class PcrsHolder:
    def __init__(self, primers):
        self.pcrs = NamedList()
        self.primers = primers

    def __len__(self):
        return len(self.pcrs)

    def __iter__(self):
        return iter(self.pcrs)

    def get_primers(self, name, fw_primer, rv_primer):
        fw = self.primers.get_default(fw_primer, "Foward of {}".format(name))
        rv = self.primers.get_default(rv_primer, "Reverse of {}".format(name))
        return fw,rv

    def add(self, pcr):
        self.pcrs.append(pcr)

    def get(self, name):
        return self.pcrs.get(name)
    
    def used_primers(self, pcrs):
        ret = set()
        for p in self.pcrs:
            if p.products:
                ret.add(p.fw)
                ret.add(p.rv)
        return ret

class BaseBlock:
    def __init__(self, title):
        self.title = title

    def svg_genome(self, t):
        pass

    def svg_transcript(self, t, transcript):
        pass

    def write_html(self, b, subfs):
        pass

class PcrsBlock(BaseBlock):
    def __init__(self, template, primers, title='Genome PCR'):
        self.pcrsholder = PcrsHolder(primers)
        self.template = template
        super().__init__(title)

    def __iter__(self):
        return iter(self.pcrsholder)

    def add(self, name, fw, rv):
        fw, rv = self.pcrsholder.get_primers(name, fw, rv)
        self.pcrsholder.add(PCR(name, self.template.seq, fw, rv))

    def svg_genome(self, t):
        t.add_pcrs_track(self.title, self.pcrsholder.pcrs)

    def html_content(self, b, toc, subfs):
        w = b.get_writer()
        if not len(self.pcrsholder):
            return

        with report.section(b, toc, self.title):
            for pcr in self.pcrsholder:
                with report.section(b, toc, "{} ({})".format(pcr.name,self.title), klass='pcr'):
                    pcr.write_html(w)

class BsPcrsBlock(PcrsBlock):
    """
    Bisulfite Sequence sample usually contains Genomic DNA as well.
    so, targets  are bisulite-methyl-DNA, unmethyl-DNA and Genomic DNA.
    """
    def __init__(self, template, primers):
        super().__init__(template, primers, 'Bisulfite PCR')

    def add(self, name, fw, rv):
        fw, rv = self.pcrsholder.get_primers(name, fw, rv)
        self.pcrsholder.add(BisulfitePCR(name, self.template.seq, fw, rv))

    def svg_genome(self, t):
        t.add_bs_pcrs_track(self.title, self.pcrsholder.pcrs)

class RtPcrsBlock(PcrsBlock):
    """
    RT-PCR sample usually contains Genomic DNA as well.
    so, targets of RT-PCR are all transcripts and Genomic DNA.
    """
    def __init__(self, template, primers):
        super().__init__(template, primers, 'RT-PCR')

    def add(self, name, fw, rv):
        fw, rv = self.pcrsholder.get_primers(name, fw, rv)
        self.pcrsholder.add(RtPCR(name, self.template.seq, self.template.transcripts, fw, rv))

    def svg_genome(self, t):
        t.add_pcrs_track(self.title, [p.genome_pcr() for p in self.pcrsholder.pcrs])

    def svg_transcript(self, t, transcript):
        t.add_pcrs_track('RT-PCR', [p.transcript_pcr(transcript) for p in self.pcrsholder.pcrs])
