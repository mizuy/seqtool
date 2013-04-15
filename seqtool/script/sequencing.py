from Bio import SeqIO

import os
from seqtool.nucleotide import to_seq
from seqtool.nucleotide.cpg import bisulfite_conversion_ambiguous, c2t_conversion
from seqtool.nucleotide import sw
#from Bio import SeqIO
from seqtool.view.baseseq_renderer import BaseseqRenderer
from ..util.subfs import SubFileSystem
from ..util import xmlwriter
from ..view.css import seqview_css
from ..nucleotide.pcr import Primer

SCORE_THRESHOLD = 1.5
# TODO: .scf file support (4peaks output)
# cf.
# http://staden.sourceforge.net/manual/formats_unix_3.html#SEC3
# format
# http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-pysrc.html
# abi = SeqIO.read(filename, "abi")

class TemplateCandidate(object):
    def __init__(self):
        self._templates = []
        self.primers = []

    def add_template(self, name, seq):
        self._templates.append((name,seq))

    def add_bisulfite_template(self, name, seq):
        self.add_template(name, seq)
        self.add_template(name+' BS+', bisulfite_conversion_ambiguous(seq, True))
        self.add_template(name+' BS-', bisulfite_conversion_ambiguous(seq, False))

    def add_c2t_template(self, name, seq):
        self.add_template(name, seq)
        self.add_template(name+' C2T+', c2t_conversion(seq, True))
        self.add_template(name+' C2T-', c2t_conversion(seq, False))

    def __iter__(self):
        for name, seq in self._templates:
            for sense in [True,False]:
                if sense:
                    s = seq
                else:
                    s = seq.reverse_complement()

                yield name+' '+('sense' if sense else 'asense'), s

    def load_fasta(self, filename):
        for record in SeqIO.parse(open(filename,'r'), "fasta"):
            name, sep, conv = record.description.partition(':')
            conv = conv.strip()

            #print name, conv

            if conv=='C2T':
                self.add_c2t_template(name, record.seq)
            elif conv=='BS':
                self.add_bisulfite_template(name, record.seq)
            elif conv=='Primer':
                self.primers.append(Primer(name, record.seq))
            else:
                self.add_template(name, record.seq)

    def alignments(self, target):
        ret = []
        for name, template in self:
            al = sw.Alignment(target, template)
            p,q = al.aseq0.location
            tempname = '{} {} {}'.format(name, al.aseq1.location, al.score_text())
            ret.append((al, p, q, tempname))
        return ret

class SeqFile(object):
    def __init__(self, name, filename, template_candidate):
        self.name = name
        self.filename = filename
        self.tc = template_candidate
        self.seq = ''.join(i.strip() for i in open(filename,'rU').readlines())
        self._alignment = self.tc.alignments(self.seq)

    def svg(self):
        render = BaseseqRenderer(to_seq(self.seq), False)
        for primer in self.tc.primers:
            render.add_primer(primer)

        for al, p, q, tempname in self._alignment:
            comp,gap = al.compare_bar(1)

            if al.score_density() < SCORE_THRESHOLD:
                continue
            render.add_alignment(tempname, p, q, [comp, gap])

        return render.track(len(self.seq)+10).svg()

    def html(self, b):
        #b.h3(self.name)
        for al, p, q, tempname in self._alignment:
            b.h4(tempname)
            if al.score_density() < SCORE_THRESHOLD:
                continue
            with b.pre:
                b.text(al.text_local())


class SequencingAnalysis(object):
    def __init__(self):
        self.tc = TemplateCandidate()
        self._seqfiles = []

    def load_fasta(self, fasta_file):
        self.tc.load_fasta(fasta_file)

    def load_seqfile(self, seqfile):
        self._seqfiles.append(SeqFile(os.path.basename(seqfile), seqfile, self.tc))

    def write_html(self, outputp):
        with SubFileSystem(outputp.dir, outputp.prefix) as subfs:

            with open(outputp.path,'w') as output:
                html = xmlwriter.XmlWriter(output)
                b = xmlwriter.builder(html)
                with b.html:
                    with b.head:
                        with b.style(type='text/css'):
                            b.text(seqview_css)
                with b.body:
                    b.h2('Sequencing Analysis:')
                    b.h3('Alignment View')
                    for sf in self._seqfiles:
                        filename = sf.name+'.svg'
                        subfs.write(filename, sf.svg())
                        link = subfs.get_link_path(filename)

                        b.h3(sf.name)
                        with b.div(cls='products'):
                            with b.a(href=link):
                                b.img(src=link, width='1000px')
                            sf.html(b)

