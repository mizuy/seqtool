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

class SequencingAnalysis(object):
    def __init__(self):
        self.tc = TemplateCandidate()
        self._seqfiles = []

    def load_fasta(self, fasta_file):
        self.tc.load_fasta(fasta_file)

    def load_seqfile(self, seqfile):
        target = ''.join(i.strip() for i in open(seqfile,'rU').readlines())
        self._seqfiles.append((seqfile, target))

    def svg_seq(self, name, seq_result):
        render = BaseseqRenderer(to_seq(seq_result), False)
        for primer in self.tc.primers:
            render.add_primer(primer)

        for name, template in self.tc:
            al = sw.Alignment(seq_result, template)
            p,q = al.aseq0.location
            #print name, al.aseq0.location, al.aseq1.location
            tempname = '{} {} {}'.format(name, al.aseq1.location, al.score_text())
            comp,gap = al.compare_bar(1)
            print name
            print al.text_local()
            print al.aseq0.mid
            print al.aseq1.mid

            if al.score_density() < SCORE_THRESHOLD:
                continue
            render.add_alignment(tempname, p, q, [comp, gap])

        return render.track(800).svg()

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
                b.h2('Sequencing Analysis:')
                for name, seq in self._seqfiles:
                    filename = os.path.basename(name)+'.svg'
                    link = subfs.get_link_path(filename)
                    b.h3(name)
                    with b.a(href=link):
                        b.img(src=link, width='1000px')

                    svg = self.svg_seq(filename, seq)
                    subfs.write(filename, svg)

        subfs.finish()
