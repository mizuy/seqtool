import os
from Bio import SeqIO

from ..nucleotide import to_seq, sw
from ..nucleotide.cpg import bisulfite_conversion, c2t_conversion
from ..nucleotide.primer import Primer
from ..util import report
from ..view.baseseq_renderer import BaseseqRenderer
from .. import format
from ..format.abi import AbiFormat
from ..format.render import SvgPeaksAlignment


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
        self.add_template(name+' BS+', bisulfite_conversion(seq, True))
        self.add_template(name+' BS-', bisulfite_conversion(seq, False))

    def add_c2t_template(self, name, seq):
        self.add_template(name, seq)
        self.add_template(name+' C2T+', c2t_conversion(seq, True))
        self.add_template(name+' C2T-', c2t_conversion(seq, False))

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
        for name, template in self._templates:
            for sense in [True,False]:
                if sense:
                    name = name+'(original)'
                    st = template
                else:
                    name = name+'(reverse)'
                    st = str(to_seq(template).reverse_complement())
                al = sw.Alignment(target, st)
                p,q = al.aseq0.location
                tempname = '{} {} {}'.format(name, al.aseq1.location, al.score_text())
                ret.append((al, p, q, tempname))
        return ret

class SeqFile(object):
    def __init__(self, name, filename, template_candidate):
        self.name = name
        self.filename = filename
        self.tc = template_candidate
        base,ext = os.path.splitext(self.filename)
        if ext=='.seq':
            self.abi = None
            self.seq = ''.join(i.strip() for i in open(filename,'rU').readlines())
        elif ext=='.ab1':
            self.abi = AbiFormat(filename)
            self.seq = self.abi.get_sequence()
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

    def svg_peak(self):
        render = SvgPeaksAlignment()

        seq_loc = self.abi.get_sequence_loc()

        for al, p, q, tempname in self._alignment:
            if al.score_density() < SCORE_THRESHOLD:
                continue

            locs = list(al.get_mid_loc(seq_loc))
            render.add_text(tempname, start=locs[0])
            render.add_text_loc(al.aseq1.mid_gap, locs)
            render.add_text_loc(al.match_bar(), locs)
            render.add_text_loc(al.aseq0.mid_gap, locs)

        render.add_text_loc(self.abi.get_sequence(), seq_loc)
        render.add_peaks(200, self.abi.get_peaks(), seq_loc)

        return render.svg()

    def html(self, b):
        #b.h3(self.name)
        for al, p, q, tempname in self._alignment:
            b.h4(tempname)
            
            if al.score_density() < SCORE_THRESHOLD:
                continue

            with b.div(cls='indent'):
                m = al.correspondance_map()
                ms = al.correspondance_str()

                with b.pre:
                    for k,v in list(m.items()):
                        b.text(str(k)+':{'+', '.join('{}:{}'.format(kk,vv) for kk,vv in list(v.items()))+'}')

                for k,v in list(ms.items()):
                    with b.span:
                        b.text(str(k)+': '+v)
                        b.br()

                with b.pre:
                    b.text(al.text_local())


class SequencingAnalysis(object):
    def __init__(self):
        self.tc = TemplateCandidate()
        self._seqfiles = []
        self.name = 'Sequencing Analysis Result'

    def load_fasta(self, fasta_file):
        self.tc.load_fasta(fasta_file)

    def load_seqfile(self, seqfile):
        self._seqfiles.append(SeqFile(os.path.basename(seqfile), seqfile, self.tc))

    def write_html(self, outputp):
        report.write_html(outputp, self.name, self.html_content)

    def html_content(self, b, toc, subfs):
        b.h2('Sequencing Analysis:')
        b.h3('Alignment View')
        for sf in self._seqfiles:
            fs = [(sf.name+'.svg', sf.svg())]
            if sf.abi:
                fs.append((sf.name+'_peak.svg', sf.svg_peak()))

            b.h3(sf.name)
            with b.div(cls='products'):
                for filename, content in fs:
                    link = subfs.get_link_path(filename)
                    subfs.write(filename, content)
                    
                    with b.a(href=link):
                        b.img(src=link, width='1000px')
                sf.html(b)

