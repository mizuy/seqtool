import os
from Bio import SeqIO

from ..nucleotide import to_seq
from ..nucleotide.alignment import make_alignment
from ..nucleotide.cpg import bisulfite_conversion, c2t_conversion
from ..nucleotide.primer import Primer
from ..util import report
from ..view.baseseq_renderer import BaseseqRenderer
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
                al = make_alignment(target,st)
                p,q = al.aseq0.location
                tempname = '{} {}'.format(name, al.aseq1.location)
                ret.append((al, p, q, tempname))
        return ret

class SequencingResult(object):
    def __init__(self, name, filename, template_candidate):
        self.name = name
        self.filename = filename
        self.tc = template_candidate
        base,ext = os.path.splitext(self.filename)
        if ext=='.seq':
            self.abi = None
            self.seq = to_seq(''.join(i.strip() for i in open(filename,'rU').readlines()))
        elif ext=='.ab1':
            self.abi = AbiFormat(filename)
            self.seq = to_seq(self.abi.view.get_sequence())

    def html_content(self, b, toc, subfs):
        b.h3(self.name)
        for sense in [True,False]:
            seq = self.seq if sense else self.seq.reverse_complement()

            name = self.name+('_sense_' if sense else '_antisense_')
            b.h4('Sense' if sense else 'Reverse Complement')

            alignments = self.tc.alignments(seq)
            
            def svg_compare():
                render = BaseseqRenderer(to_seq(seq), False)
                for primer in self.tc.primers:
                    render.add_primer(primer)

                for al, p, q, tempname in alignments:
                    if al.score_density() < SCORE_THRESHOLD:
                        continue

                    comp,gap = al.compare_bar()
                    render.add_alignment(tempname, p, q, [comp, gap])

                return render.track(len(seq)+10).svg()
                
            def svg_peak():
                view = self.abi.get_view(sense)
                render = SvgPeaksAlignment()

                seq_loc = view.get_location()

                for al, p, q, tempname in alignments:
                    if al.score_density() < SCORE_THRESHOLD:
                        continue

                    u,d = al.get_common_first_last_length()
                    tlocs = list(al.get_loc(seq_loc,u,d))
                    render.add_text(tempname, start=tlocs[0])
                    render.add_text_loc(al.aseq1.local(u,d), tlocs)
                    
                    mlocs = list(al.get_loc(seq_loc,0,0))
                    render.add_text_loc(al.match_bar(), mlocs)
                    render.add_text_loc(al.aseq0.mid_gap, mlocs)

                render.add_text_loc(view.get_sequence(), seq_loc)
                render.add_peaks(200, view.get_peaks(), seq_loc)

                return render.svg()
                

            fs = [(name+'_compare.svg', svg_compare())]
            if self.abi:
                fs.append((name+'_peak.svg', svg_peak()))

            with b.div(klass='products'):
                with b.ul:
                    for filename, content in fs:
                        link = subfs.get_link_path(filename)
                        subfs.write(filename, content)

                        with b.li:
                            with b.a(href=link):
                                b.text(filename)

                for al, p, q, tempname in alignments:
                    if al.score_density() < SCORE_THRESHOLD:
                        continue

                    b.h4(tempname)
                    with b.div(cls='indent'):
                        cm = al.reversed().correspondance_map()

                        with b.pre:
                            #b.text(cm.text_map())
                            b.text(cm.text_str())

                        with b.pre:
                            b.text(al.text_local())


class SequencingAnalysis(object):
    def __init__(self):
        self.tc = TemplateCandidate()
        self._results = []
        self.name = 'Sequencing Analysis'

    def load_fasta(self, fasta_file):
        self.tc.load_fasta(fasta_file)

    def load_sequencingfile(self, seqfile):
        self._results.append(SequencingResult(os.path.basename(seqfile), seqfile, self.tc))

    def write_html(self, outputp):
        report.write_html(outputp, self.name, self.html_content)

    def html_content(self, b, toc, subfs):
        for result in self._results:
            result.html_content(b, toc, subfs)

def sequencing_alignment(template_file, sequencing_files, outputp):
    p = SequencingAnalysis()
    p.load_fasta(template_file)
    for seqf in sequencing_files:
        print(' processing:', seqf)
        p.load_sequencingfile(seqf)
    p.write_html(outputp)
