import os
from Bio import SeqIO
import io
from ..nucleotide import to_seq
from ..nucleotide.alignment import make_alignment
from ..nucleotide.cpg import bisulfite_conversion, c2t_conversion
from ..nucleotide.primer import Primer
from ..util import report
from ..view.baseseq_renderer import BaseseqRenderer
from ..format.abi import AbiFormat
from ..format.render import SvgPeaksAlignment
from ..util.debug import report_exceptions
from collections import defaultdict
from ..util import DefaultOrderedDict
from .bsa_figure import BsaFigure
# TODO: .scf file support (4peaks output)
# cf.
# http://staden.sourceforge.net/manual/formats_unix_3.html#SEC3
# format
# http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-pysrc.html
# abi = SeqIO.read(filename, "abi")
        
class NamedLists(DefaultOrderedDict):
    def __init__(self):
        super().__init__(list)

    def map(self, func):
        ret = NamedLists()
        for key, values in self.items():
            ret[key] = [func(v) for v in values]
        return ret

def values_flatten(d):
    for key, values in d.items():
        yield from values

        
def get_(a0, a1):
    return a0 if a0.score_density() > a1.score_density() else a1
    
class PairSeqAlignment:
    def __init__(self, seq0, seq1):
        self.seq0 = seq0
        self.seq1 = seq1

        self.s00 = make_alignment(seq0.seq, seq1.seq)
        self.s01 = make_alignment(seq0.seq, seq1.rcs)
        self.s10 = self.s01.reverse_complement() #make_alignment(seq0.rcs, seq1.seq)
        self.s11 = self.s00.reverse_complement() #make_alignment(seq0.rcs, seq1.rcs)
        
    def get_alignments_seq0(self, sense):
        """
        return best alignments for seq0[sense]
        """
        if sense:
            return get_(self.s00, self.s01)
        else:
            return get_(self.s10, self.s11)

    def get_alignments_seq1(self, sense):
        if sense:
            return get_(self.s00, self.s10)
        else:
            return get_(self.s01, self.s11)


class NamedSequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        self.rcs = str(to_seq(self.seq).reverse_complement())

    def reversed(self):
        return Target(self.name+'(reversed)', self.rcs)

class Template(NamedSequence):
    pass

class Target(NamedSequence):
    @classmethod
    def from_filename(self, filename):
        base,ext = os.path.splitext(filename)
        if ext=='.seq':
            abi = None
            seq = to_seq(''.join(i.strip() for i in open(filename,'rU').readlines()))
        elif ext=='.ab1':
            abi = AbiFormat(filename)
            seq = to_seq(abi.view.get_sequence())
        return Target(os.path.basename(filename), seq, abi, abi.get_view(True))
        
    def __init__(self, name, seq, abi, view):
        super().__init__(name, seq)
        self.abi = abi
        self.view = view

    def reversed(self):
        return Target(self.name+'(reversed)', self.rcs, self.abi, self.abi.get_view(False))

class Matrix:
    def __init__(self, targets_dict, templates_dict):
        self.targets = targets_dict
        self.templates = templates_dict
        
        self.alignments = {}
        for k,ss in self.targets.items():
            for s in ss:
                print('calculating: {}'.format(s.name))
                for l,tt in self.templates.items():
                    for t in tt:
                        self.alignments[(s,t)] = PairSeqAlignment(s, t)

    def each_target(self):
        return self.targets.map( \
                    lambda s: ViewTarget(s, \
                                         self.templates.map(lambda t: self.alignments[(s,t)] )))
    def each_template(self):
        return self.templates.map( \
                    lambda t: ViewTemplate(t, \
                                           self.targets.map(lambda s: self.alignments[(s,t)] )))

    def html_content_each_target(self, b, toc, subfs):
        for keyname, views in self.each_target().items():
            with report.section(b, toc, keyname):
                for view in views:
                    view.html_content(b, toc, subfs)


    def html_content_each_template(self, b, toc, subfs):
        for keyname, views in self.each_template().items():
            with report.section(b, toc, keyname):
                for view in views:
                    view.html_content(b, toc, subfs)
        

    def html_content(self, b, toc, subfs):
        l = list(values_flatten(self.each_template()))
        if l:
            first = l[0]
            content = first.bs_result_file()
            bf = BsaFigure()
            bf.readfp(io.StringIO(content), 'bs_result.txt')
            with b.a(href=subfs.write('bs_result.txt', content)):
                b.text("bisulfite sequencing analysis result text")
            with b.a(href=subfs.write('bs_result.svg', bf.svg(False))):
                b.text("bisulfite sequencing analysis result svg")
                
        self.html_content_each_target(b, toc, subfs)

                    
class ViewTarget:
    def __init__(self, target, templates_dict):
        self.target = target
        self.templates = templates_dict
                
    def html_content(self, b, toc, subfs):
        with report.section(b, toc, self.target.name):
            for sense in [True,False]:
                with report.section(b, toc, 'Sense' if sense else 'Reverse Complement'):

                    target = self.target if sense else self.target.reversed()
                    seq = target.seq
                    name = target.name+('_sense_' if sense else '_antisense_')

                    fs = []
                    if self.target.abi:
                        fs.append((name+'_peak.svg', self.svg_peak(self.target.view, sense)))

                    if fs:
                        with b.ul:
                            for filename, content in fs:
                                with b.li:
                                    with b.a(href=subfs.write(filename, content)):
                                        b.text(filename)

                    for keyname, alignments in self.templates.items():
                        with report.section(b, toc, keyname):
                            for pal in alignments:
                                al = pal.get_alignments_seq0(sense)
                                if not al.criteria():
                                    continue

                                with report.section(b, toc, pal.seq1.name):
                                    cm = al.reversed().correspondance_map()

                                    with b.pre:
                                        #b.text(cm.text_map())
                                        b.text(cm.text_str())
                                        b.text(cm.get_bsa_result())

                                    with b.pre:
                                        b.text(al.text_shrinked())


    def svg_peak(self, view, sense):
        render = SvgPeaksAlignment()
        seq_loc = view.get_location()

        for keyname, alignments in self.templates.items():
            for pal in alignments:
                al = pal.get_alignments_seq0(sense)
                tempname = pal.seq1.name

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

    def __repr__(self):
        return 'ViewTarget({})'.format(self.target.name)
        
class ViewTemplate:
    def __init__(self, template, targets):
        self.template = template
        self.targets = targets

    def __repr__(self):
        return 'ViewTemplate({})'.format(self.template.name)
                
    def html_content(self, b, toc, subfs):
        with report.section(b, toc, self.template.name):
            for keyname, alignments in self.targets.items():
                with report.section(b, toc, keyname):

                    for pal in alignments:
                        al = pal.get_alignments_seq1(True)
                        if not al.criteria():
                            continue

                        with report.section(b, toc, pal.seq0.name):
                            cm = al.reversed().correspondance_map()

                            with b.pre:
                                #b.text(cm.text_map())
                                b.text(cm.text_str())

                            with b.pre:
                                b.text(al.text_shrinked())

    def bs_result_file(self):
        to = io.StringIO()
        print('settings:', file=to)
        print('    template:', self.template.seq, file=to)
        print('', file=to)
        print('results:', file=to)
        for keyname, alignments in self.targets.items():

            print('    {}:'.format(keyname), file=to)
            
            for pal in alignments:
                al = pal.get_alignments_seq1(True)
                r = al.reversed().get_match_result()
                print('        {}: {}'.format(pal.seq0.name, r), file=to)
            print('\n', file=to)
        return to.getvalue()


def load_templates_fasta(filename):
    ret = NamedLists()

    for record in SeqIO.parse(open(filename,'r'), "fasta"):
        name, sep, conv = record.description.partition(':')
        conv = conv.strip()
        seq = record.seq

        #print name, conv

        if conv=='C2T':
            ret[name] = [Template(name, seq),
                         Template(name+' C2T+', c2t_conversion(seq, True)),
                         Template(name+' C2T-', c2t_conversion(seq, False))]
        elif conv=='BS':
            ret[name] = [Template(name, seq),
                         Template(name+' BS+', bisulfite_conversion(seq, True)),
                         Template(name+' BS-', bisulfite_conversion(seq, False))]
        elif conv=='Primer':
            ret['primer'].append(Template(name,seq))
        else:
            ret['others'].append(Template(name,seq))
    return ret

def sequencing_alignment(template_file, sequencing_files, outputp):
    return sequencing_alignment_dir(template_file, {'top': sequencing_files}, outputp)

def sequencing_alignment_dir(template_file, sequencing_files_dict, outputp):
    print('loading: {}'.format(template_file))
    templates = load_templates_fasta(template_file)

    targets = NamedLists()
    for keyname, seqfs in sequencing_files_dict.items():
        for seqf in seqfs:
            print('loading: {}'.format(seqf))
            targets[keyname].append(Target.from_filename(seqf))
    
    m = Matrix(targets, templates)

    report.write_html(outputp, 'Sequencing Analysis', m.html_content)

