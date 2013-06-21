from ..nucleotide.cpg import bisulfite_conversion
from . import svg

__all__ = ['BaseseqRenderer']

def reverse(seq):
    return str(seq)[::-1]

def complement(seq):
    return seq.complement()


class Annotation:
    def __init__(self, name, left, right):
        self.name = name
        self.left = left
        self.right = right

    def svg_item_filtered(self, p, q):
        if p < self.left or 

    def svg_item(self):
        w = svg.font_width()
        lt = self.left*w
        rt = self.right*w
        h = svg.font_height()

        b = svg.SvgItemsFixedHeight(h)
        b.add(svg.SvgText(self.name, lt, 0))
        b.add(svg.SvgRect(lt, 0, rt-lt, h, style='fill:none;', stroke='black'))
        return b

class AnnotationTexts(Annotation):
    def __init__(self, name, left, right, texts=[]):
        super().__init__(name, left, right)
        self.texts = texts

    def svg_item(self):
        w = svg.font_width()
        lt = self.left*w
        rt = self.right*w

        t = svg.SvgItemsVStack()
        t.add(svg.SvgText(self.name, lt, 0))
        for s in self.texts:
            t.add(svg.SvgText(s, lt, 0))

        b = svg.SvgItemsFixedHeight(t.rect.height)
        b.add(t)
        b.add(svg.SvgRect(lt, 0, rt-lt, t.rect.height, style='fill:none;', stroke='black'))
        return b

class AnnotationPrimerAnneal:
    def __init__(self, pta):
        self.pta = pta
        self.name = pta.primer.name
        
        self.matched = (pta.primer_match, pta.left, pta.right)
        self.left = pta.left

        if not pta.full:
            if pta.strand:
                al = pta.left - pta.adapter_length
                ar = pta.left
                self.left = al
            else:
                al = pta.right
                ar = pta.right + pta.adapter_length
                self.left = pta.left
            self.adapter = (pta.primer_adapter, al, ar)

    def svg_item(self):
        w = svg.font_width()
        h = svg.font_height()

        ps = svg.SvgItemsFixedHeight(h)
        s,l,r = self.matched
        ps.add(svg.SvgText(s, l*w, 0))

        if not self.pta.full:
            s,l,r = self.adapter
            ps.add(svg.SvgBoundbox(svg.SvgText(s, l*w, 0)))

        t = svg.SvgItemsVStack()
        t.add(svg.SvgText(self.name, self.left*w, 0))
        t.add(ps)

        return svg.SvgBoundbox(t)


class AnnotatedSeqTrack(svg.SvgMatrix):
    def __init__(self):
        self.gen = svg.SvgItemGenerator(1, 1)
        super(AnnotatedSeqTrack,self).__init__()

    def add_sequence(self, name, index, seq):
        self.add_named(name , index, " 5'-", self.gen.text(seq), "-3'")
        self.add_named('' , '', " 3'-", self.gen.text(complement(seq)), "-5'")

    def add_named(self, name, index, pre, track, post):
        self.add_row([self.gen.text(name+" "),
                     self.gen.text(" {0}: ".format(str(index))),
                     self.gen.text(pre), track, self.gen.text(post)],
                     ['right', 'right', 'right', None, 'left'])

    def add(self, track):
        self.add_row([None, None, None, track, None])

    def add_padding(self, height):
        self.add(svg.SvgItemsFixedHeight(height))

    def add_hline(self, length, w, height=5):
        t = svg.SvgItemsFixedHeight(height)
        t.add(self.gen.hline(0, length, height/2, **{'stroke-width':str(w)}))
        self.add(t)

    def add_annotation(self, annos, reverse):
        if not annos:
            return

        u = svg.SvgItemsVFree(reverse)
        for a in annos:
            u.add(a.svg_item())
        self.add(u)

    def add_doublestrand(self, ds):
        self.add_padding(8)
        self.add_annotation(ds.pos, True)
        self.add_sequence(ds.name, 0, ds.seq)
        self.add_annotation(ds.neg, False)
        self.add_padding(8)


class DoubleStrand:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

        self.pos = []
        self.neg = []

no_conversions = [("", lambda x: x)]

bisulfite_conversions = [
    ("BS(+)", lambda x: bisulfite_conversion(x, sense=True)),
    ("", lambda x: x),
    ("BS(-)", lambda x: bisulfite_conversion(x, sense=False)),
]


class BaseseqRenderer:
    def __init__(self, seq, bisulfite=False):
        self.length = len(seq)

        conversions = no_conversions
        if bisulfite:
            conversions = bisulfite_conversions

        self.dss = []
        for name, conv in conversions:
            self.dss.append(DoubleStrand(name, conv(seq)))

    def add_primer(self, primer):
        for ds in self.dss:
            pp,pc = primer.search(ds.seq)
            for a in pp:
                ds.pos.append(AnnotationPrimerAnneal(a))
            for a in pc:
                ds.neg.append(AnnotationPrimerAnneal(a))


    def add_restriction_batch(self, restriction_batch):
        for ds in self.dss:
            for enzyme, locs in list(restriction_batch.search(ds.seq).items()):
                #if len(locs)>1:
                #    continue
                for l in locs:
                    # TODO: pretty prent restriction cut pattern
                    p = l -1 - enzyme.fst5
                    ds.pos.append(AnnotationTexts(str(enzyme), p, p+len(enzyme.site),[enzyme.site]))

    def add_alignment(self, name, p, q, bars):
        for ds in self.dss:
            ds.pos.append(AnnotationTexts(ds.name, p, q, bars))

    def track(self, width):
        l = self.length
        w = svg.font_width()

        t = AnnotatedSeqTrack()

        t.add_hline(l*w, 3)

        count = 0
        for ds in self.dss:
            if count != 0:
                t.add_hline(l*w, 0.5)
            count += 1

            t.add_doublestrand(ds)

        t.add_hline(l*w, 3)

        p = svg.SvgItemsWrapping(width*w, t)

        return p



