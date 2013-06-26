
from ..nucleotide.cpg import bisulfite_conversion
from . import svg
from .rectangle import Rectangle,Line
from ..nucleotide import reverse, complement

__all__ = ['BaseseqRenderer']

class Annotation:
    def __init__(self, name, left, right):
        self.name = name
        self.range = Line(left,right)
        w = svg.font_width()
        self._vr = self.range.scale(w)

    def svg_item(self, r=None):
        if r and not self.range.has_intersect(r):
            return None

    @classmethod
    def get_track(self, annotations, reverse, p, q):
        u = svg.SvgItemsVFree(reverse)
        for a in annotations:
            item = a.svg_item(Line(p,q))
            if item:
                u.add(item)
        w = svg.font_width()
        return svg.SvgClipWidth(0, w*(q-p), svg.SvgTranslate(-w*p, 0, u))

class AnnotationTexts(Annotation):
    def __init__(self, name, left, right, texts=[]):
        super().__init__(name, left, right)
        self.texts = texts

    def svg_item(self, r=None):
        if r and not self.range.has_intersect(r):
            return None

        t = svg.SvgItemsVStack()
        t.add(svg.SvgText(self.name, self._vr.start, 0))
        for s in self.texts:
            t.add(svg.SvgText(s, self._vr.start, 0))

        return svg.SvgBoundbox(svg.SvgExpandWidth(self._vr.length, t))

class AnnotationPrimerAnneal(Annotation):
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
                # todo complementary for - strand...
                al = pta.right
                ar = pta.right + pta.adapter_length
                self.left = pta.left
            self.adapter = (pta.primer_adapter, al, ar)

        self.right = pta.right

        super().__init__(pta.primer.name, self.left, self.right)

    def svg_item(self, r=None):
        if r and not self.range.has_intersect(r):
            return None
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


class NamedStack(svg.SvgMatrix):
    def __init__(self):
        super().__init__()

    def add_seq(self, name, pre, seq, post):
        self.add_row([svg.SvgText(name), svg.SvgText(pre), svg.SvgText(seq), svg.SvgText(post)])

    def add_named(self, name, t):
        self.add_row([svg.SvgText(name), svg.SvgBase(), t, svg.SvgBase()])

    def add(self, t):
        self.add_row([svg.SvgBase(), svg.SvgBase(), t, svg.SvgBase()])


class DoubleStrand:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

        self.pos = []
        self.neg = []

    def get_track(self, p, q):
        t = NamedStack()

        # padding
        t.add(svg.SvgItemsFixedHeight(8))

        # positive strand annotations
        t.add(Annotation.get_track(self.pos, True, p, q))

        # sequences
        s = self.seq[p:q]
        t.add_seq(self.name, "5'-", s, "-3'")
        t.add_seq(self.name, "3'-", complement(s), "-5'")

        # negative strand annotations
        t.add(Annotation.get_track(self.neg, False, p, q))

        # padding
        t.add(svg.SvgItemsFixedHeight(8))

        return t

    def add_primer(self, primer):
        pp,pc = primer.search(self.seq)
        for a in pp:
            self.pos.append(AnnotationPrimerAnneal(a))
        for a in pc:
            self.neg.append(AnnotationPrimerAnneal(a))

    def add_restriction_batch(self, rb):
        for enzyme, locs in list(rb.search(self.seq).items()):
            #if len(locs)>1:
            #    continue
            for l in locs:
                # TODO: pretty prent restriction cut pattern
                p = l -1 - enzyme.fst5
                self.pos.append(AnnotationTexts(str(enzyme), p, p+len(enzyme.site),[enzyme.site]))



NO_CONVERSIONS = [("", lambda x: x)]


BISULFITE_CONVERSIONS = [
    ("BS(+)", lambda x: bisulfite_conversion(x, sense=True)),
    ("", lambda x: x),
    ("BS(-)", lambda x: bisulfite_conversion(x, sense=False)),
]


def iter_step(width, start, end):
    for x in range(start, end, width):
        yield x, min(x+width, end)


class BaseseqRenderer:
    def __init__(self, seq, bisulfite=False):
        self.length = len(seq)

        conversions = NO_CONVERSIONS
        if bisulfite:
            conversions = BISULFITE_CONVERSIONS

        self.doublestrands = []
        for name, conv in conversions:
            self.doublestrands.append(DoubleStrand(name, conv(seq)))

    def add_primer(self, primer):
        for ds in self.doublestrands:
            ds.add_primer(primer)

    def add_restriction_batch(self, restriction_batch):
        for ds in self.doublestrands:
            ds.add_restriction_batch(restriction_batch)

    def add_alignment(self, name, p, q, bars):
        for ds in self.doublestrands:
            ds.pos.append(AnnotationTexts(ds.name, p, q, bars))

    def track(self, width):
        l = self.length
        w = svg.font_width()

        t = svg.SvgItemsVStack()

        for p,q in iter_step(width, 0, l):
            # TODO aligned Hlinebox
            t.add(svg.SvgHlineBox(0, width*w, 10))

            count = 0
            for ds in self.doublestrands:
                if count != 0:
                    t.add(svg.SvgHlineBox(0, width*w, 0.1))
                count += 1

                t.add(ds.get_track(p,q))

        t.add(svg.SvgHlineBox(0, width*w, 10))

        return svg.SvgPadding(10,10,t)



