
from ..nucleotide.cpg import bisulfite_conversion
from ..util import svg
from ..util.rectangle import Rectangle,Line
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
        
        self.matched = (pta.display_match, pta.left, pta.right)
        self.left = pta.leftmost
        self.right = pta.right

        if not pta.full:
            self.adapter = (pta.display_adapter, pta.a_l, pta.a_r)


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

    def add_subline(self, width):
        self.add(svg.SvgHlineBox(0, width, 10, **{'stroke-width':0.1, 'stroke-dasharray':'5 5'}))

    def add_line(self, width):
        self.add(svg.SvgHlineBox(0, width, 10))

class DoubleStrand:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

        self.pos = []
        self.neg = []

    def assign_track(self, named_track, p, q):
        t = named_track

        # padding
        t.add(svg.SvgItemsFixedHeight(8))

        # positive strand annotations
        t.add(Annotation.get_track(self.pos, True, p, q))

        # sequences
        s = self.seq[p:q]
        t.add_seq(self.name, "5'-", s, "-3'")
        t.add_seq('', "3'-", complement(s), "-5'")

        # negative strand annotations
        t.add(Annotation.get_track(self.neg, False, p, q))

        # padding
        t.add(svg.SvgItemsFixedHeight(8))

    def add_primer(self, primer):
        pp,pc = primer.search(self.seq, template_ambiguous = True)
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

def iter_sep(iterable, sep, firstend = True):
    count =  0
    for i in iterable:
        if count == 0:
            if firstend:
                sep()
        else:
            sep()
        yield i
        count += 1

    if firstend:
        sep()
        

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
            ds.pos.append(AnnotationTexts(name, p, q, bars))

    def track_partial(self, start, end, width = None):
        start = max(start, 0)
        end = min(end, self.length)

        width = width or (end - start)
        w = svg.font_width()
        ww =  width * w

        t = NamedStack()

        for p,q in iter_sep(iter_step(width, start, end), lambda:t.add_line(ww)):
            for ds in iter_sep(self.doublestrands, lambda:t.add_subline(ww), False):
                ds.assign_track(t, p, q)

        return svg.SvgPadding(20,20,t)


    def track(self, width = None):
        return self.track_partial(0, self.length, width)
