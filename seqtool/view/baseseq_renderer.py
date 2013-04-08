from __future__ import absolute_import

from ..nucleotide.cpg import bisulfite_conversion
from . import svg

__all__ = ['BaseseqRenderer']


def reverse(seq):
    return str(seq)[::-1]

def complement(seq):
    return reverse(seq.reverse_complement())

def iter_step(width, start, end):
    for x in range(start, end, width):
        yield x, min(x+width, end)

class Annotation(object):
    def __init__(self, name, left, right, lt_closed=True, rt_closed=True):
        self.name = name
        self.left = left
        self.right = right
        self.lt_closed = lt_closed
        self.rt_closed = rt_closed

    def overlap_range(self, left, right):
        """return True if the self overlap the range [left:right]"""
        return (self.right > left) and (right > self.left)

    def cut_range(self, left, right):
        """return intersected anntation of the self and the range [left:right]"""
        if not self.overlap_range(left, right):
            return None
        l = max(left, self.left)
        lo = not (self.left < left)
        r = min(right, self.right)
        ro = not (right < self.right)
        return Annotation(self.name, l, r, lo, ro)

    def shift(self, x):
        return Annotation(self.name, self.left+x, self.right+x, self.lt_closed, self.rt_closed)

class Annotations(object):
    def __init__(self):
        self._items = []

    def __iter__(self):
        return iter(self._items)

    def add(self, name, left, right, lt_closed=True, rt_closed=True):
        self._items.append(Annotation(name, left, right, lt_closed, rt_closed))

    def cut_range(self, left, right):
        """
        get annotations within the range [left:right]
        """
        for i in self:
            c = i.cut_range(left, right)
            if c:
                yield c.shift(-left)

def svg_primer_bar(name, lt, rt):
    n = svg.SvgText(name, lt, 0, color='black')
    r = svg.SvgRect(lt, 0, rt-lt, svg.font_height(), style='fill:none;')
    b = svg.SvgItemsFixedHeight(8)
    b.add(n)
    b.add(r)
    return b

class AnnotatedSeqTrack(svg.SvgMatrix):
    def __init__(self):
        self.gen = svg.SvgItemGenerator(1, 1)
        super(AnnotatedSeqTrack,self).__init__()

    def svg_seq(self, seq, p, q):
        r = svg.SvgItemsVStack()
        r.add(svg.SvgText(str(seq[p:q])))
        r.add(svg.SvgText(complement(seq)[p:q]))
        return r

    def add_sequence(self, name, index, seq):
        self.add_named(name , index, " 5'-", self.gen.text(seq), "-3'")
        self.add_named(name , index, " 3'-", self.gen.text(complement(seq)), "-5'")

    def add_named(self, name, index, pre, track, post):
        self.add_row([self.gen.text(name+" "),
                     self.gen.text(" {0}: ".format(str(index))),
                     self.gen.text(pre), track, self.gen.text(post)],
                     ['right', 'right', 'right', None, 'left'])

    def add(self, track):
        self.add_row([None, None, None, track, None])

    def add_padding(self, height):
        self.add(svg.SvgItemsFixedHeight(height))

    def add_hline(self, length, height=5, **kwargs):
        t = svg.SvgItemsFixedHeight(height)
        t.add(self.gen.hline(0, length, height/2, **kwargs))
        self.add(t)

    def add_annotation(self, annos):
        w = svg.font_width()
        u = svg.SvgItemsVFree()

        flag = False
        for a in annos:
            u.add(svg_primer_bar(a.name, a.left*w, a.right*w))
            flag = True
        if flag:
            self.add(u)

class DoubleStrand(object):
    def __init__(self, seq):
        self.seq = seq
        self.positive = str(seq)
        self.negative = complement(seq)

        self.pos_anno = Annotations()
        self.neg_anno = Annotations()

class BaseseqRenderer(object):
    def __init__(self, seq):
        self.primers = []
        self.regions = []

        self.length = len(seq)
        self.seq = DoubleStrand(seq)
        self.bs_pos = DoubleStrand(bisulfite_conversion(seq, True))
        self.bs_neg = DoubleStrand(bisulfite_conversion(seq.reverse_complement(), True).reverse_complement())

        self.dss = [("BS(+)", self.bs_pos),
                     ("", self.seq),
                     ("BS(-)", self.bs_neg)]

    def add_primer(self, primer):
        self.primers.append(primer)

        for name, ds in self.dss:
            pp,pc = primer.search(ds.seq)
            for a in pp:
                ds.pos_anno.add(primer.name, a.left, a.right)
            for a in pc:
                ds.neg_anno.add(primer.name, a.left, a.right)


    def add_restriction_batche(self, restriction_batch):
        for name, ds in self.dss:
            for enzyme, locs in restriction_batch.search(ds.seq).items():
                if len(locs)>1:
                    continue
                for l in locs:
                    ds.pos_anno.add(str(enzyme), l, l+len(enzyme.site))

    def add_region(self, name, p, q):
        self.regions.append((name,p,q))

    def track(self, width=200):

        t = AnnotatedSeqTrack()

        t.add_hline(width*svg.font_width(), color='black', **{'stroke-width':'3'})

        for p,q in iter_step(width, 0, self.length):

            count = 0
            for ds_name, ds in self.dss:
                if count != 0:
                    t.add_hline(width*svg.font_width(), color='gray')
                count += 1

                t.add_padding(8)
                t.add_annotation(ds.pos_anno.cut_range(p, q))
                t.add_sequence(ds_name, p, ds.seq[p:q])
                t.add_annotation(ds.neg_anno.cut_range(p, q))
                t.add_padding(8)

            t.add_hline(width*svg.font_width(), color='black', **{'stroke-width':'3'})

        return t
