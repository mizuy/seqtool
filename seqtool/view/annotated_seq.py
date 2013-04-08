from __future__ import absolute_import

from ..nucleotide.cpg import bisulfite_conversion
from . import svg

__all__ = ['Annotation', 'AnnotatedSeq']

def join(l, sep):
    for x in l[:-1]:
        yield x
        yield sep
    yield l[-1]

def reverse(seq):
    return str(seq)[::-1]


def negative(seq):
    return "3'-{0}-5'".format(reverse(seq))

def positive(seq):
    return "5'-{0}-3'".format(seq)

def complement(seq):
    return reverse(seq.reverse_complement())

def c_to_t(seq):
    muta = seq.tomutable()
    for i,c in enumerate(muta):
        if c=='C':
            muta[i] = 'T'
    return muta.toseq()



def iter_step(width, start, end):
    for x in range(start, end, width):
        yield x, min(x+width, end)

def step_location(x, step):
    rx = x / step
    xx = x % step

    if rx == 0:
        return xx
    elif rx == 1:
        return step + 1 + xx
    else:
        return rx * (step + 1) -1 + xx

def location(x, width, step):
    rx1 = x / width
    x1 = x % width

    return (rx1, step_location(x1,step))

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

"""
    b = svg.SvgItemsFixedHeight(8)
    b.add(svg.SvgHbar(lt, rt, svg.font_height()/2, 4, color='gray'))

    v = svg.SvgItemsVStack()

    v.add(n)
    v.add(b)
    return v
"""
class AnnotatedSeqTrack(svg.SvgMatrix):
    def __init__(self):
        self.gen = svg.SvgItemGenerator(1, 1)
        super(AnnotatedSeqTrack,self).__init__()

    def svg_seq(self, seq, p, q):
        r = svg.SvgItemsVStack()
        r.add(svg.SvgText(str(seq[p:q]), 0, 0))
        r.add(svg.SvgText(complement(seq)[p:q], 0, 0))
        return r

    def add_sequence(self, name, index, seq):
        self.add_named(name , index, " 5'-", self.gen.text(seq), "-3'")
        self.add_named(name , index, " 3'-", self.gen.text(complement(seq)), "-5'")

    def add_named(self, name, index, pre, track, post):
        self.add_row([self.gen.text(name+' '),
                     self.gen.text(str(index)+': '),
                     self.gen.text(pre), track, self.gen.text(post)],
                     ['right', 'right', 'right', None, 'left'])

    def add(self, track):
        self.add_row([None, None, None, track, None])

    def add_padding(self, height):
        self.add(svg.SvgItemsFixedHeight(height))

    def add_hline(self, length, height=5):
        t = svg.SvgItemsFixedHeight(height)
        t.add(self.gen.hline(0, length, height/2))
        self.add(t)

    def add_annotation(self, annos):
        w = svg.font_width()
        u = svg.SvgItemsVFree()
        for a in annos:
            u.add(svg_primer_bar(a.name, a.left*w, a.right*w))
        self.add(u)

class DoubleStrand(object):
    def __init__(self, seq):
        self.seq = seq
        self.positive = str(seq)
        self.negative = complement(seq)

        self.pos_anno = Annotations()
        self.neg_anno = Annotations()

class AnnotatedSeq(object):
    def __init__(self, seq):
        self.primers = []
        self.regions = []
        self.res = []

        self.length = len(seq)
        self.seq = DoubleStrand(seq)
        self.bs_pos = DoubleStrand(bisulfite_conversion(seq, True))
        self.bs_neg = DoubleStrand(bisulfite_conversion(seq.reverse_complement(), True).reverse_complement())

        self.dss = [("BS+", self.bs_pos),
                     ("origin", self.seq),
                     ("BS-", self.bs_neg)]

    def add_primer(self, primer):
        self.primers.append(primer)

        for name, ds in self.dss:
            pp,pc = primer.search(ds.seq)
            for a in pp:
                ds.pos_anno.add(primer.name, a, a+len(primer))
            for a in pc:
                ds.neg_anno.add(primer.name, a, a+len(primer))


    def add_restriction_enzymes(self, enzyme):
        self.res.append(enzyme)

    def add_region(self, name, p, q):
        self.regions.append((name,p,q))

    def track(self, width=100):

        t = AnnotatedSeqTrack()

        for p,q in iter_step(width, 0, self.length):

            for ds_name, ds in self.dss:
                t.add_padding(8)
                t.add_annotation(ds.pos_anno.cut_range(p, q))
                t.add_sequence(ds_name, p, ds.seq[p:q])
                t.add_annotation(ds.neg_anno.cut_range(p, q))
                t.add_padding(8)

            t.add_hline(width*svg.font_width())

        return t


'''
    def write_text(self, outfp, width, step, start):
        length = len(self.seq)

        for p,q in iter_step(width, 0, length):
            isize = len(str(length))
            index = str(p+start).rjust(isize)
            space = ' '.rjust(isize)

            l = loc(width)


            def loc(x):
                return step_location(x, step)

            def print_bar(r,s):
                assert p<=r<=s<=q
                loc(r), loc(s)

            def print_seq(seq):
                return ' '.join(join(seq[pp:qq] for pp, qq in iter_step(step, p, q-p), ' '))

            # print positive

            print >>outfp,index,"5'-{0}-3'".format(print_seq(self.pos_seq[p:q]))
            print >>outfp,index,"3'-{0}-5'".format(print_seq(self.neg_seq[p:q]))

            # print negative

            print >>outfp,space,"---{0}---".format("-"**l)


class AnnotatedLines(object):
    def __init__(self, seq, width=100, start=0):
        self.seq = seq
        self.aseq = AnnotatedSeq(seq)

        self.width = width

        x = 0
        while x < length:
            self.seqline.append(SeqLine(seq[x:x+width]))
            x += width

    def get_location(self, x):
        return (x / self.width, x % self.wdith)


class AnnotatedSequenceTrack(NamedTracks):
    def __init__(self, seq, width=100, step=10):
        self.gen = SvgItemGenerator(1, 1)

    def add_seqline(self):
        self.seqline.append()
        self.add_row
'''
