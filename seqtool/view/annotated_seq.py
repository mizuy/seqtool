


def reverse(seq):
    return str(seq)[::-1]

def complement(seq):
    return seq.complement()
    #reverse(seq.reverse_complement())

def iter_step(width, start, end):
    for x in range(start, end, width):
        yield x, min(x+width, end)


class Annotation(object):
    def __init__(self, name, left, right, sequences=[], lt_closed=True, rt_closed=True):
        self.name = name
        self.left = left
        self.right = right
        self.sequences = sequences
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
        seq = [s[l-self.left:r-self.left] for s in self.sequences]
        return Annotation(self.name, l, r, seq, lo, ro)

    def shift(self, x):
        return Annotation(self.name, self.left+x, self.right+x, self.sequences, self.lt_closed, self.rt_closed)

class PrimerAnnealAnno(Annotation):
    def __init__(self):
        pass

class Annotations(object):
    def __init__(self):
        self._items = []

    def __iter__(self):
        yield from self._items

    def add(self, anno):
        self._items.append(anno)

    def cut_range(self, left, right):
        """
        get annotations within the range [left:right]
        """
        for i in self:
            c = i.cut_range(left, right)
            if c:
                yield c.shift(-left)

class AnnotatedSeq:
    def __init__(self, seq):
        self.seq = seq
        self.pos_anno = Annotations()
        self.neg_anno = Annotations()

    def add_primer(self, primer):
        self.primers.append(primer)

        pp,pc = primer.search(self.seq)
        for a in pp:
            name = primer.name + ('' if a.full else '(partial)')
            self.pos_anno.add(Annotation(name, a.left, a.right))
        for a in pc:
            name = primer.name + ('' if a.full else '(partial)')
            self.neg_anno.add(Annotation(name, a.left, a.right))

    def add_restriction_batch(self, restriction_batch):
        for enzyme, locs in list(restriction_batch.search(ds.seq).items()):
            #if len(locs)>1:
            #    continue
            for l in locs:
                # TODO: pretty prent restriction cut pattern
                p = l -1 - enzyme.fst5
                ds.pos_anno.add(Annotation(str(enzyme), p, p+len(enzyme.site)))

    def add_alignment(self, name, p, q, bars):
        self.seq.pos_anno.add(Annotation(name, p, q, bars))

    def track(self, width=160):
        width = width or self.length

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


class AnnotatedSeqConv:
    def __init__(self, seq, conversions):
        self.aseqs = []
        for name, conv in conversions:
            self.aseqs.append(name, AnnotaedSeq(conv(seq)))

bisulfite_conversions = [
    ("BS(+)", lambda x: bisulfite_conversion(x, sense=True)),
    ("", lambda x: x),
    ("BS(-)", lambda x: bisulfite_conversion(x, sense=False)),
]

class BisulfiteAnnotatedSeqConv:
    def __init__(self, seq):
        super().__init__(seq, bisulfite_conversions)

