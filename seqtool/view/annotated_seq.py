

def reverse(seq):
    return str(seq)[::-1]

def complement(seq):
    return seq.complement()
    #reverse(seq.reverse_complement())

class SeqWrappingTrack(svg.SvgMatrix):
    def __init__(self, length, width):
        self.gen = svg.SvgItemGenerator(1, 1)
        super().__init__()

        self.length = length
        self.width = width
        self.tracks = []
        for p,q in iter_step(self.width, 0, self.length):
            stack = svg.SvgItemsVStack()
            group = svg.SvgTranslate(-p, 0, stack) # clip-path

            self.tracks.append(stack)


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
            lt = a.left*w
            rt = a.right*w
            h = svg.font_height()

            if not a.sequences:
                b = svg.SvgItemsFixedHeight(h)
                b.add(svg.SvgText(a.name, lt, 0))
                b.add(svg.SvgRect(lt, 0, rt-lt, h, style='fill:none;'))
            else:
                t = svg.SvgItemsVStack()
                t.add(svg.SvgText(a.name, lt, 0))
                for s in a.sequences:
                    t.add(svg.SvgText(s, lt, 0))

                b = svg.SvgItemsFixedHeight(t.rect.height)
                b.add(t)
                b.add(svg.SvgRect(lt, 0, rt-lt, t.rect.height, style='fill:none;'))
            u.add(b)

            flag = True
        if flag:
            self.add(u)

class AnnotatedSeqTrack(WrappingTrack):
    def add_sequence(self, name, index, seq):
        self.add_named(name , index, " 5'-", self.gen.text(seq), "-3'")
        self.add_named(name , index, " 3'-", self.gen.text(complement(seq)), "-5'")

def iter_step(width, start, end):
    for x in range(start, end, width):
        yield x, min(x+width, end)


class Annotation:
    def __init__(self, name, left, right, sequences=[]):
        self.name = name
        self.left = left
        self.right = right
        self.sequences = sequences

    def draw(self, stack, l, r):
        pass

class PrimerAnno(Annotation):
    pass

class RestrictionEnzymeAnno(Annotation):
    pass

class Annotations:
    def __init__(self):
        self._items = []

    def __iter__(self):
        yield from self._items

    def add(self, anno):
        self._items.append(anno)


class AnnotatedSeq:
    def __init__(self, seq):
        self.seq = seq
        self.anno = Anntations()
        self.pos_anno = Annotations()
        self.neg_anno = Annotations()

    def add_primer(self, primer):
        pp,pc = primer.search(self.seq)
        for a in pp:
            name = primer.name + ('' if a.full else '(partial)')
            self.pos_anno.add(PrimerAnno(name, a.left, a.right))
        for a in pc:
            name = primer.name + ('' if a.full else '(partial)')
            self.neg_anno.add(PrimerAnno(name, a.left, a.right))

    def add_restriction_batch(self, restriction_batch):
        for enzyme, locs in restriction_batch.search(ds.seq).items():
            #if len(locs)>1:
            #    continue
            for l in locs:
                # TODO: pretty prent restriction cut pattern
                p = l -1 - enzyme.fst5
                ds.anno.add(RestrictionEnzymeAnno(str(enzyme), p, p+len(enzyme.site)))

    def add_alignment(self, name, p, q, bars):
        self.seq.pos_anno.add(Annotation(name, p, q, bars))

    def track(self, width=160):
        width = width or self.length

        t = AnnotatedSeqTrack()

        for p,q in iter_step(width, 0, self.length):
            t = t.get_track(p,q)
            
            for a in self.pos_anno:
                a.draw(t, p, q)
            
            for a in self.anno:
                a.draw(t, p, q)
            
            t.add_sequence(ds_name, p, ds.seq[p:q])

            for a in self.neg_anno:
                a.draw(t, p, q)

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

