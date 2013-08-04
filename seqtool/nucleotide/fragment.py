
STRAND_POS = True
STRAND_NEG = False

class FragmentEnd(object):
    def __init__(self, strand=None, seq=''):
        """
        strand is True or False
        seq is always sense strand sequence
        """
        self.strand = strand
        self.seq = seq

class Fragment(object):
    """
    Double strand linear DNA Fragment with blunt or sticky end.
    """
    def __init__(self, body, p5=FragmentEnd(), p3=FragmentEnd()):
        self.p5 = p5
        self.body = body
        self.p3 = p3

    def __repr__(self):
        ret = self.p5.seq \
            + ('_' if self.p5.strand else '^') \
            + self.body \
            + ('_' if self.p3.strand else '^') \
            + self.p3.seq
        return "Fragment({})".format(ret)

    def digestion(self, ri):
        return []

    def ligation(self, flag):
        return []

