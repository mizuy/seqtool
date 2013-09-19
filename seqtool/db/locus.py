__all__ = ['Locus','StrandPos']

class StrandPos(object):
    def __init__(self, sense, lower, higher):
        """
        sense: True or False
        always, lower < higher
        for sense,     lower,higher = start,end
        for antisense, lower,higher = end,start

        chr     --------------------------------------------------
                     lower                               higher

        sense 5' -------------------------------------------------> 3'
                     start                               stop

        anti  3' <------------------------------------------------- 5'
                     stop                                start
        """
        assert(lower <= higher)
        self.sense = sense
        self.lower = lower
        self.higher = higher

    def __repr__(self):
        return 'StrandPos({}, {}, {})'.format(self.sense, self.lower, self.higher)

    def expand(self, upstream, downstream):
        """
        >>> StrandPos(True, 2000, 4000).expand(100, 200)
        StrandPos(True, 1900, 4200)
        >>> StrandPos(False, 2000, 4000).expand(100, 200)
        StrandPos(False, 1800, 4100)
        """
        if self.sense:
            # lower - higher
            lower = self.lower - upstream
            higher = self.higher + downstream
        else:
            # higher - lower
            lower = self.lower - downstream
            higher = self.higher + upstream

        return StrandPos(self.sense, lower, higher)

    def antisense(self):
        """
        >>> StrandPos(True, 2000, 4000).antisense()
        StrandPos(False, 2000, 4000)
        """
        return StrandPos(not self.sense, self.lower, self.higher)

    @property
    def start(self):
        """
        return location of 5' end

        >>> StrandPos(True, 2000, 4000).start
        2000
        >>> StrandPos(False, 2000, 4000).start
        4000
        """
        if self.sense:
            return self.lower
        else:
            return self.higher
    
    @property
    def stop(self):
        """
        return location of 3' end

        >>> StrandPos(True, 2000, 4000).stop
        4000
        >>> StrandPos(False, 2000, 4000).stop
        2000
        """
        if self.sense:
            return self.higher
        else:
            return self.lower
        

    def index_5(self, i):
        """
        return location of 5' end + i

        >>> StrandPos(True, 2000, 4000).index_5(100)
        2100
        >>> StrandPos(False, 2000, 4000).index_5(100)
        3900
        """
        if self.sense:
            return self.lower+i
        else:
            return self.higher-i

    def index_3(self, i):
        """
        return location of 5' end + i

        >>> StrandPos(True, 2000, 4000).index_3(100)
        4100
        >>> StrandPos(False, 2000, 4000).index_3(100)
        1900
        """
        if self.sense:
            return self.higher+i
        else:
            return self.lower-i
        
    def rel_5(self, i):
        """
        >>> StrandPos(True, 2000, 4000).rel_5(100)
        -1900
        >>> StrandPos(False, 2000, 4000).rel_5(100)
        3900
        """
        if self.sense:
            return i - self.lower
        else:
            return self.higher - i

class Locus(object):
    @classmethod
    def from_strandpos(self, chromosome, strandpos):
        return Locus(chromosome, strandpos.sense, strandpos.lower, strandpos.higher)
        
    def __init__(self, chromosome, sense, lower, higher):
        self.chrom = chromosome
        self.pos = StrandPos(sense, lower, higher)

    def __repr__(self):
        return "Locus('{}', {})".format(self.chrom, self.pos)

    def expand(self, upstream, downstream):
        """
        >>> Locus('chr1', True, 2000, 4000).expand(1000,2000)
        Locus('chr1', StrandPos(True, 1000, 6000))
        >>> Locus('chr1', False, 2000, 4000).expand(1000,2000)
        Locus('chr1', StrandPos(False, 0, 5000))
        """
        return Locus.from_strandpos(self.chrom, self.pos.expand(upstream, downstream))

    def antisense(self):
        return Locus.from_strandpos(self.chrom, self.pos.antisense)
