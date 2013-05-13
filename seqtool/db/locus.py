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

    def expand(self, upstream, downstream):
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
        print 'ANTISENSE ???'
        return StrandPos(not self.sense, self.lower, self.higher)

    @property
    def start(self):
        """
        return location of 5' end
        """
        if self.sense:
            return self.lower
        else:
            return self.higher
    
    @property
    def stop(self):
        """
        return location of 3' end
        """
        if self.sense:
            return self.higher
        else:
            return self.lower
        

    def index_5(self, i):
        """
        return location of 5' end + i
        """
        if self.sense:
            return self.lower+i
        else:
            return self.higher-i

    def index_3(self, i):
        """
        return location of 5' end + i
        """
        if self.sense:
            return self.higher+i
        else:
            return self.lower-i
        
    def rel_5(self, i):
        if self.sense:
            return i - self.lower
        else:
            return self.higher - i

class Locus(object):
    def __init__(self, chromosome, sense, lower, higher):
        self.chrom = chromosome
        self.pos = StrandPos(sense, lower, higher)

    def expand(self, upstream, downstream):
        p = self.pos.expand(upstream, downstream)
        return Locus(self.chrom, p.sense, p.lower, p.higher)

    def antisense(self):
        p = self.pos.antisense()
        return Locus(self.chrom, p.sense, p.lower, p.higher)
