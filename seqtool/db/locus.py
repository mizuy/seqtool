__all__ = ['Locus']

class Locus(object):
    def __init__(self, chromosome, strand, lower, higher):
        """
        chromosome: chr?
        strand: True or False
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
        self.chromosome = chromosome
        self.strand = strand
        self.lower = lower
        self.higher = higher

    def expand(self, upstream, downstream):
        if self.strand:
            # lower - higher
            lower = self.lower - upstream
            higher = self.higher + downstream
        else:
            # higher - lower
            lower = self.lower - downstream
            higher = self.higher + upstream

        return Locus(self.chromosome, self.strand, lower, higher)

    def antisense(self):
        return Locus(self.chromosome, not self.strand, lower, higher)

    @property
    def start(self):
        """
        return location of 5' end
        """
        if self.strand:
            return self.lower
        else:
            return self.higher
    
    @property
    def stop(self):
        """
        return location of 3' end
        """
        if self.strand:
            return self.higher
        else:
            return self.lower
        

    def index_5(self, i):
        """
        return location of 5' end + i
        """
        if self.strand:
            return self.lower+i
        else:
            return self.higher-i

    def index_3(self, i):
        """
        return location of 5' end + i
        """
        if self.strand:
            return self.higher+i
        else:
            return self.lower-i
        
    def rel_5(self, i):
        if self.strand:
            return i - self.lower
        else:
            return self.higher - i
