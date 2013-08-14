
__all__=['PrimerCondition','NORMAL_PCR','BISULFITE_PCR']


class Condition(object):
    def __init__(self, weight, optimal, minimum, maximum):
        self.optimal = optimal
        self.weight = weight
        self.minimum = minimum
        self.maximum = maximum
    def score(self, value):
        return self.weight * abs(value-self.optimal)
    def bound(self, value):
        return self.minimum <= value <= self.maximum

class PrimerCondition(object):
    def __init__(self):
        self.primer_length = Condition(0.5, 23., 10, 30)
        self.gc        = Condition(1.0, 50., 30, 70)
        self.tm            = Condition(1.0, 60., 45, 70)
        self.sa            = Condition(0.1, 0., 0, 20)
        self.sea           = Condition(0.2, 0., 0, 10)
        self.pa            = Condition(0.1, 0., 0, 20)
        self.pea           = Condition(0.2, 0., 0, 10)

    def bound_primer(self, p):
        return self.primer_length.bound(len(p)) \
              and self.gc.bound(p.gc_ratio) \
              and self.tm.bound(p.melting_temperature()) \
              and self.sa.bound(p.sa.score) \
              and self.sea.bound(p.sea.score)

    def bound_primerpair(self, pp):
        return self.bound_primer(pp.fw) \
            and self.bound_primer(pp.rv) \
            and self.pa.bound(pp.pa.score) \
            and self.pea.bound(pp.pea.score)
        
    def score_primer(self, p):
        return self.primer_length.score(len(p)) \
              + self.gc.score(p.gc_ratio) \
              + self.tm.score(p.melting_temperature()) \
              + self.sa.score(p.sa.score) \
              + self.sea.score(p.sea.score)

    def score_primerpair(self, pp):
        return self.score_primer(pp.fw) \
            + self.score_primer(pp.rv) \
            + self.pa.score(pp.pa.score) \
            + self.pea.score(pp.pea.score)

    

class PrimerConditionBisulfite(PrimerCondition):
    def __init__(self):
        super(PrimerConditionBisulfite,self).__init__()
        self.gc        = Condition(1.0, 30., 0, 60)

NORMAL_PCR = PrimerCondition()
BISULFITE_PCR = PrimerConditionBisulfite()
