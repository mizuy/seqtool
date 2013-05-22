
class Templates(object):
    def __init__(self):
        self._list = NamedList()

    def append(self, name, seq):
        self._list[name] = seq

    def items(self):
        return list(self._list.items())


class Template(object):
    def bisulfite_conversion(self):
        BisulfiteTemplate(self.name, self.seq)

class BisulfiteTemplates(Templates):
    def __init__(self, name, seq):
        self.seq = seq
        self.sense = bisulfite_conversion_ambigusous(self.template, sense=True)
        self.asense = bisulfite_conversion_ambigusous(self.template, sense=False)

    def __iter__(self):
        yield self.name, self.seq
        yield self.name+' BS+', self.sense
        yield self.name+' BS-', self.asense

class RTTemplates(object):
    pass
