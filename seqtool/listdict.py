from collections import OrderedDict

class ListDict(object):
    def __init__(self):
        self._dict = OrderedDict()

    def append(self, value):
        if value.name in self._dict:
            raise ValueError('already exist: %s'%value.name)
        self._dict[value.name] = value

    def __len__(self):
        return len(self._dict)

    def __getitem__(self, name):
        return self._dict[name]

    def __iter__(self):
        return self._dict.itervalues()

