from collections import OrderedDict

class NamedList(object):
    def __init__(self, get_name=lambda x:x.name):
        self._dict = OrderedDict()
        self._get_name = get_name

    def append(self, value):
        name = self._get_name(value)
        if name in self._dict:
            raise ValueError('already exist: %s'%name)
        self._dict[name] = value

    def __len__(self):
        return len(self._dict)

    def __getitem__(self, name):
        v = self.get(name)
        if v:
            return v
        else:
            raise KeyError("No such item: {0}".format(name))

    def get(self, name):
        v = self._dict.get(name)
        if v:
            return v
        else:
            return None

    def __iter__(self):
        return iter(self._dict.values())

    def items(self):
        return iter(self._dict.items())

class DefaultNamedList(NamedList):
    def __init__(self, default, get_name=lambda x:x.name):
        """
        class Exsample(object):
            def __init__(self, name):
                self.name = name
        DefaultNamedList(Example)
        """
        self._default = default
        super(DefaultNamedList, self).__init__(get_name)

    def get(self, name):
        v = self._dict.get(name)
        if v:
            return v
        else:
            d = self._default(name)
            self._dict[name] = d
            return d

