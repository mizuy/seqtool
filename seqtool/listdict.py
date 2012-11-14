
class ListDict(object):
    def __init__(self):
        self._list = []
        self._dict = {}

    def append(self, value):
        if value.name in self._dict:
            raise ValueError('already exist: %s'%value.name)
        self._dict[value.name] = value
        self._list.append(value)

    def __len__(self):
        return len(self._list)

    def __getitem__(self, name):
        return self._dict[name]

    def __iter__(self):
        return iter(self._list)

