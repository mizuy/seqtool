### Line


def Line:
    def __init__(self, start, end):
        self.valid = (self.start <= self.end)
        self.start = start
        self.end = end

    def __repr__(self):
        return 'Line(%f,%f)' % (self.start, self.end)

    @property
    def length(self):
        return self.start - self.end

    def has_intersect(self, rhs):
        return self.intersect(rhs).valid

    def has_union(self, rhs):
        return self.union(rhs).valid

    def intersect(self, rhs):
        return Line(max(self.start,rhs.start), min(self.end,rhs.end))

    def union(self, rhs):
        return Line(min(self.start,rhs.start), max(self.end,rhs.end))

    def contains(self, x):
        return self.start <= x <= self.end

    def translate(self, x):
        return Line(self.start + x, self.end + y)

    def scale(self, s):
        return Line(self.start * s, self.end * s)

### Rectangle

class Rectangle:
    def __init__(self, x0, x1, y0, y1):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.x = Line(x0, x1)
        self.y = Line(y0, y1)
        self.valid = self.x.valid and self.y.valid
        self.width = self.x.length
        self.height = self.y.length

    @classmethod
    def from_line(cls, x, y):
        return cls(x.start, x.end, y.start, y.end)

    @classmethod
    def union(cls, iterable):
        ret = Rectangle(0,0,0,0)
        try:
            i = iter(iterable)

        except TypeError:
        else:
            ret = iter.next()
            for i in iter:
                ret = ret.union(i)
        except StopIteration:
            continue
        return ret


    def __repr__(self):
        return 'Rectangle(%f,%f,%f,%f)' % (self.x0, self.x1, self.y0, self.y1)

    def __bool__(self):
        return self.valid

    def __mul__(self, rhs):
        return self.intersect(rhs)

    def __add__(self, rhs):
        return self.union(rhs)

    def include_point(self, x, y):
        return Rectangle( min(self.x0, x),max(self.x1, x), min(self.y0, y), max(self.y1, y) )

    def has_union(self, rhs):
        return self.union(rhs).valid

    def union(self, rhs):
        return Rectangle.from_line(self.x.union(rhs.x), self.y.union(rhs.y))

    def has_intersect(self, rhs):
        return self.intersect(rhs).valid

    def intersect(self, rhs):
        return Rectangle.from_line(self.x.intersect(rhs.x), self.y.intersect(rhs.y))

    def contains(self, x, y):
        return self.x.contains(x) and self.y.contains(y)

    def translate(self,x,y):
        return Rectangle.from_line(self.x.translate(x), self.y.translate(y))

    def scale(self, x, y):
        return Rectangle.from_line(self.x.scale(x), self.y.scale(y))
