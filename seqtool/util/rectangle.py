### Line


class Line:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __repr__(self):
        return 'Line({},{})'.format(self.start, self.end)

    @property
    def length(self):
        return self.end - self.start

    @property
    def valid(self):
        return self.start < self.end

    def has_intersect(self, rhs):
        return self.intersect(rhs).valid

    def has_union(self, rhs):
        return self.union(rhs).valid

    def intersect(self, rhs):
        """
        >>> Line(0,100).intersect(Line(20,400))
        Line(20,100)
        """
        return Line(max(self.start,rhs.start), min(self.end,rhs.end))

    def union(self, rhs):
        """
        >>> Line(0,100).union(Line(20,400))
        Line(0,400)
        """
        return Line(min(self.start,rhs.start), max(self.end,rhs.end))

    def contains(self, x):
        """
        >>> Line(0,100).contains(50)
        True
        >>> Line(0,100).contains(-10)
        False
        """
        return self.start <= x <= self.end

    def translate(self, x):
        return Line(self.start + x, self.end + x)

    def scale(self, s):
        return Line(self.start * s, self.end * s)

### Rectangle


class Rectangle:
    def __init__(self, x0, x1, y0, y1):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1

    @property
    def x(self):
        return Line(self.x0, self.x1)

    @property
    def y(self):
        return Line(self.y0, self.y1)

    @property
    def valid(self):
        return self.x.valid and self.y.valid

    @property
    def width(self):
        return self.x.length

    @property
    def height(self):
        return self.y.length

    @classmethod
    def from_line(cls, x, y):
        return cls(x.start, x.end, y.start, y.end)

    @classmethod
    def union_all(cls, iterable):
        """
        >>> Rectangle.union_all([])
        Rectangle(0,0,0,0)
        >>> Rectangle.union_all([Rectangle(0,10,100,300)])
        Rectangle(0,10,100,300)
        >>> Rectangle.union_all([Rectangle(10,20,100,300), Rectangle(5,100,0,400)])
        Rectangle(5,100,0,400)
        """
        ret = Rectangle(0,0,0,0)
        try:
            it = iter(iterable)
            ret = next(it)
            for i in it:
                ret = ret.union(i)
        except TypeError:
            pass
        except StopIteration:
            pass
        return ret


    def __repr__(self):
        return 'Rectangle({},{},{},{})'.format(self.x0, self.x1, self.y0, self.y1)

    def __bool__(self):
        return self.valid

    def __mul__(self, rhs):
        return self.intersect(rhs)

    def __add__(self, rhs):
        return self.union(rhs)

    def include_point(self, x, y):
        return Rectangle( min(self.x0, x), max(self.x1, x), min(self.y0, y), max(self.y1, y) )

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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
