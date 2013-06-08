

from ..util import xmlwriter
import itertools
from io import StringIO
from functools import reduce

DEBUG = False

SVG_HEADER = '''<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" 
         "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
'''

SVG_CSS = '''
line{
    stroke: black;
    stroke-width: 1;
}
path{
    stroke: black;
    stroke-width: 1;
}
text{
    fill: black;
    font: 10pt Monaco;
}
.graphline{
    stroke: red;
    stroke-width:0.5;
    stroke-dasharray:30,10;
}
'''

def ffmt(x):
    return '{:.2f}'.format(x) if x%1 else '{}'.format(int(x))

### Rectangle

class Rectangle(object):
    def __init__(self, x0, x1, y0, y1):
        self.valid = (x0 <= x1 and y0 <= y1)
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.width = x1-x0
        self.height = y1-y0

    def __repr__(self):
        return 'Rectangle(%f,%f,%f,%f)' % (self.x0, self.x1, self.y0, self.y1)

    def __bool__(self):
        return self.valid

    def __mul__(self, rhs):
        return self.intersect(rhs)

    def __add__(self, rhs):
        return self.add(rhs)

    def include_point(self, x, y):
        return Rectangle( min(self.x0, x),max(self.x1, x),min(self.y0, y),max(self.y1, y) )

    def add(self, rhs):
        return Rectangle(
            min(self.x0,rhs.x0),
            max(self.x1, rhs.x1),
            min(self.y0, rhs.y0),
            max(self.y1, rhs.y1) )

    def contains(self, x, y):
        return self.x0 <= x <= self.x1 and self.y0 <= y <= self.y1

    def intersect(self, rhs):
        return Rectangle(max(self.x0,rhs.x0),
            min(self.x1, rhs.x1),
            max(self.y0, rhs.y0),
            min(self.y1, rhs.y1))

    def intersects(self, rhs):
        return not ( self.x0 > rhs.x1 or self.x1 < rhs.x0
            or self.y0 > rhs.y1 or self.y1 < rhs.y0 )

    def translate(self,x,y):
        return Rectangle(self.x0+x, self.x1+x, self.y0+y, self.y1+y)

    def scale(self, x, y):
        return Rectangle(self.x0 * x, self.x1 * x, self.y0 * y, self.y1 * y)

    def intersects_x(self, rhs):
        return not ( self.x0 > rhs.x1 or self.x1 < rhs.x0 )

    def intersects_y(self, rhs):
        return not ( self.y0 > rhs.y1 or self.y1 < rhs.y0 )

### Base

class SvgBase(object):
    def __init__(self):
        self._rect = Rectangle(0,0,0,0)

    @property
    def rect(self):
        return self._rect

    def draw(self, b):
        pass

    def _style(self, kwargs):
        ret = {}
        #ret['stroke'] = 'black'
        for key, value in list(kwargs.items()):
            if key=='color':
                ret['stroke'] = value
                ret['style'] = 'fill:%s;'%value
            else:
                ret[key] = value
        return ret

    def svg_node(self, width=None, height=None):
        w = self.rect.width
        h = self.rect.height
        width = width or w
        height = height or h

        buff = StringIO()
        writer = xmlwriter.XmlWriter(buff)
        b = xmlwriter.builder(writer)
        with b.svg(xmlns="http://www.w3.org/2000/svg",
                   width=ffmt(width), height=ffmt(height),
#                   viewBox="0 0 %d %d"%(w, h),
                   preserveAspectRatio='none',
                   **{"xmlns:xlink":"http://www.w3.org/1999/xlink"}):
            with b.style(type='text/css'):
                writer.write('<![CDATA[')
                writer.write(self.svg_css())
                writer.write(']]>')
            self.draw(b)
        return buff.getvalue()

    def svg(self, width=None, height=None):
        return SVG_HEADER + self.svg_node(width,height)

    def svg_css(self):
        return SVG_CSS

class SvgNone(SvgBase):
    def __init__(self):
        super(SvgNone,self).__init__()

### wrappers

def g_translate(b,x,y):
    return b.g(transform='translate(%.2f,%.2f)'%(x,y))
def g_scale(b,x,y):
    return b.g(transform='scale(%.2f,%.2f)'%(x,y))
def g_translate_scale(b,x,y,sx,sy):
    return b.g(transform='translate(%.2f,%.2f) scale(%.2f,%.2f)'%(x,y,sx,sy))

class SvgTranslate(SvgBase):
    def __init__(self, x, y, child):
        self._rect = child.rect.translate(x,y)
        self.x = x
        self.y = y
        self.child = child
    
    def draw(self, b):
        with g_translate(b,self.x,self.y):
            self.child.draw(b)

class SvgScale(SvgBase):
    def __init__(self, x, y, child):
        self._rect = child.rect.scale(x,y)
        self.x = x
        self.y = y
        self.child = child
    
    def draw(self, b):
        with g_scale(b, self.x,self.y):
            self.child.draw(b)

class SvgTranslateScale(SvgBase):
    def __init__(self, x, y, sx, sy, child):
        self._rect = child.rect.translate(x,y).scale(sx,sy)
        self.x = x
        self.y = y
        self.sx = sx
        self.sy = sy
        self.child = child
    
    def draw(self, b):
        with g_translate_scale(b, self.x,self.y,self.sx,self.sy):
            self.child.draw(b)

class SvgPadding(SvgBase):
    def __init__(self, x, y, child):
        self._rect = Rectangle(child.rect.x0, child.rect.x1+2*x,child.rect.y0, child.rect.y1+2*y)
        self.x = x
        self.y = y
        self.child = child
    
    def draw(self, b):
        with g_translate(b, self.x, self.y):
            self.child.draw(b)

### Containers

class SvgItems(SvgBase):
    def __init__(self):
        super(SvgItems,self).__init__()
        self.items = []

    def add(self, item):
        self.items.append(item)

    @property
    def rect(self):
        return reduce(lambda a,b: a+b, (i.rect for i in self.items)) if self.items else Rectangle(0,0,0,0)

    def draw(self, b):
        for c in self.items:
            c.draw(b)

class SvgItemsFixedHeight(SvgItems):
    def __init__(self, height):
        super(SvgItemsFixedHeight,self).__init__()
        self.height = height

    @property
    def rect(self):
        if self.items:
            r = reduce(lambda a,b: a+b, (i.rect for i in self.items))
            r.y1 = max(r.y1, self.height)
            return r
        else:
            return Rectangle(0,0,0,self.height)

class SvgItemsVStack(SvgItems):
    """
    each items are located at different lines.
    +----------------+
    | item1          |
    |                |
    +--+-------------+---+
       | item2           |
    +--+----------+------+
    | item3       |
    +-------------+
    """
    def __init__(self):
        super(SvgItemsVStack, self).__init__()

    def _item_heights(self):
        return sum(item.rect.height for item in self.items)

    @property
    def rect(self):
        x0 = min(item.rect.x0 for item in self.items)
        x1 = max(item.rect.x1 for item in self.items)
        h = sum(item.rect.y1 for item in self.items)
        #h = sum(item.rect.height for item in self.items)
        return Rectangle(x0,x1,0,h)

    def draw(self, b):
        y = 0
        for item in self.items:
            with g_translate(b, 0, y):
                item.draw(b)
            y += item.rect.height

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(fillvalue=fillvalue, *args)

class SvgItemsVFree(SvgItems):
    """
    move items along Y axis automaticaly not to intersects other items.

                                +-------+
            +-------+           | item3 |
            | item1 |           +-------+
            +-------+       +--------+
                            | item2  |
                            +--------+
    """
    def __init__(self):
        super(SvgItemsVFree,self).__init__()
        self.trans_y = []

    def _find_freespace(self, rects, rect):
        for r in rects:
            pass
            #print r, rect, r.intersects_x(rect)
        conflicts = [(r.y0, r.y1) for r in rects if r.intersects_x(rect)]
        conflicts.sort()
        if not conflicts:
            return 0
        lconflicts = [0] + list(itertools.chain.from_iterable(conflicts))

        y = 0
        for y0, y1 in grouper(2, lconflicts, None):
            y = max(y0, y)
            if y1==None:
                return y

            h = y1-y

            if h < 0:
                continue

            if h >= rect.height:
                return y
        raise Exception("Implemantation Error")

    def _calc_rect(self):
        rects = []
        trans_ys = []
        for item in self.items:
            y0 = self._find_freespace(rects, item.rect)
            trans_y = y0 - item.rect.y0
            rects.append(item.rect.translate(0, trans_y))
            trans_ys.append(trans_y)

        return reduce(Rectangle.__add__,rects,Rectangle(0,0,0,0)), trans_ys

    @property
    def rect(self):
        rect, tys = self._calc_rect()
        return rect

    def draw(self, b):
        rect, tys = self._calc_rect()

        for i,item in enumerate(self.items):
            with g_translate(b, 0, tys[i]):
                item.draw(b)

class SvgMatrix(SvgBase):
    def __init__(self):
        super(SvgMatrix,self).__init__()
        self.rows = []
        self.aligns = []

    def add_row(self, fixeditems, align=None):
        """fixed items are list. default item is SvgBase() if item is None"""
        replaced = [x if x else SvgNone() for x in fixeditems]
        self.rows.append(replaced)
        self.aligns.append(align)

    def _width_height(self):
        #num_rows = len(self.rows)
        num_columns = 0
        for row in self.rows:
            num_columns = max(num_columns, len(row))
        col_widths = [0 for x in range(num_columns)]
        row_heights = []

        for row in self.rows:
            row_heights.append(max(c.rect.y1 for c in row))
            for i,c in enumerate(row):
                col_widths[i] = max(col_widths[i], c.rect.x1)
        return col_widths, row_heights

    @property
    def rect(self):
        col_w, row_h = self._width_height()
        return Rectangle(0, sum(col_w), 0, sum(row_h))

    def get_align(self, ix, iy):
        try:
            return self.aligns[iy][ix] or 'left'
        except:
            return 'left'

    def draw(self, b):
        col_w, row_h,  = self._width_height()

        iy = 0
        y = 0
        for row in self.rows:
            x = 0
            ix = 0
            for cell in row:
                align = self.get_align(ix,iy)
                if align=='right' or align=='r':
                    tx = col_w[ix] - cell.rect.width
                elif align=='center' or align=='c':
                    tx = (col_w[ix] - cell.rect.width) / 2
                else:
                    tx = 0

                if not isinstance(cell, SvgNone):
                    with g_translate(b,x+tx,y):
                        cell.draw(b)
                        #b.test(repr(cell))
                x += col_w[ix]
                ix += 1
            y += row_h[iy]
            iy += 1

### Item

class SvgLine(SvgBase):
    def __init__(self, x0, x1, y0, y1, **kwargs):
        self.x0 = x0
        self.x1 = x1
        self.y0 = y0
        self.y1 = y1
        self.kwargs = self._style(kwargs)
        self._rect = Rectangle(x0, x1, y0, y1)

    def draw(self, b):
        b.line(x1=ffmt(self.x0), x2=ffmt(self.x1), y1=ffmt(self.y0), y2=ffmt(self.y1), **self.kwargs)

class SvgVline(SvgLine):
    def __init__(self, x, y0, y1, **kwargs):
        super().__init__(x,x,y0,y1,**kwargs)

class SvgHline(SvgLine):
    def __init__(self, x0, x1, y, **kwargs):
        super().__init__(x0,x1,y,y,**kwargs)

class SvgRect(SvgBase):
    def __init__(self, x, y, width, height, **kwargs):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.kwargs = self._style(kwargs)
        self._rect = Rectangle(x, x+width, y, y+height)

    def draw(self,b):
        b.rect(x=ffmt(self.x), y=ffmt(self.y), width=ffmt(self.width), height=ffmt(self.height), **self.kwargs)

class SvgLines(SvgBase):
    def __init__(self, lines, **kwargs):
        self.lines = lines
        self.kwargs = self._style(kwargs)

        x0,y0,x1,y1 = self.lines[0]
        for xx0,yy0,xx1,yy1 in self.lines[1:]:
            x0 = min(x0, xx0, xx1)
            y0 = min(y0, yy0, yy1)
            x1 = max(x1, xx0, xx1)
            y1 = max(y1, yy0, yy1)
        self._rect = Rectangle(x0, x1, y0, y1)

    def draw(self, b):
        d = " ".join("M{} {} l{} {}".format(ffmt(x0),ffmt(y0),ffmt(x1-x0),ffmt(y1-y0)) for x0,y0,x1,y1 in self.lines)
        b.path(d=d, **self.kwargs)

"""
class SvgPath(SvgBase):
    def __init__(self, x, y, **kwargs):
        self.d = []
        self.kwargs = self._style(kwargs)
        self._rect = Rectangle(x, x, y, y)
        self.add_p("M", x,y)

    def add_p(self, i, x, y):
        self.d.append((i, [x,y]))
        self._rect = self._rect.include_point(x,y)

    def m(self, x, y):
        self.add_p('m',x,y)
    def M(self, x, y):
        self.add_p('M',x,y)
    def L(self, x, y):
        self.add_p('L',x,y)
    def l(self, x, y):
        self.add_p('l',x,y)

    def z(self):
        self.d.append(('z', []))

    # TODO other instruction supports

    def draw(self,b):
        d = " ".join("{} {}".format(i, ' '.join(v)) for i, v in self.d)
        b.path(d=d, **self.kwargs)

    def scale(self, sx, sy):
        ret = SvgPath()
        for i, (x,y) in self.d:
            ret.add_p(i, x/sx, y/sy)
        return ret
"""

class SvgCircle(SvgBase):
    def __init__(self, x, y, radius, **kwargs):
        self.x = x
        self.y = y
        self.radius = radius
        self.kwargs = self._style(kwargs)
        self._rect = Rectangle(x-radius, x+radius, y-radius, y+radius)

    def draw(self,b):
        b.circle(cx=ffmt(self.x), cy=ffmt(self.y), r=ffmt(self.radius), **self.kwargs)

class SvgHbar(SvgRect):
    def __init__(self, start, end, y, thick, **kwargs):
        super().__init__(start, y-thick/2, end-start, thick, **kwargs)

DEFAULT_FONTSIZE = 12

def font_width(fontsize=DEFAULT_FONTSIZE):
    return fontsize*0.6
def font_height(fontsize=DEFAULT_FONTSIZE):
    return fontsize*1.1

class SvgText(SvgBase):
    def __init__(self, text, x, y, anchor='start', fontsize=DEFAULT_FONTSIZE):
        self.text = text

        self.w = font_width(fontsize) * len(text)
        self.h = font_height(fontsize)

        # dont use svg attribute, anchor='middle'
        if anchor == 'middle':
            x -= self.w/2.

        # Firefox dont support textLength. so you need to tell all the location of each characters respectively.
        # self.x = ' '.join(str(x+i*font_width(fontsize)) for i in range(len(text)))
        self.x = x
        self.y = y
        self.textLength = self.w

        self.style = {}
        if fontsize != DEFAULT_FONTSIZE:
            self.style['fontsize'] = fontsize

        self._rect = Rectangle(x, x+self.w, y, y+self.h)

    def draw(self,b):
        with b['text'](x=ffmt(self.x), y=ffmt(self.y+self.h), textLength=ffmt(self.textLength), **self.style):
            b.text(self.text)
        if DEBUG:
            b.rect(x=ffmt(self.rect.x0), y=ffmt(self.rect.y0), width=ffmt(self.rect.width), height=ffmt(self.rect.height), stroke='red', style='fill:none;')

class SvgGraphline(SvgItemsFixedHeight):
    def __init__(self, height, values, bars=[], scalex=1., width=None, **kwargs):
        super(SvgGraphline,self).__init__(height)

        step = max(1, int(1./scalex)) # for smaller file size.
        self.points = [(i*scalex,height*(1.-values[i])) for i in range(0, len(values), step)]
        self.bars = bars
        self.kwargs = self._style(kwargs)
        self.width = width or len(values)

        self._rect = Rectangle(0, len(values), 0, height)

    def draw(self, b):
        p = ','.join(["%.2f %.2f"%(x,y) for (x,y) in self.points])
        b.polyline(points=p, **self.kwargs)
        for bar in self.bars:
            t = self.height * (1.-bar)
            b.line(x1=0, x2=ffmt(self.width), y1=ffmt(t), y2=ffmt(t), klass='graphline')


class SvgItemGenerator(object):
    def __init__(self, scalex, scaley):
        self.sx = scalex
        self.sy = scaley

    def line(self, x0, x1, y0, y1, **kwargs):
        x0 /= self.sx
        x1 /= self.sx
        y0 /= self.sy
        y1 /= self.sy
        return SvgLine(x0, x1, y0, y1, **kwargs)

    def lines(self, lines, **kwargs):
        lines = [(x0/self.sx, y0/self.sy, x1/self.sx, y1/self.sy) for x0,y0,x1,y1 in lines]
        return SvgLines(lines, **kwargs)

    def vline(self, x, y0, y1, **kwargs):
        x /= self.sx
        y0 /= self.sy
        y1 /= self.sy
        return SvgVline(x, y0, y1, **kwargs)

    def hline(self, x0, x1, y, **kwargs):
        x0 /= self.sx
        x1 /= self.sx
        y /= self.sy
        return SvgHline(x0, x1, y, **kwargs)

    def rect(self, x, y, width, height, **kwargs):
        x /= self.sx
        width /= self.sx
        y /= self.sy
        height /= self.sy
        return SvgRect(x, y, width, height, **kwargs)

    def hbar(self, start, end, y, thick, **kwargs):
        start /= self.sx
        end /= self.sx
        y /= self.sy
        thick /= self.sy
        return SvgHbar(start, end, y, thick, **kwargs)

    def text(self, text, x=0, y=0, fontsize=DEFAULT_FONTSIZE, anchor='start'):
        x /= self.sx
        y /= self.sy
        return SvgText(text, x, y, fontsize=fontsize, anchor=anchor)

    def graphline(self, height, values, bars=[], width=None, **kwargs):
        height /= self.sy
        return SvgGraphline(height, values, bars, scalex=self.sx, width=width, **kwargs)



