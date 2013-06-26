

from ..util import xmlwriter
import itertools
from io import StringIO
from functools import reduce

from .rectangle import Line, Rectangle

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
rect{
    fill: none;
    stroke: black;
}
'''

def fm(x):
    """
    >>> fm(0.5)
    0.5
    >>> fm(0.0555)
    0.05
    >>> fm(45.)
    45
    """
    return '{:.2f}'.format(x) if x%1 else '{}'.format(int(x))


### Base

class SvgBase(object):
    def __init__(self):
        self._rect = Rectangle(0,0,0,0)

    @property
    def rect(self):
        return self._rect

    def draw(self, b):
        pass

    def draw_defs(self, b):
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
        """
        For embeded svg
        """
        w = self.rect.width
        h = self.rect.height
        width = width or w
        height = height or h

        buff = StringIO()
        writer = xmlwriter.XmlWriter(buff)
        b = xmlwriter.builder(writer)
        with b.svg(xmlns="http://www.w3.org/2000/svg",
                   width=fm(width), height=fm(height),
#                   viewBox="0 0 %d %d"%(w, h),
                   preserveAspectRatio='none',
                   **{"xmlns:xlink":"http://www.w3.org/1999/xlink"}):
            with b.style(type='text/css'):
                writer.write('<![CDATA[')
                writer.write(self.svg_css())
                writer.write(']]>')
            with b.defs:
                self.draw_defs(b)
            self.draw(b)
        return buff.getvalue()

    def svg(self, width=None, height=None):
        return SVG_HEADER + self.svg_node(width,height)

    def svg_css(self):
        return SVG_CSS

class SvgNone(SvgBase):
    def __init__(self):
        super().__init__()

### wrappers

def g_translate(b,x,y):
    return b.g(transform='translate(%.2f,%.2f)'%(x,y))
def g_scale(b,x,y):
    return b.g(transform='scale(%.2f,%.2f)'%(x,y))
def g_translate_scale(b,x,y,sx,sy):
    return b.g(transform='translate(%.2f,%.2f) scale(%.2f,%.2f)'%(x,y,sx,sy))

class SvgParent(SvgBase):
    def __init__(self, children=[]):
        self.children = children
        super().__init__()

    def draw_defs(self, b):
        for c in self.children:
            c.draw_defs(b)

    def draw(self, b):
        for c in self.children:
            c.draw(b)

    @property
    def rect(self):
        return Rectangle.union_all(i.rect for i in self.children)

class SvgParentSingle(SvgParent):
    def __init__(self, child):
        super().__init__([child])
        self.child = child

    @property
    def rect(self):
        return self.child.rect

class SvgTranslate(SvgParentSingle):
    def __init__(self, x, y, child):
        super().__init__(child)
        self.x = x
        self.y = y
    
    def draw(self, b):
        with g_translate(b,self.x,self.y):
            super().draw(b)

    @property
    def rect(self):
        return super().rect.translate(self.x, self.y)

class SvgClipWidth(SvgParentSingle):
    num = 0
    def __init__(self, x0, x1, child):
        super().__init__(child)
        self.x0 = x0
        self.x1 = x1
        self.num = SvgClipWidth.num
        self.id = 'clip{}'.format(self.num)
        SvgClipWidth.num += 1

    def draw_defs(self, b):
        if self.rect.valid:
            with b.clipPath(id=self.id):
                b.rect(x=fm(self.rect.x0-1), y=fm(self.rect.y0-1), width=fm(self.rect.width+2), height=fm(self.rect.height+2))
        super().draw_defs(b)

    def draw(self, b):
        if self.rect.valid:
            with b.g(**{'clip-path':'url(#{})'.format(self.id)}):
                super().draw(b)

    @property
    def rect(self):
        r = super().rect
        return Rectangle(max(self.x0,r.x0), min(self.x1,r.x1), r.y0, r.y1)

class SvgScale(SvgParentSingle):
    def __init__(self, x, y, child):
        super().__init__(child)
        self.x = x
        self.y = y
    
    def draw(self, b):
        with g_scale(b, self.x,self.y):
            super().draw(b)

    @property
    def rect(self):
        return super().rect.scale(self.x, self.y)

class SvgTranslateScale(SvgParentSingle):
    def __init__(self, x, y, sx, sy, child):
        super().__init__(child)
        self.x = x
        self.y = y
        self.sx = sx
        self.sy = sy
    
    def draw(self, b):
        with g_translate_scale(b, self.x,self.y,self.sx,self.sy):
            super().draw(b)

    @property
    def rect(self):
        return super().rect.translate(self.x,self.y).scale(self.sx,self.sy)

class SvgPadding(SvgParentSingle):
    def __init__(self, x, y, child):
        super().__init__(child)
        self.x = x
        self.y = y
    
    def draw(self, b):
        with g_translate(b, self.x, self.y):
            super().draw(b)

    @property
    def rect(self):
        r = super().rect
        return Rectangle(r.x0, r.x1+2*self.x,r.y0, r.y1+2*self.y)

class SvgBoundbox(SvgParentSingle):
    def __init__(self, child, **kwargs):
        super().__init__(child)
        self.kwargs = self._style(kwargs)

    def draw(self,b):
        r = self.rect
        b.rect(x=fm(r.x0), y=fm(r.y0), width=fm(r.width), height=fm(r.height), **self.kwargs)
        super().draw(b)

class SvgExpandWidth(SvgParentSingle):
    def __init__(self, width, child, **kwargs):
        super().__init__(child)
        self.kwargs = self._style(kwargs)
        self.expand = width

    @property
    def rect(self):
        r = super().rect
        r.x1 = r.x0 + max(r.width, self.expand)
        return r

class SvgExpandHeight(SvgParentSingle):
    def __init__(self, height, child, **kwargs):
        super().__init__(child)
        self.kwargs = self._style(kwargs)
        self.expand = height

    @property
    def rect(self):
        r = super().rect
        r.y1 = r.y0 + max(r.height, self.expand)
        return r

### Containers

class SvgItems(SvgParent):
    def __init__(self):
        super().__init__([])

    def add(self, item):
        self.children.append(item)

class SvgItemsFixedHeight(SvgItems):
    def __init__(self, height):
        super().__init__()
        self.height = height

    @property
    def rect(self):
        if self.children:
            r = Rectangle.union_all(i.rect for i in self.children)
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
        super().__init__()

    def _item_heights(self):
        return sum(item.rect.height for item in self.children)

    @property
    def rect(self):
        if not self.children:
            return Rectangle(0,0,0,0)
        x0 = min(item.rect.x0 for item in self.children)
        x1 = max(item.rect.x1 for item in self.children)
        h = sum(item.rect.y1 for item in self.children)
        #h = sum(item.rect.height for item in self.children)
        return Rectangle(x0,x1,0,h)

    def draw(self, b):
        y = 0
        for item in self.children:
            with g_translate(b, 0, y):
                item.draw(b)
            y += item.rect.height

def grouper(n, iterable, fillvalue=None):
    """
    >>> grouper(3, 'ABCDEFG', 'x')
    "ABC DEF Gxx"
    """
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
    def __init__(self, reverse=False):
        super().__init__()
        self.reverse=reverse

    def _find_freespace(self, rects, rect):
        conflicts = [(r.y0, r.y1) for r in rects if r.x.has_intersect(rect.x)]
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
        for item in self.children:
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
        height = rect.y1

        for i,item in enumerate(self.children):
            y = tys[i]
            if not self.reverse:
                yy = y
            else:
                yy = height-(y+item.rect.y1)
            with g_translate(b, 0, yy):
                item.draw(b)

class SvgMatrix(SvgBase):
    def __init__(self):
        super().__init__()
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

    def draw_defs(self, b):
        for row in self.rows:
            for cell in row:
                if not isinstance(cell, SvgNone):
                    cell.draw_defs(b)

### Wrapping

def iter_step(width, start, end):
    for x in range(start, end, width):
        yield x, min(x+width, end)

class SvgItemsWrapping(SvgParentSingle):
    unique_counter = 0

    def __init__(self, width, child):
        super().__init__(child)

        self.width = width
        self.id = 'siw_{}'.format(SvgItemsWrapping.unique_counter)
        SvgItemsWrapping.unique_counter += 1

    @property
    def length(self):
        return super().rect.x1

    @property
    def height(self):
        return super().rect.y1

    @property
    def rect(self):
        l = self.length
        h = self.height
        c = l / self.width
        if l % self.width:
            c += 1
        return Rectangle(0, self.width, 0, c * h)

    def iter_step(self):
        c = 0
        y = 0
        h = self.height
        for p,q in iter_step(int(self.width), 0, int(self.length)):
            r = Rectangle(0, self.width, y, y+h)
            yield p,q,r,'siw_{}_{}'.format(p,q)
            c += 1
            y += h

    def draw_defs(self, b):
        super().draw_defs(b)

        with b.g(id=self.id):
            super().draw(b)

        for p,q,r,bid in self.iter_step():
            with b.clipPath(id=bid):
                b.rect(x=fm(r.x0), y=fm(r.y0), width=fm(r.width), height=fm(r.height))

    def draw(self, b):
        for p,q,r,bid in self.iter_step():
            #with b.g(**{'clip-path':'url(#{})'.format(bid)}):
            with g_translate(b, -p, r.y0):
                #super().draw(b)
                b.use(**{'xlink:href':'#{}'.format(self.id)})

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
        b.line(x1=fm(self.x0), x2=fm(self.x1), y1=fm(self.y0), y2=fm(self.y1), **self.kwargs)

class SvgVline(SvgLine):
    def __init__(self, x, y0, y1, **kwargs):
        super().__init__(x,x,y0,y1,**kwargs)

class SvgHline(SvgLine):
    def __init__(self, x0, x1, y, **kwargs):
        super().__init__(x0,x1,y,y,**kwargs)

class SvgHlineBox(SvgExpandHeight):
    def __init__(self, x0, x1, height):
        super().__init__(height, SvgHline(x0, x1, height/2.))


class SvgRect(SvgBase):
    def __init__(self, x, y, width, height, **kwargs):
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self.kwargs = self._style(kwargs)
        self._rect = Rectangle(x, x+width, y, y+height)

    def draw(self,b):
        b.rect(x=fm(self.x), y=fm(self.y), width=fm(self.width), height=fm(self.height), **self.kwargs)

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
        d = " ".join("M{} {} l{} {}".format(fm(x0),fm(y0),fm(x1-x0),fm(y1-y0)) for x0,y0,x1,y1 in self.lines)
        b.path(d=d, **self.kwargs)

class SvgCircle(SvgBase):
    def __init__(self, x, y, radius, **kwargs):
        self.x = x
        self.y = y
        self.radius = radius
        self.kwargs = self._style(kwargs)
        self._rect = Rectangle(x-radius, x+radius, y-radius, y+radius)

    def draw(self,b):
        b.circle(cx=fm(self.x), cy=fm(self.y), r=fm(self.radius), **self.kwargs)

class SvgHbar(SvgRect):
    def __init__(self, start, end, y, thick, **kwargs):
        super().__init__(start, y-thick/2, end-start, thick, **kwargs)

DEFAULT_FONTSIZE = 12

def font_width(fontsize=DEFAULT_FONTSIZE):
    return fontsize*0.6
def font_height(fontsize=DEFAULT_FONTSIZE):
    return fontsize*1.1

class SvgText(SvgBase):
    def __init__(self, text, x=0, y=0, anchor='start', fontsize=DEFAULT_FONTSIZE):
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
        with b['text'](x=fm(self.x), y=fm(self.y+self.h), textLength=fm(self.textLength), **self.style):
            b.text(self.text)
        if DEBUG:
            b.rect(x=fm(self.rect.x0), y=fm(self.rect.y0), width=fm(self.rect.width), height=fm(self.rect.height), stroke='red', style='fill:none;')

class SvgGraphline(SvgItemsFixedHeight):
    def __init__(self, height, values, bars=[], scalex=1., width=None, **kwargs):
        super().__init__(height)

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
            b.line(x1=0, x2=fm(self.width), y1=fm(t), y2=fm(t), klass='graphline')


class SvgItemGenerator:
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



