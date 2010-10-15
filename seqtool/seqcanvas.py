import Image, ImageDraw, ImageFont, sys
import math,os

def isnum(n):
    """ check if n is a number """
    return type(n) in (int, float, long)
class Point(list):
    def __init__(self, x, y):
        super(Point,self).__init__((x,y))

    def __pos__(self):
        return Point(self.x, self.y)
    def __neg__(self):
        return Point(-self.x, -self.y)
    def __add__(self, p):
        return Point(self.x+p.x, self.y+p.y)
    def __sub__(self,a):
        return self.__add__(-a)
    def __repr__(self):
        return "Point(%s,%s)"%(self.x,self.y)
    def __item__(self,i):
        if i==0:
            return self.x
        elif i==1:
            return self.y
        else:
            raise "Point: index error"
    def tuple(self):
        return (self.x, self.y)
    def __div__(self,v):
        if isnum(v):
            return Point(1.*self.x/v,1.*self.y/v)
        else:
            raise TypeError()
    def __mul__(self,v):
        if isnum(v):
            return Point(self.x*v,self.y*v)
        else:
            raise TypeError, 'Point must be multiplied with number, but, %s'%v

    def __getattr__(self, v):
        if v=='x':
            return self[0]
        elif v=='y':
            return self[1]
        else:
            raise IndexError()
    def map(self, f):
        return Point(f(self.x),f(self.y))

class CanvasBase(object):
    def set_height(self, height):
        raise NotImplementedError
    def draw_text(self, xy, text, pxy=Point(0,0), color='#000', bgcolor='#fff'):
        raise NotImplementedError
    def draw_line(self, xy0, xy1, pxy0=Point(0,0), pxy1=Point(0,0), color='#000'):
        raise NotImplementedError
    def draw_hbar(self, xy, style, color='#000'):
        raise NotImplementedError
    def draw_textbar(self, x0, x1, y, text='', textstyle=0, color='#000'):
        raise NotImplementedError

class CanvasImage(CanvasBase):
    font = ImageFont.truetype(os.path.join(os.path.dirname(__file__),'type_writer.ttf'), 8)
    
    @classmethod
    def font_size(cls):
        return cls.text_size('A')
    @classmethod
    def text_size(cls,text):
        return Point(*cls.font.getsize(text))

    def __init__(self, size, bgcolor='#fff'):
        self.font_size = CanvasImage.font_size()
        self.size = size
        self.bgcolor = bgcolor
        self.image = Image.new("RGB",self.size.map(int).tuple(),self.bgcolor)
        self.draw = ImageDraw.Draw(self.image)
    def save(self, filename, format):
        self.image.save(filename,format)
        
    def resize(self, newsize):
        if self.size == newsize:
            return
        self.size = newsize
        i = Image.new("RGB", newsize.map(int).tuple(), self.bgcolor)
        i.paste(self.image, (0,0))
        self.image = i
        self.draw = ImageDraw.Draw(self.image)

    def expand(self, y):
        newy = self.size.y
        while y+1 > newy:
            newy = 2*newy
        self.set_height(newy)

    def set_height(self, height):
        self.resize(Point(self.size.x, height))
        
    def draw_text(self, xy, text, pxy=Point(0,0), color='#000',bgcolor='#fff'):
        xy0 = xy+pxy
        self.expand(xy0.y)
        if bgcolor!=self.bgcolor:
            xy1 = xy0+self.text_size(text)
            self.draw.rectangle((xy0-Point(1,1)).tuple() + (xy1+Point(1,1)).tuple(), fill=bgcolor)
        self.draw.text(xy, text, font=self.font, fill=color)
        
    def draw_line(self, xy0, xy1, pxy0=Point(0,0), pxy1=Point(0,0), color='#000'):
        self.expand(xy0.y)
        self.expand(xy1.y)
        self.draw.line((xy0+pxy0).tuple()+(xy1+pxy1).tuple(),fill=color)

    def draw_hbar(self, xy, style, color='#000'):
        x,y = xy
        self.expand(y)
        if style==0: #'|'
            self.draw.line((x,y-3,x,y+3),fill=color)
        elif style<0: #'<'
            self.draw.line((x,y,x+3,y+3),fill=color)
            self.draw.line((x,y,x+3,y-3),fill=color)
        elif style>0: #'>'
            self.draw.line((x,y,x-3,y+3),fill=color)
            self.draw.line((x,y,x-3,y-3),fill=color)
        else:
            pass
            
    def draw_textbar(self, x0, x1, y, text='', textstyle=0, color='#000'):
        self.expand(y)
        '''
        draw a horizontal bar
        lstyle/rstyle: None or 1 or 2, 1 means ARROW, 2 means bar
        textstyle: -1 for text above the bar, 0 for text within the bar, 1 for text below the bar
        '''
        if text:
            ts = self.text_size(text)
            width,height = ts.x, ts.y
            c = (x1+x0)/2
            t0 = c-width/2
            t1 = c+width/2
            ty = y-height/2
            if textstyle==0:
                if x1-x0 >= width+2:
                    self.draw.text((t0,ty),text, font=self.font, fill=color)
                    self.draw.line((x0,y,t0-1,y), fill=color)
                    self.draw.line((t1+1,y,x1,y), fill=color)
                else:
                    # dont write too long text
                    self.draw.line((x0,y,x1,y),fill=color)
            else:
                tty = ty-height if textstyle<0 else ty+height
                self.draw.text((t0,tty),text, font=self.font, fill=color)
                self.draw.line((x0,y,x1,y),fill=color)
        else:
            self.draw.line((x0,y,x1,y),fill=color)

class CanvasRenderer(CanvasBase):
    def __init__(self):
        pass
    def set_target(self, renderer):
        self.renderer = renderer
    def transform(self, xy):
        raise NotImplementedError, "CanvasRenderer.transform must be implemented in derived classes"
    def set_height(self, height):
        a = self.transform(Point(0,height))
        return self.renderer.set_height(a.y)
        
    def draw_text(self, xy, text, pxy=Point(0,0), color='#000', bgcolor='#fff'):
        assert isinstance(xy, Point)
        assert isinstance(pxy, Point)
        return self.renderer.draw_text(self.transform(xy), text, pxy, color, bgcolor)
    def draw_line(self, xy0, xy1, pxy0=Point(0,0), pxy1=Point(0,0), color='#000'):
        assert isinstance(xy0, Point)
        assert isinstance(xy1, Point)
        assert isinstance(pxy0, Point)
        assert isinstance(pxy1, Point)
        return self.renderer.draw_line(self.transform(xy0), self.transform(xy1), pxy0, pxy1, color)

    def draw_hbar(self, xy, style, color='#000'):
        return self.renderer.draw_hbar(self.transform(xy), style, color)
    def draw_textbar(self, x0, x1, y, text='', textstyle=0, color='#000'):
        a = self.transform(Point(x0,y))
        b = self.transform(Point(x1,y))
        return self.renderer.draw_textbar(a.x, b.x, a.y, text, textstyle, color)

class CanvasRendererScaleTransform(CanvasRenderer):
    def __init__(self, scale_x, scale_y):
        self.scale_x = scale_x
        self.scale_y = scale_y
        super(CanvasRendererScaleTransform, self).__init__()

    def transform(self, xy):
        assert isinstance(xy, Point)
        return Point(self.scale_x*xy[0],self.scale_y*xy[1])

class SeqCanvasBase(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.vcolumn = end-start
        self.mcolumn = self.vcolumn
        self.renderer = None
        self.canvas = None

    def save(self,filename,format,height=None):
        if height:
            self.renderer.set_height(height)
        self.canvas.save(filename,format)

    def get_inner_index(self, xy):
        return xy-Point(self.start, 0)
    def get_mapped_index(self, i):
        return i

    def put_char(self, xy, char, color='#000', bgcolor='#fff'):
        self.renderer.draw_text(self.get_mapped_index(self.get_inner_index(xy)), char, color=color, bgcolor=bgcolor)

    def show_primer(self, x, y, primer, color, name='', forward=True):
        if forward:
            self.put_bar(x, x+len(primer)-1, y, color, primer.name, rarrow=True)
        else:
            self.put_bar(x, x+len(primer)-1, y, color, primer.name, larrow=True)

    def show_pcr(self, start, end, y, primer_fw, primer_rv, color, name=''):
        self.put_line(start+len(primer_fw)-1, end-len(primer_rv), y, color, name)
        self.show_primer(start, y, primer_fw, '#f00', primer_fw.name, forward=True)
        self.show_primer(end-len(primer_rv), y, primer_rv, '#00f', primer_rv.name, forward=False)

    def put_yline(self, x, y0, y1, color='#000'):
        i0 = self.get_inner_index(Point(x,y0))
        i1 = self.get_inner_index(Point(x,y1))
        c0 = self.get_mapped_index(i0)
        c1 = self.get_mapped_index(i1)
        self.renderer.draw_line(c0,c1,color=color)

    def put_line(self, x0, x1, y, color='#000', name=''):
        i0 = self.get_inner_index(Point(x0,y))
        i1 = self.get_inner_index(Point(x1,y))
        c0 = self.get_mapped_index(i0)
        c1 = self.get_mapped_index(i1)
        self.renderer.draw_textbar(c0.x, c1.x, c0.y, name, textstyle=self.textstyle, color=color)

    def put_bar(self, x0, x1, y, color='#000', name='', larrow=False, rarrow=False):
        i0 = self.get_inner_index(Point(x0,y))
        i1 = self.get_inner_index(Point(x1,y))
        c0 = self.get_mapped_index(i0)
        c1 = self.get_mapped_index(i1)
        self.renderer.draw_hbar(c0, style=-1 if larrow else 0, color=color)
        self.renderer.draw_hbar(c1, style=+1 if rarrow else 0, color=color)
        self.put_line(x0, x1, y, color, name)


class SeqCanvasFold(SeqCanvasBase):
    def __init__(self, start, end, vrow, mcolumn):
        '''
        vrow: number of characters of virtual rows
        mcolumn: number of characters of real column

        geometry:
          input index: x,y
                 start <= x < end, 0<=y<vrow
          innner index: ix, iy
                 0 <= ix < vcolumn, 0<=iy<vrow
          mapped index: mx, my
                 0 <= mx < mcolumn, 0 <= my < mrow
        '''
        super(SeqCanvasFold,self).__init__(start, end)
        self.vrow = vrow
        self.mcolumn = mcolumn
        self.mrow = (vrow) * ((self.vcolumn / mcolumn) + (1 if self.vcolumn%mcolumn>0 else 0))

    def get_mapped_index(self, i):
        return Point(i.x%self.mcolumn, (i.x/self.mcolumn * self.vrow + i.y))

    def put_line(self, x0, x1, y, color='#000', name=''):
        i0 = self.get_inner_index(Point(x0,y))
        i1 = self.get_inner_index(Point(x1,y))

        ''' bar and names '''
        namec = len(name)

        xx0 = i0.x%self.mcolumn
        yy0 = i0.x/self.mcolumn
        xx1 = i1.x%self.mcolumn
        yy1 = i1.x/self.mcolumn
        for yy in range(yy0,yy1+1):
            my = yy * self.vrow + i0.y
            s0 = Point(xx0 if yy==yy0 else 0, my)
            s1 = Point(xx1 if yy==yy1 else self.mcolumn-1, my)

            self.renderer.draw_textbar(s0.x, s1.x, my, name, textstyle=self.textstyle, color=color)

class SeqCanvasOverview(SeqCanvasBase):
    def __init__(self, start, end):
        super(SeqCanvasOverview,self).__init__(start,end)
        self.renderer = CanvasRendererScaleTransform(1, CanvasImage.font_size().y*2)
        self.size = self.renderer.transform(Point(end-start, 20))
        self.canvas = CanvasImage(self.size)
        self.renderer.set_target(self.canvas)
        self.textstyle=+1

class SeqCanvas(SeqCanvasFold):
    def __init__(self, start, end, vrow, mcolumn):
        super(SeqCanvas,self).__init__(start,end,vrow,mcolumn)
        fs = CanvasImage.font_size()
        self.renderer = CanvasRendererScaleTransform(fs.x, fs.y)
        self.size = self.renderer.transform(Point(self.mcolumn,self.mrow))
        self.canvas = CanvasImage(self.size)
        self.renderer.set_target(self.canvas)
        self.textstyle=0

