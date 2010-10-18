from xml.sax.saxutils import quoteattr,escape

__author__ = "MIZUGUCHi Yasuhiko"
__copyright__ = "Copyright (C) 2010 MIZUGUCHI Yasuhiko"
__license__ = "MIT"
__version__ = "1.0"

class HtmlWriter(object):
    def __init__(self,output,xhtml=True,indent=True):
        self._indent = indent
        self._output = output
        self._stack = []
        self._fpush = False
        if xhtml:
            self._ctag = self._ctag_xhtml
        else:
            self._ctag = self._ctag_html

    def _top(self):
        return self._stack[-1]
    def _write(self,text):
        self._output.write(text)
    def _nl(self):
        if not self._indent:
            return
        self._output.write('\n')

    def write(self,t):
        self._write(t)

    def writel(self,t):
        self._s()
        self._fpush = False
        self._write(t)
        self._nl()

    def text(self,t):
        self.write(escape(t))
    def textl(self,t):
        self.writel(escape(t))

    def attr_hook(self,attrv):
        if attrv=='cls':
            return 'class'
        return attrv

    def _starttag(self,tag,**attrs):
        return ' '.join([str(tag)]+['%s=%s'%(escape(self.attr_hook(str(k))),quoteattr(str(v))) for k,v in attrs.items() if k])

    def _stag(self,tag,**attrs):
        self._write('<%s>'%self._starttag(tag,**attrs))
    def _etag(self,tag):
        self._write('</%s>'%tag)
    def _ctag_xhtml(self,tag,**attrs):
        self._write('<%s/>'%self._starttag(tag,**attrs))
    def _ctag_html(self,tag,**attrs):
        self._write('<%s>'%self._starttag(tag,**attrs))

    def _s(self):
        if not self._indent:
            return
        if self._fpush:
            self._write('\n')
        self._write(' '*len(self._stack))

    def _e(self):
        if not self._indent:
            return
        if not self._fpush:
            self._write(' '*(len(self._stack)-1))

    def push(self,tag,**attrs):
        self._s()
        self._fpush = True
        self._stack.append(tag)
        self._stag(tag,**attrs)

    def pop(self):
        self._e()
        self._fpush = False
        tag = self._top()
        self._etag(tag)
        self._nl()
        self._stack = self._stack[:-1]

    def insertc(self,tag,**attrs):
        self._s()
        self._fpush = False
        self._ctag(tag,**attrs)
        self._nl()
        
    def insert(self,tag,text,**attrs):
        self._s()
        self._fpush = False
        self._stag(tag,**attrs)
        self._write(escape(text))
        self._etag(tag)
        self._nl()

class ListWriter:
    def __init__(self,htmlwriter):
        self.w = htmlwriter
        self.deep = 0
    def move(self,level):
        assert level > 0
        if self.deep == level:
            self.w.pop()
            self.w.push('li')
        elif self.deep < level:
            for i in range(level-self.deep):
                self.w.push('ul')
                self.w.push('li')
        else:
            for i in range(self.deep-level):
                self.w.pop()
                self.w.pop()
            self.w.pop()
            self.w.push('li')
        self.deep = level
    def finish(self):
        if self.deep:
            for i in range(self.deep):
                self.w.pop()
                self.w.pop()

class OUListWriter:
    def __init__(self,htmlwriter):
        self.w = htmlwriter
        self.deep = 0
        self.olstack = {}
        
    def move(self,level,ol=False):
        assert level > 0

        if level in self.olstack:
            samemark = not (self.olstack[level]^ol)
        else:
            samemark = True
        mark = {True:'ol',False:'ul'}[ol]
        self.olstack[level] = ol
            
        if self.deep == level:
            if samemark:
                self.w.pop()
                self.w.push('li')
            else:
                self.w.pop()
                self.w.pop()
                self.w.push(mark)
                self.w.push('li')
        elif self.deep < level:
            for i in range(level-self.deep):
                self.w.push(mark)
                self.w.push('li')
        else:
            for i in range(self.deep-level):
                self.w.pop()
                self.w.pop()
            if samemark:
                self.w.pop()
                self.w.push('li')
            else:
                self.w.pop()
                self.w.pop()
                self.w.push(mark)
                self.w.push('li')
        self.deep = level
    def finish(self):
        if self.deep:
            for i in range(self.deep):
                self.w.pop()
                self.w.pop()


class builder(object):
    def __init__(self, writer):
        self.w = writer
        self.queue = None

    def __getattr__(self, name):
        self._clear_queue()
        e = element(name, self)
        self.queue = e
        return e
    __getitem__ = __getattr__
    
    def _clear_queue(self):
        if self.queue:
            self.w.insertc(self.queue.name, **self.queue.attrs)
            self.queue = None

    def text(self, text):
        self.w.text(text)

    def get_writer(self):
        return self.w

class element(object):
    def __init__(self, name, builder):
        self.name = name
        self.builder = builder
        self.w = builder.w
        self.attrs = {}

    def __call__(self, value=None, **kargs):
        self.attrs.update(kargs)
        if value:
            with self:
                self.w.text(value)
        else:
            return self
        
    def __enter__(self):
        self.builder.queue = None
        self.w.push(self.name, **self.attrs)
    def __exit__(self, type, value, tb):
        self.builder._clear_queue()
        self.w.pop()
    
if __name__=='__main__':
    import sys
    w = WikiWriter(sys.stdout)
    w.push('body')
    w.push('body')
    w.push('a')
    w.pop()
    w.push('body')
    w.insertc('br')
    w.insert('title','title')
    w.pop()
    w.text('hoganona')
    w.text('hdkkaskdfn')
    w.text('hoganona3')
    w.a('hogehoge',href='http://&&Jlkjasdf.lkjalsdf==++||/lsls;&amp;',cls='tt',id='hogehoge',lalal='asdf')
    o = OUListWriter(w)
    o.move(1)
    w.text('a')
    o.move(2)
    w.text('bbbb1')
    w.text('bbbb2')
    w.text('bbbb3')
    o.move(1)
    w.text('a')
    o.move(1)
    w.text('a')
    o.move(1)
    o.move(2)
    w.text('bbbb1')
    w.text('bbbb2')
    w.text('bbbb3')
    w.text('a')
    o.move(1)
    w.text('a')
    o.finish()
    w.pop()
    w.pop()
    
    w = HtmlWriter(sys.stdout)
    b = builder(w)
    with b.html(hoge='hogehoge'):
        with b.head(div='div'):
            b.meta(name='meta name')
        b.hoge('anana',style='hoge')
        b.hoge('anananana')
        b.br
        b.hr
        b.hoge

