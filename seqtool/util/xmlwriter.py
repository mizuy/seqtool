from xml.sax.saxutils import quoteattr,escape

__author__ = "MIZUGUCHi Yasuhiko"
__copyright__ = "Copyright (C) 2010 MIZUGUCHI Yasuhiko"
__license__ = "MIT"
__version__ = "1.0"

class XmlWriter(object):
    def __init__(self, output, indent=True):
        self._indent = indent
        self._output = output
        self._stack = []
        self._fpush = False

    def _top(self):
        return self._stack[-1]
    def _write(self,text):
        self._output.write(text)
    def _nl(self):
        if not self._indent:
            return
        self._output.write('\n')

    def write(self, t):
        self._write(t)

    def writel(self, t):
        self._s()
        self._fpush = False
        self._write(t)
        self._nl()

    def text(self, t):
        self.write(escape(str(t)))
        
    def textl(self, t):
        self.writel(escape(str(t)))

    def attr_hook(self, attrv):
        if attrv=='cls':
            return 'class'
        return attrv

    def _starttag(self, tag, **attrs):
        return ' '.join([str(tag)]+['%s=%s'%(escape(self.attr_hook(str(k))),quoteattr(str(v))) for k,v in attrs.items() if k])

    def _stag(self, tag, **attrs):
        self._write('<%s>'%self._starttag(tag,**attrs))
    def _etag(self, tag):
        self._write('</%s>'%tag)
    def _ctag(self, tag, **attrs):
        self._write('<%s/>'%self._starttag(tag,**attrs))

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

    def push(self, tag, **attrs):
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

    def insertc(self, tag, **attrs):
        self._s()
        self._fpush = False
        self._ctag(tag,**attrs)
        self._nl()
        
    def insert(self, tag, text, **attrs):
        self._s()
        self._fpush = False
        self._stag(tag,**attrs)
        self._write(escape(str(text)))
        self._etag(tag)
        self._nl()

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
        self.w.text(str(text))

    def write_raw(self, text):
        self.w.write(text)

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

class ListWriter:
    def __init__(self, writer):
        self.w = writer
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
    def __init__(self, writer):
        self.w = writer
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

