from xml.sax.saxutils import quoteattr,escape
import io

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
        if attrv=='cls' or attrv=='klass':
            return 'class'
        return attrv

    def _starttag(self, tag, **attrs):
        return ' '.join([str(tag)]+['%s=%s'%(escape(self.attr_hook(str(k))),quoteattr(str(v))) for k,v in list(attrs.items()) if k])

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

    def push(self):
        self.move(self.deep+1)
    def pop(self):
        self.move(self.deep-1)

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

class toc_builder(object):
    """
    toc = toc_builder('toc','toc')
    with toc.section('chapter 1'):
        with toc.section('section 1'):
            brabra
        with toc.section('section 4'):
            bra
            with toc.section('subsection 1'):
                pass
    toc.write(w)

    <ul>
        <li>chapter 1</li>
        <li>chapter 1</li>
        <ul>
            <li>chapter 1</li>
        </ul>
    </ul>
    """
    def __init__(self, title, anchor=None):
        self.top = toc_element(title, anchor, self)
        self.stack = [self.top]

    def write(self, w):
        w.push('ul')
        self.top.html(w)
        w.pop()

    def section(self, title, anchor=None):
        e = toc_element(title, anchor, self)
        self.stack[-1].add(e)
        return e

    def push(self, e):
        self.stack.append(e)

    def pop(self):
        self.stack = self.stack[:-1]
    
class toc_element(object):
    def __init__(self, name, anchor, parent):
        self.name = name
        self.anchor = anchor or name
        self.parent = parent
        self.children = []

    def add(self, child):
        self.children.append(child)

    def __enter__(self):
        self.parent.push(self)

    def __exit__(self, type, value, tb):
        self.parent.pop()

    def html(self, w):
        w.push('li')
        w.insert('a', self.name, href='#'+self.anchor)
        w.pop()
        if self.children:
            w.push('ul')
            for child in self.children:
                child.html(w)
            w.pop()

class switchable_output(object):
    def __init__(self):
        self._l = [io.StringIO()]

    def insert(self, output):
        self._l.append(output)
        self._l.append(io.StringIO())

    def write(self, text):
        self._l[-1].write(text)

    def getvalue(self):
        return ''.join(v.getvalue() for v in self._l)

