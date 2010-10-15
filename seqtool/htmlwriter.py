from xml.sax.saxutils import quoteattr,escape

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
        return ' '.join([str(tag)]+['%s=%s'%(escape(self.attr_hook(str(k))),quoteattr(str(v))) for k,v in attrs.items() if k and v])

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


class WikiWriter(HtmlWriter):
    def __init__(self,output,xhtml=True):
        super(WikiWriter,self).__init__(output,xhtml)
        self.slink = 0

    def _stag(self,tag,**attrs):
        if tag.strip().lower()=='a':
            self.slink+=1
            if self.slink>1:
                return
        super(WikiWriter,self)._stag(tag,**attrs)

    def _etag(self,tag):
        if tag.strip().lower()=='a':
            self.slink -= 1
            if self.slink>0:
                return
        super(WikiWriter,self)._etag(tag)

    # utility functions
    def hr(self):
        self.insertc('hr')
    def br(self):
        self.insertc('br')
        
    def hr_full(self):
        self.insertc('hr',cls='full')

    def meta(self,name,content):
        self.insertc('meta',name=name,content=content)

    def a(self,text,href,**attrs):
        attrs.update(href=href)
        self.insert('a',text,**attrs)

    def aimg(self,url,img,title,width=None,height=None):
        self.push('a',href=img)
        self.img(img,title,width,height)
        self.pop()

    def _style(self,width=None,height=None):
        r = ''
        if width:
            r += 'width:%s;'%width
        if height:
            r += 'height:%s;'%height
        return r;

    def img(self,img,title,width=None,height=None):
        self.insertc('img',src=img,alt=title,title=title,style=self._style(width,height))

    def space(self):
        self.writel('&nbsp;')

    def empathis(self,text):
        self.insert('span',text,cls='empathis')
    def underline(self,text):
        self.insert('span',text,cls='underline')
    def deleteline(self,text):
        self.insert('span',text,cls='deleteline')

    def link_edit(self,link_command):
        self.a('[EDIT]',href=link_command,cls='link_edit')
        # link_command are not to be escaped.
        # self.insert('a','[EDIT]','href="%s"'%link_command,cls='link_edit')
        
    def link_wiki(self,text,href):
        self.a(text,href=href,cls='wiki')
    def link_wikinotfound(self,text,href):
        self.a(text,href=href,cls='wikinotfound')
    def link_footnote(self,text,href,iden=None):
        self.a(text,href=href,cls='footnote',iden=iden)
    def link_external(self,text,href):
        self.a(text,href=href,cls='external')
    def anchor(self,name,aname):
        self.a(name,href='#'+aname,id=aname)

    def errormsg(self,msg):
        self.br()
        self.insert('span','ERROR: '+msg,cls='error')
        self.br()

    def insert_comment_box(self,name,value=""):
        self.insertc('input',type='text', name=name, size="70", value=value)

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
    
