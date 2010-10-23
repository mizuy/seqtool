# encoding:utf-8
from seqtool import xmlwriter
from nose.tools import *

from StringIO import StringIO

def test_xmlwriter():
    s = StringIO()
    w = xmlwriter.XmlWriter(s)
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
    w.insert('a','hogehoge',href='http://&&Jlkjasdf.lkjalsdf==++||/lsls;&amp;',cls='tt',id='hogehoge',lalal='asdf')

    w.pop()
    w.pop()

    eq_(s.getvalue(), \
"""<body>
 <body>
  <a></a>
  <body>
   <br/>
   <title>title</title>
  </body>
hoganonahdkkaskdfnhoganona3  <a lalal="asdf" href="http://&amp;&amp;Jlkjasdf.lkjalsdf==++||/lsls;&amp;amp;" id="hogehoge" class="tt">hogehoge</a>
 </body>
</body>
""")

def test_listwriter():
    s = StringIO()
    w = xmlwriter.XmlWriter(s)

    w.push('html')
    w.push('body')

    o = xmlwriter.OUListWriter(w)
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
    eq_(s.getvalue(), \
"""<html>
 <body>
  <ul>
   <li>a
    <ul>
     <li>bbbb1bbbb2bbbb3</li>
    </ul>
   </li>
   <li>a</li>
   <li>a</li>
   <li>
    <ul>
     <li>bbbb1bbbb2bbbb3a</li>
    </ul>
   </li>
   <li>a</li>
  </ul>
 </body>
</html>
""")

def test_unicode():
    "note: do not use cStringIO.StringIO, which can not handle unicode"
    s = StringIO()
    w = xmlwriter.XmlWriter(s)

    w.push('body')
    w.insert('title',u'タイトル')
    w.pop()

    eq_(s.getvalue(), \
u"""<body>
 <title>タイトル</title>
</body>
""")
    

def test_builder():
    buff = StringIO()
    w = xmlwriter.XmlWriter(buff)
    b = xmlwriter.builder(w)
    with b.html(hoge='hogehoge'):
        with b.head(div='div'):
            b.meta(name='meta name')
        b.hoge('anana',style='hoge')
        b.hoge('anananana')
        b.br
        b.hr
        b.hoge
        b.foo(style='fill:color;')
        with b.bar(**{'xmlns:my':'http://.....'}):
            b.text('something')
            with b.div:
                b.text('something')
                
    eq_(buff.getvalue(), \
'''<html hoge="hogehoge">
 <head div="div">
  <meta name="meta name"/>
 </head>
 <hoge style="hoge">anana</hoge>
 <hoge>anananana</hoge>
 <br/>
 <hr/>
 <hoge/>
 <foo style="fill:color;"/>
 <bar xmlns:my="http://.....">something
  <div>something</div>
 </bar>
</html>
''')
