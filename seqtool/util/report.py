from ..util import xmlwriter
from ..util.subfs import SubFileSystem
from ..util.dirutils import Filepath
import io
import contextlib

REPORT_CSS = '''
body{font-family: monospace}
.images{}
.image{border: solid 1px;}
.template{margin-left: 1em;}
.indent{margin-left: 3em;}
.pcr{margin: 1em; padding: 1em;}
.products{margin-left: 2em;}
.length{margin-left: 5em;}
.copybox{margin-left:4em;}
.section{margin: 1em; padding: 1em;}
.primerpairtable{ font-family: monospace }
'''

@contextlib.contextmanager
def section(b, tb, title, klass='section'):
    anchor = 'TODO'
    with tb.section(title, anchor):
        b.h3(title, anchor=anchor)
        with b.div(klass=klass):
            yield


def write_html(outputp, title, html_content):
    """
    def html_content(self, b, toc, subfs)
    """
    assert(isinstance(outputp, Filepath))
    subfs = SubFileSystem(outputp.dir, outputp.prefix)

    out = xmlwriter.switchable_output()
    html = xmlwriter.XmlWriter(out)
    toc_out = io.StringIO()
    b = xmlwriter.builder(html)
    tb = xmlwriter.toc_builder(title)
    with b.html:
        with b.head:
            with b.style(type='text/css'):
                b.text(REPORT_CSS)
            b.title(title)
    with b.body:
        # table of contents
        b.h2('Table Of Contents')
        with b.div(klass='toc'):
            with b.ul:
                out.insert(toc_out)
        with b.div(klass='main'):
            html_content(b, tb, subfs)

    subfs.finish()

    tb.write(xmlwriter.XmlWriter(toc_out))

    with open(outputp.path,'w') as output:
        output.write(out.getvalue())

'''
                count = 0
                for child in self.html_items():
                    count += 1
                    name = child.title or '%s'%count
                    subsubfs = subfs.get_subfs(name)
                    child.html_section(b, tb, subsubfs)
'''
