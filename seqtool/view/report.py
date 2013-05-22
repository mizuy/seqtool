

from ..util import xmlwriter

from ..util.namedlist import DefaultNamedList

from ..util.subfs import SubFileSystem
from ..util.dirutils import Filepath

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

.primerpairtable{ font-family: monospace }
'''

class ReportEntityBase(object):
    def anchor(self):
        return "#{}".format(self.title.replace(' ','_'))

    def title(self):
        pass

    def write_html(self, b, subfs):
        pass

class Report(object):
    def __init__(self):
        self.entities = []

    def append(self, title, entity):
        self.entities.append(title,entity)

    def write_html(self, outputp):
        subfs = SubFileSystem(outputp.dir, outputp.prefix)

        with open(outputp.path,'w') as output:
            html = xmlwriter.XmlWriter(output)
            b = xmlwriter.builder(html)
            with b.html:
                with b.head:
                    with b.style(type='text/css'):
                        b.text(REPORT_CSS)
            with b.body:
                count = 0
                with b.ul:
                    for gt in self.entities:
                        with b.li:
                            b.a(href=gt.anchor, gt.title)

                for gt in self.entities:
                    count += 1
                    name = gt.name or '%s'%count
                    subsubfs = subfs.get_subfs(name)
                    gt.write_html(b, subsubfs)

        subfs.finish()

    def load_seqv(self, filename):
        e = SeqviewEntity.load_seqv(filename)
        self.append(e)
        return e

    def load_gene(self,name, gene_id):
        e = SeqviewEntity.create_gene(name, gene_id)
        self.append(e)
        return e
