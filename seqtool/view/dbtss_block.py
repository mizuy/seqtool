from . import block
from ..db.dbtss import TissuesetLocus

class DbtssBlock(block.BaseBlock):
    def __init__(self, template):
        self.template = template
        self.tss = None
        self.tsl = None

        self.tsss = []
        self.tss_name_counter = 1
    
    def set_tissueset(self, tissues):
        if not self.template.locus:
            print("No Locus Defined.")
            return False
        self.tsl = TissuesetLocus(tissues, self.template.locus)
        return True

    def add_default_tss(self, tss_name="Assumed TSS"):
        p = self.template.transcript_start_site
        start,end = p-200, p+200
        self.tsss.append((tss_name, start, end))

    def add_tss(self, start, end, name=None):
        assert(not not start and not not end)
        if not name:
            name = "%s TSS No. %s" % (self.name, self.tss_name_counter)
            self.tss_name_counter += 1
        self.tsss.append((name, start, end))
    
    def tss_count_csv(self):
        content = ''
        for name, start, stop in self.tsss:
            content += ','.join([name] + self.tsl.count_tags(start,stop)) + '\n'
        return content

    def svg_genome(self, t):
        assert(self.template)
        length = len(self.template)

        if self.tsl:
            for r in self.tsl:
                t.add_hline(length)
                t.add_dbtss_track(r, self.tsl.maxtag, self.template.seq)
            t.add_hline(length)

    def write_tss_count_csv(self, subfs):
        content = ''
        content += ', '.join(['range \\ tissue']+[t.name for t in self.tss])
        content += '\n'
        for name,start,end in self.tss_count:
            content += ', '.join([name]+[str(t.count_range(start,end)) for t in self.tss])
            content += '\n'
        subfs.write('tss.csv', content)
