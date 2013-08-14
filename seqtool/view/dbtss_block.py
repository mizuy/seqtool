from . import block
from .. import db

class DbtssBlock(block.BaseBlock):
    def __init__(self, template, tissues):
        self.template = template
        self.gsl = db.dbtss.get_genomeset_locus(tissues, self.template.locus)
    
    def svg_genome(self, t):
        assert(self.template)
        length = len(self.template)

        if self.gsl:
            for r in self.gsl:
                t.add_hline(length)
                t.add_dbtss_track(r, self.gsl.maxtag, self.template.seq)
            t.add_hline(length)
