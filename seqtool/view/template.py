

from Bio import SeqIO
from io import StringIO
from ..util.memoize import memoize
from ..nucleotide.cpg import bisulfite

class BaseTemplate(object):
    def __init__(self):
        pass

    def __len__(self):
        return len(self.seq)

    @property
    def locus(self):
        return None
    
    @property
    def name(self):
        return ""

    @property
    def description(self):
        return ""

    @property
    def seq(self):
        raise NotImplementedError();

    @property
    @memoize
    def bs_met(self):
        return bisulfite(self.seq, True)

    @property
    @memoize
    def bs_unmet(self):
        return bisulfite(self.seq, False)

    @property
    def transcripts(self):
        return []

    @property
    @memoize
    def transcript_start_site(self):
        return 0
        
    @property
    def features(self):
        return []

class SequenceTemplate(BaseTemplate):
    def __init__(self, sequence):
        self.seq_ = sequence

    @property
    def seq(self):
        return self.seq_

class Transcript(object):
    def __init__(self, name, seq, feature):
        self.name = name
        self.seq = seq
        self.feature = feature

class GenbankTemplate(BaseTemplate):
    def __init__(self, genbank_content, locus=None):
        #with prompt('loading genbank...', '...done.'):
        self.genbank_ = SeqIO.read(StringIO(genbank_content), "genbank")
        self.locus_ = locus

    @property
    def locus(self):
        return self.locus_

    @property
    def genbank(self):
        return self.genbank_

    @property
    def name(self):
        return self.name_

    @property
    def description(self):
        return self.genbank_.description

    @property
    def seq(self):
        return self.genbank_.seq

    @property
    def transcripts(self):
        for f in self.features:
            if f.type=='mRNA':
                name = f.qualifiers['product'][0]
                yield Transcript(name, f.extract(self.genbank_.seq), f)
    
    @property
    @memoize
    def transcript_start_site(self):
        v = [int(t.location.start) for t in self.features if t.type=='mRNA']
        return min(v) if v else 0
        
    @property
    def features(self):
        return self.genbank_.features
