

from Bio import Seq
from collections import defaultdict
from . import PPrintSequence, no_stop_in_frame
from ..util.memoize import memoize
from ..util import xmlwriter
from .cpg import bisulfite, bisulfite_conversion, cpg_sites, count_cpg, bisulfite_conversion_unambiguous

from .primer import PrimerPair, Primer

MAX_PRODUCT_SIZE = 1500

__all__ = ['PCR', 'BisulfitePCR', 'RtPCR']

class PCRProduct:
    fcount = 0
    
    def __init__(self, i, j):
        """

        template:
                     i.left  i.right       j.left j.right   
                         |     |             |    |         
                    |---------->                            
                         |------------------------|         
                                             <----------|
        seq:        0   v[0]  v[1]          v[2] v[3]  v[4]

        """
        assert(i.template == j.template)
        assert(i <= j)

        self.i = i
        self.j = j

        self.fw = i.primer
        self.rv = j.primer

        self.template = i.template

        self.start_0 =  i.left - i.adapter_length
        self.start = i.left
        self.start_i = i.right
        self.end_i = j.left
        self.end = j.right
        self.end_0 =  j.right + j.adapter_length

        self.head = self.fw.seq[:-self.i.length]
        self.fw_3 = self.fw.seq[-self.i.length:]
        self.middle = self.template[self.start_i:self.end_i]
        self.rv_3 = self.rv.seq.reverse_complement()[:self.j.length]
        self.tail = self.rv.seq.reverse_complement()[self.j.length:]

        if self.head:
            self.fwp = Primer('{}({}-mer)'.format(i.primer.name, len(self.fw_3)),self.fw_3)
        else:
            self.fwp = i.primer

        if self.tail:
            self.rvp = Primer('{}({}-mer)'.format(j.primer.name, len(self.rv_3)),self.rv_3.reverse_complement())
        else:
            self.rvp = j.primer

        self.partial_match = not not (self.head or self.tail)
        self.primerpair = PrimerPair(self.fwp, self.rvp)

        self.parts = [self.head, self.fw_3, self.middle, self.rv_3, self.tail]
        self.seq = self.head + self.fw_3 + self.middle + self.rv_3 + self.tail
        #self.seq = self.head + self.template[self.start:self.end] + self.tail
        v = []
        k = 0
        for part in self.parts:
            k += len(part)
            v.append(k)
        self.v = v
        self.partsv = [(0,v[0]), (v[0],v[1]), (v[1],v[2]), (v[2],v[3]), (v[3],v[4])]

    def cpg_sites(self):
        return cpg_sites(self.template, (self.start_i, self.end_i))

    def __repr__(self):
        return "PCRProduct(%s -> %s: %s)"%(self.fw.name, self.rv.name, self.seq)

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return str(self.seq)

    def write_html(self, b, toc, subfs, callback):
        w = b.get_writer()

        with b.div:
            if self.partial_match:
                self.primerpair.write_html(w, pair_values=False, annealings=False)

            seq = self.seq

            if len(self) < MAX_PRODUCT_SIZE:
                cm = PPrintSequence(seq)

                parts = [(100,100,100),
                         (255, 0, 0),
                         (0, 0, 0),
                         (0, 0, 255),
                         (100,100,100) ]

                for i, (rr, gg, bb) in enumerate(parts):
                    p,q = self.partsv[i]
                    for i in range(p,q):
                        cm.add_color(i, rr,gg,bb)

                for i in cpg_sites(seq, self.partsv[2]):
                    cm.add_underbar(i)
                    cm.add_underbar(i+1)

                cm.write_html(w)
                callback(self)
            else:
                w.text('sequence ommitted. product length = {}'.format(len(self)))

            with b.span(**{'class':'length'}):
                b.text('length=%d, CpG=%d'%(len(seq), count_cpg(self.middle)))

            filename = 'pcr_product_{}.seq.txt'.format(PCRProduct.fcount)
            PCRProduct.fcount += 1
            subfs.write(filename, str(self.seq))
            b.a('.seq file',href=subfs.get_link_path(filename))

class PCRProducts:
    def __init__(self, products):
        self.products = list(products)

    def __len__(self):
        return len(self.products)

    def write_html(self, b, toc, subfs, callback = None):
        w = b.get_writer()
        if len(self.products) > 0:
            with b.div(klass='products'):
                for c in self.products:
                    c.write_html(b, toc, subfs, callback)
        else:
            with b.div(klass='products'):
                b.p('no products')

    def __iter__(self):
        yield from self.products

class PCR:
    def __init__(self, name, template, primer_fw, primer_rv):
        """
        >>> fw = Primer('Fw','TTTCAGCAAGGACTGGTCTTT')
        >>> rv = Primer('Rv','CACAACTTTCAGCAGCTTACAAA')
        >>> t = Seq.Seq('TTTCAGCAAGGACTGGTCTTTCTATCTCTTGTACTACACTGAATTCACCCCCACTGAAAAAGATGAGTATGCCTGCCGTGTGAACCATGTGACTTTGTCACAGCCCAAGATAGTTAAGTGGGGTAAGTCTTACATTCTTTTGTAAGCTGCTGAAAGTTGTG')
        >>> pcr = PCR('B2M', t, fw, rv)
        """
        assert(isinstance(template, Seq.Seq))
        assert(isinstance(name, str))
        assert(isinstance(primer_fw, Primer))
        assert(isinstance(primer_rv, Primer))
        self.name = name
        self.template = template
        self.primers = PrimerPair(primer_fw, primer_rv)

        self.fw = self.primers.fw
        self.rv = self.primers.rv
        self.pair_annealing = self.primers.pair_annealing
        self.pair_end_annealing = self.primers.pair_end_annealing

        self.products = PCRProducts(self._search(self.template))

    @property
    @memoize
    def primer_score(self):
        return self.primers.score

    def _search(self, template, template_ambiguous=False):
        fpp,fpc = self.fw.search(template,template_ambiguous)
        rpp,rpc = self.rv.search(template,template_ambiguous)

        for i in fpp:
            for j in fpc: # fw -> fw
                if i <= j:
                    yield PCRProduct(i,j)
            for j in rpc: # fw -> rv
                if i <= j:
                    yield PCRProduct(i,j)
        for i in rpp:
            for j in fpc: # rv -> fw
                if i <= j:
                    yield PCRProduct(i,j)
            for j in rpc: # rv -> rv
                if i <= j:
                    yield PCRProduct(i,j)

    def debugprint(self):
        print('%s: score=%.2f'%(self.name, self.primer_score()))
        self.primers.debugprint()
        for c in self.products:
            print('product: len=%d, detectable CpG=%d'%(len(c),c.detectable_cpg()))
            print(c.seq)

    def write_html(self, b, toc, subfs, callback = None):
        self.primers.write_html(b.get_writer())
        self.products.write_html(b, toc, subfs, callback)

class PCRBand:
    """
    summary of PCR products which has exactly same location and different template.
    """
    def __init__(self):
        self.match_str = ""
        self.products = []

    def set_product(self, abbr, product):
        self.match_str += abbr
        self.products.append(product)

    def get_product(self):
        return self.products[0]

"""
list of tuples
(full-name, abbr-name, conversion_function)
"""
bisulfite_conversions = [
    ("Bisulfite-Treated (+)", "+", lambda x: bisulfite_conversion(x, sense=True)),
    ("Bisulfite-Treated (-)", "-", lambda x: bisulfite_conversion(x, sense=False)),
#    ("Bisulfite-Treated Sense Methyl", "M", lambda x: bisulfite(x, methyl=True, sense=True)),
#    ("Bisulfite-Treated Sense Unmethyl", "U", lambda x: bisulfite(x, methyl=False, sense=True)),
#    ("Bisulfite-Treated Antisense Methyl", "m", lambda x: bisulfite(x, methyl=True, sense=False)),
#    ("Bisulfite-Treated Antisense Unmethyl", "u", lambda x: bisulfite(x, methyl=False, sense=False)),
    ("Genome", "G", lambda x: x),
]

class PCRconv(PCR):
    """
    PCR for converted templates.
    """
    def __init__(self, name, original_template, primer_fw, primer_rv, conversions=[]):
        super().__init__(name, original_template, primer_fw, primer_rv)
        self.conversions = conversions
        self.converted_temp = []

        for name, abbr, conv in conversions:
            temp = conv(original_template)
            products = PCRProducts(self._search(temp,True))
            self.converted_temp.append((name, abbr, temp, products))

    def write_html(self, b, toc, subfs, callback = None):
        w = b.get_writer()
        self.primers.write_html(w)

        for name, abbr, temp, products in self.converted_temp:
            b.h5('template = {}'.format(name))
            products.write_html(b, toc, subfs, callback)

    def bands(self):
        ret = defaultdict(lambda :PCRBand())
        for name, abbr, temp, products in self.converted_temp:
            for p in products:
                ret[(p.start,p.end)].set_product(abbr,p)

        return list(ret.values())

class BisulfitePCR(PCRconv):
    def __init__(self, name, genome, primer_fw, primer_rv):
        super().__init__(name, genome, primer_fw, primer_rv, bisulfite_conversions)

    def get_bs_products(self):
        return self.converted_temp[0][3]

class PCRmulti(PCR):
    """
    PCR for multiple different templates
    """
    def __init__(self, name, templates, primer_fw, primer_rv):

        super().__init__(name, templates[0][1], primer_fw, primer_rv)
        self.templates = templates
        self.multi_temp = []

        for tname, seq in templates:
            pcr = PCR(name, seq, primer_fw, primer_rv)
            self.multi_temp.append((tname, pcr))

    def write_html(self, b, toc, subfs, callback = None):
        w = b.get_writer()
        self.primers.write_html(w)

        for tname, pcr in self.multi_temp:
            b.h5('template = {}'.format(tname))
            pcr.products.write_html(b, toc, subfs, callback)

class RtPCR(PCRmulti):
    def __init__(self, name, genome, transcripts, primer_fw, primer_rv):
        templates = [('genome',genome)] + [(t.name,t.seq) for t in transcripts]
        super().__init__(name, templates, primer_fw, primer_rv)

    def genome_pcr(self):
        return self.multi_temp[0][1]

    def transcript_pcr(self, transcript):
        return PCR(self.name, transcript.seq, self.fw, self.rv)
