from __future__ import absolute_import

from Bio import Seq
from collections import defaultdict
from . import ColorMap, pprint_sequence_html, no_stop_in_frame
from ..util.memoize import memoize
from ..util import xmlwriter
from .cpg import bisulfite, cpg_sites, count_cpg

from .primer import PrimerPair

__all__ = ['PCR']

class PCRProduct(object):
    def __init__(self, i, j):
        assert(i.template == j.template)
        assert(i <= j)

        self.i = i
        self.j = j

        self.template = i.template

        self.start = i.left
        self.start_i = i.right
        self.end_i = j.left
        self.end = j.right

        self.fw = i.primer
        self.rv = j.primer

        self.head = self.fw.seq[:-self.i.length]
        self.fw_3 = self.fw.seq[-self.i.length:]
        self.middle = self.template[self.start_i:self.end_i]
        self.rv_3 = self.rv.seq.reverse_complement()[:self.j.length]
        self.tail = self.rv.seq.reverse_complement()[self.j.length:]

        self.seq = self.head + self.template[self.start:self.end] + self.tail

    def num_cpg(self):
        return len(self.cpg_sites())

    def cpg_sites(self):
        return cpg_sites(self.template, (self.start_i, self.end_i))

    def __repr__(self):
        return "PCRProduct(%s -> %s: %s)"%(self.fw.name, self.rv.name, self.seq)

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return str(self.seq)

    def write_html(self, w):
        b = xmlwriter.builder(w)

        with b.div:
            cm = ColorMap()

            parts = [(self.head, 200, 0, 0),
                     (self.fw_3, 0, 200, 0),
                     (self.middle, 0, 0, 0),
                     (self.rv_3, 0, 0, 200),
                     (self.tail, 200, 0, 0) ]

            allseq = Seq.Seq('')
            k = 0
            for seq, rr, gg, bb in parts:
                l = len(seq)
                for i in range(k, k+l):
                    cm.add_color(i, rr,gg,bb)
                k += l
                allseq += seq

            for i in cpg_sites(allseq):
                cm.add_color(i, 255, 0, 0)
                cm.add_color(i+1, 255, 0, 0)

            assert( len(allseq) == len(self.seq) )

            pprint_sequence_html(w, allseq, cm.get_color)
            with b.span(**{'class':'length'}):
                b.text('length=%d, CpG=%d, no stop=%s'%(len(allseq), count_cpg(allseq), no_stop_in_frame(allseq)))

            #with b.textarea(cols='10', rows='1', cls='copybox'):
            #    w.write(str(self.seq))



class PCRBand(object):
    def __init__(self, pcr):
        self.pcr = pcr
        self.bs_pos_met = None
        self.bs_pos_unmet = None
        self.bs_neg_met = None
        self.bs_neg_unmet = None
        self.origin = None

    def set_bs_product(self, methyl, sense, product):
        if sense:
            if methyl:
                self.bs_pos_met = product
            else:
                self.bs_pos_unmet = product
        else:
            if methyl:
                self.bs_neg_met = product
            else:
                self.bs_neg_unmet = product

    def get_bs_product(self, methyl, sense):
        if sense:
            return self.bs_pos_met if methyl else self.bs_pos_unmet
        else:
            return self.bs_neg_met if methyl else self.bs_neg_unmet

    def get_product(self):
        return self.origin or self.bs_pos_met or self.bs_pos_unmet or self.bs_neg_met or self.bs_neg_unmet

    def match_str(self):
        r = ''
        if self.bs_pos_met:
            r += "M"
        if self.bs_pos_unmet:
            r += "U"
        if self.bs_neg_met:
            r += "m"
        if self.bs_neg_unmet:
            r += "u"
        if self.origin:
            r += "G"
        return r

class PCR(object):
    def __init__(self, name, template, primer_fw, primer_rv):
        self.name = name
        self.template = template
        self.primers = PrimerPair(primer_fw, primer_rv)

    @property
    def fw(self):
        return self.primers.fw
    @property
    def rv(self):
        return self.primers.rv

    @property
    @memoize
    def pair_annealing(self):
        return self.primers.pair_annealing

    @property
    @memoize
    def pair_end_annealing(self):
        return self.primers.pair_end_annealing

    @property
    @memoize
    def primer_score(self):
        return self.primers.score

    @property
    @memoize
    def products(self):
        return list(self._search(self.template))

    def bs_products(self, methyl, sense):
        return list(self._search(bisulfite(self.template, methyl=methyl, sense=sense)))

    def _search(self, template):
        fpp,fpc = self.fw.search(template)
        rpp,rpc = self.rv.search(template)

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

    @property
    @memoize
    def bisulfite_products(self):
        """
        return list of the tuple: (met_product, unmet_product, genome_product)
        products in the same position are in the same tuple.
        if not all products are exist, the other values are None.
        """
        ret = defaultdict(lambda :PCRBand(self))
        
        for methyl in [True,False]:
            for sense in [True,False]:
                for p in self.bs_products(methyl, sense):
                    ret[(p.start,p.end,p.fw,p.rv)].set_bs_product(methyl,sense,p)

        for p in self.products:
            ret[(p.start,p.end,p.fw,p.rv)].origin = p

        return list(ret.values())

    def debugprint(self):
        print '%s: score=%.2f'%(self.name, self.primer_score())
        self.primers.debugprint()
        for c in self.products:
            print 'product: len=%d, detectable CpG=%d'%(len(c),c.detectable_cpg())
            print c.seq

    def write_html(self, w):
        b = xmlwriter.builder(w)
        with b.div:
            b.h2(self.name)

            self.primers.write_html(w)
            #self.pair_annealing().write_html(w)
            #self.pair_end_annealing().write_html(w)
            for c in self.products:
                c.write_html(w)
