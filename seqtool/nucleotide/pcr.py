from __future__ import absolute_import

import re

from collections import defaultdict
from . import to_seq, melt_temp, tm_gc, ColorMap, pprint_sequence_html
from ..util.memoize import memoize
from ..util import xmlwriter
from .. import parser
from .cpg import gc_ratio, bisulfite, cpg_sites, count_cpg

__all__ = ['Primer', 'PrimerPair', 'PrimerCondition', 'PCR', 'primers_write_html', 'primers_write_csv', 'load_primer_list_file']

def annealing_score_n(x,y):
    if (x=='A' and y=='T') or (x=='T' and y=='A'):
        return 2
    elif (x=='G' and y=='C') or (x=='C' and y=='G'):
        return 4
    else:
        return 0

def annealing_score(p,q,end_annealing=False,getindex=False):
    sv = annealing_score_n
    p = str(p).upper()
    q = str(q).upper()
    w = p
    v = q[::-1]
    n = len(w)
    m = len(v)

    def ea(ss):
        ret = 0
        for n in ss:
            if n==0:
                return ret
            ret += n
        return ret
    def ea_l(ss):
        return ea(ss)
    def ea_r(ss):
        return ea(ss[::-1])
    def ea_lr(ss):
        return max(ea_l(ss),ea_r(ss))

    def max_(old, new):
        old_s,v = old
        new_s,vv = new
        if old_s<new_s:
            return new
        return old

    eav = (-1, None)
    av = (-1, None)
    if n<=m:
        assert m-n >= 0
        for k in xrange(-(n-1),m-1 +1):
            if k<=0:
                # 5'- w[0]....w[-k]....w[n-1] -3'
                #         3'- v[0].....v[n+k-1]....v[m-1] -5'
                ss = [sv(w[-k+i],v[i]) for i in xrange(n+k)]
                av = max_(av, (sum(ss),(k,ss)))
                eav = max_(eav,(ea_lr(ss),(k,ss)))
            elif k<=m-n:
                #         w[0]....w[n-1]
                # v[0]....v[k]....v[k+n-1].....v[m-1]
                ss = [sv(w[0+i],v[k+i]) for i in xrange(n)]
                av = max_(av, (sum(ss),(k,ss)))
                eav = max_(eav,(ea_r(ss),(k,ss)))
            else:
                #        w[0]...w[m-k-1]....w[n-1]
                # v[0]...v[k]...v[m-1]
                ss = [sv(w[i],v[k+i]) for i in xrange(m-k)]
                av = max_(av, (sum(ss),(k,ss)))
    else:
        assert m-n <= 0
        for k in xrange(-(n-1),m-1 +1):
            if k<=m-n:
                # w[0]....w[-k]....w[n-1]
                #         v[0].....v[n+k-1]....v[m-1]
                ss = [sv(w[-k+i],v[i]) for i in xrange(n+k)]
                av = max_(av, (sum(ss),(k,ss)))
                eav = max_(eav,(ea_lr(ss),(k,ss)))
            elif k<=0:
                # w[0]....w[k]....w[m-k-1].....w[n-1]
                #         v[0]....v[m-1]
                ss = [sv(w[k+i],v[0+i]) for i in xrange(m)]
                av = max_(av, (sum(ss),(k,ss)))
                eav = max_(eav,(ea_l(ss),(k,ss)))
            else:
                #        w[0]...w[m-k-1]....w[n-1]
                # v[0]...v[k]...v[m-1]
                ss = [sv(w[i],v[k+i]) for i in xrange(m-k)]
                av = max_(av, (sum(ss),(k,ss)))

    if not end_annealing:
        return av[0], av[1]
    else:
        return eav[0], eav[1]


class Primer(object):
    def __init__(self, name, seq):
        self.name = name
        self.seq = to_seq(seq)

    def __repr__(self):
        return "Primer(%s: %s)"%(self.name,self.seq)
    def __len__(self):
        return len(self.seq)
    def __str__(self):
        return str(self.seq)

    def reverse(self):
        return self.seq[::-1]

    @property
    @memoize
    def gc_ratio(self):
        return gc_ratio(self.seq)

    def melting_temperature(self, pcr_mix=melt_temp.DEFAULT_MIX, unmethyl=True):
        return melt_temp.melting_temperature_unmethyl(self.seq, pcr_mix, unmethyl)

    @property
    @memoize
    def self_annealing(self):
        return PrimerAnnealing(self.seq, self.seq)
    @property
    @memoize
    def self_end_annealing(self):
        return PrimerAnnealing(self.seq, self.seq, True)

    sa = self_annealing
    sea = self_end_annealing

    def search(self, template, match_length=None):
        """
        Return tuple of fowards anneal locations and reverse anneal locations.
        anneal locations are (5'location, 3'location)
        """
        def _gen_re_seq(seq):
            table = {
                'A': 'A',
                'T': 'T',
                'G': 'G',
                'C': 'C',
                'R': '[GA]',
                'Y': '[TC]',
                'M': '[AC]',
                'K': '[GT]',
                'S': '[GC]',
                'W': '[AT]',
                'H': '[ACT]',
                'B': '[GTC]',
                'V': '[GCA]',
                'D': '[GAT]',
                'N': '[GATC]',
                }
            return ''.join([table[s] for s in str(seq).upper()])

        if match_length and len(self.seq) > match_length > 0:
            l = len(self.seq) - match_length
            seq = self.seq[l:]
        else:
            seq = self.seq
        primer = seq
        cprimer = seq.reverse_complement()
        reg = re.compile('(%s)|(%s)'%(_gen_re_seq(primer),_gen_re_seq(cprimer)))

        template = str(template).upper()
        pp = []
        pc = []
        start = 0
        while True:
            m = reg.search(template, start)
            if not m:
                break
            start = m.start()+1
            length = m.end()-m.start()
            if m.group(1):
                pp.append(PrimerTemplateAnnealing(self, template, True, m.end()-1, length))
            if m.group(2):
                pc.append(PrimerTemplateAnnealing(self, template, False, m.start(), length))
        return pp,pc

    def debugprint(self):
        print "%s: %s, len=%2d, Tm=%.2f, oTm=%.2f, MarmurTm: %s, GC=%.2f, sa=%2d, sea=%2d" % self.get_table_row()
        print 'sa=%s, index=%s'%(self.sa.score, self.sa.index)
        print '\n'.join(self.sa.get_bar())
        print 'sea=%s, index=%s'%(self.sea.score, self.sea.index)
        print '\n'.join(self.sea.get_bar())


    @classmethod
    def get_table_head(cls):
        return ['name', 'sequence', 'length[bp]', 'SSDS[bp]', 'Tm[C]', 'oTm[C]', 'old Tm[C]','GC[%]', 'sa.', 'sea.']

    def get_table_row(self):
        return [str(v) for v in [self.name,
                                 "5'-%s-3'"%self.seq,
                                 len(self),
                                 self.sdss_length(),
                                 '%.2f'%self.melting_temperature(unmethyl=True),
                                 '%.2f'%self.melting_temperature(unmethyl=False),
                                 tm_gc(self.seq),
                                 '%.2f'%self.gc_ratio,
                                 self.sa.score,
                                 self.sea.score ]]
    def _score(self, bisulfite=False):
        pc = PrimerCondition(bisulfite)

        return pc.primer_length.score(len(self)) \
              + pc.gc.score(self.gc_ratio) \
              + pc.tm.score(self.melting_temperature()) \
              + pc.sa.score(self.sa.score) \
              + pc.sea.score(self.sea.score)

    @property
    @memoize
    def score(self):
        return self._score(False)

    @property
    @memoize
    def score_bisulfite(self):
        return self._score(True)


    def sdss_length(self, anneal_temp=None, pcr_mix=melt_temp.DEFAULT_MIX, unmethyl=False):
        s = str(self.seq).upper()
        if unmethyl:
            s = s.replace('R','A').replace('Y','T')
        else:
            s = s.replace('R','G').replace('Y','C')

        if not anneal_temp:
            anneal_temp = self.melting_temperature()
            
        l = len(s)
        for i in xrange(2,l):
            f = melt_temp.complex_fraction(s[l-i:-1], anneal_temp, pcr_mix)
            if f>0.01:
                break
        return i


class Condition(object):
    def __init__(self, weight, optimal, minimum, maximum):
        self.optimal = optimal
        self.weight = weight
        self.minimum = minimum
        self.maximum = maximum
    def score(self, value):
        return self.weight * abs(value-self.optimal)
    def bound(self, value):
        return self.minimum <= value <= self.maximum

class PrimerCondition(object):
    def __init__(self, bisulfite=False):
        self.primer_length = Condition(0.5, 23., 10, 30)
        if not bisulfite:
            self.gc        = Condition(1.0, 50., 30, 70)
        else:
            self.gc        = Condition(1.0, 30., 0, 60)
        self.tm            = Condition(1.0, 60., 45, 70)
        self.sa            = Condition(0.1, 0., 0, 25)
        self.sea           = Condition(0.2, 0., 0, 15)
        self.pa            = Condition(0.1, 0., 0, 25)
        self.pea           = Condition(0.2, 0., 0, 15)

class PrimerPair(object):
    def __init__(self, fw, rv):
        self.fw = fw
        self.rv = rv

    @property
    @memoize
    def pair_annealing(self):
        return PrimerAnnealing(self.fw.seq, self.rv.seq)
    @property
    @memoize
    def pair_end_annealing(self=False):
        return PrimerAnnealing(self.fw.seq, self.rv.seq, True)

    pa = pair_annealing
    pea = pair_end_annealing

    def _score(self, bisulfite):
        pc = PrimerCondition(bisulfite)
        return self.fw._score(bisulfite) + self.rv._score(bisulfite) + pc.pa.score(self.pa.score) + pc.pea.score(self.pea.score)

    @property
    @memoize
    def score(self):
        return self._score(False)

    @property
    @memoize
    def score_bisulfite(self):
        return self._score(True)

    def debugprint(self):
        print 'score=', self.score
        print 'bisulfite score=', self.score_bisulfite
        for r in [self.fw, self.rv]:
            r.debugprint()
        pa = self.pair_annealing()
        print 'pa=%s, index=%s'%(pa.score,pa.index)
        print '\n'.join(pa.get_bar())
        pea = self.pair_end_annealing()
        print 'pea=%s, index=%s'%(pea.score,pea.index)
        print '\n'.join(pea.get_bar())


    def write_html(self, w):
        b = xmlwriter.builder(w)
        # primer table
        with b.div(cls='primerpair'):
            with b.table(cls='primerpairtable', border=1):
                with b.tr:
                    for p in (Primer.get_table_head() + ['pa.', 'pea.', 'score', 'bisulfite score']):
                        b.th(p)

                with b.tr:
                    for v in self.fw.get_table_row():
                        b.td(str(v))
                    for v in [self.pa.score, self.pea.score, '%.2f'%self.score, '%.2f'%self.score_bisulfite]:
                        b.td(str(v),rowspan='2')
                with b.tr:
                    for v in self.rv.get_table_row():
                        b.td(str(v))

            b.p('Tm melting temperature(SantaLucia), oTm melting temperature for bisulfite methyl-template, sa. self annealing, sea. self end annealing, pa. pair annealing, pea. pair end annealing', style='font-size:x-small')

def primers_write_csv(primers):
    ret = ""
    ret += ", ".join(Primer.get_table_head())
    ret += "\n"

    for r in primers:
        ret += (", ".join([str(v) for v in r.get_table_row()]) + "\n")

    return ret+"\n"

def primers_write_html(w, primers):
    b = xmlwriter.builder(w)
    # primer table
    with b.div(cls='primerpair'):
        with b.table(cls='primerpairtable', border=1):
            with b.tr:
                for p in Primer.get_table_head():
                    b.th(p)

            for r in primers:
                with b.tr:
                    for v in r.get_table_row():
                        b.td(str(v))

        b.p('Tm melting temperature(SantaLucia), oTm melting temperature for bisulfite methyl-template, sa. self annealing, sea. self end annealing, pa. pair annealing, pea. pair end annealing', style='font-size:x-small')

class PrimerAnnealing(object):
    def __init__(self, p, q, end_annealing=False):
        s, (i, ss) = annealing_score(p,q,end_annealing)
        self.p = p
        self.q = q
        self.score = s
        self.scores = ss
        self.index = i

    def get_bar(self):
        i = self.index
        p = self.p
        q = self.q
        spc = ' '*abs(i)
        ss = self.scores
        bar = ''.join(['|' if s>0 else ' ' for s in ss])

        if i>0:
            return [spc+"5'-%s-3'"%p, spc+"  <"+bar+">", "3'-%s-5'"%q[::-1] ]
        else:
            return ["5'-%s-3'"%p, spc+"  <"+bar+">", spc+"3'-%s-5'"%q[::-1] ]

    def write_html(self, w):
        w.push('div',style='annealing')
        w.push('p','score=%s, index=%s'%(self.score,self.index))
        w.push('pre')
        w.text('\n'.join(self.get_bar()))
        w.pop()
        w.pop()


class PrimerTemplateAnnealing(object):
    def __init__(self, primer, template, strand, loc_3p, length):
        """
        strand == True

            5      3
            ------->
            |       |
           left    right

        strand == False

            3      5
            <-------
            |       |
           left    right

        """
        self.primer = primer
        self.template = template
        self.strand = strand
        self.length = length
        self.loc_3p = loc_3p

        if self.strand:
            self.loc_5p = self.loc_3p - self.length + 1
            self.left =  self.loc_5p
            self.right = self.loc_3p + 1
        else:
            self.loc_5p = self.loc_3p + self.length - 1
            self.left =  self.loc_3p
            self.right = self.loc_5p + 1

    def __le__(self, rhs):
        return (self.left <= rhs.left) and (self.right <= rhs.right)


class PCRProduct(object):
    def __init__(self, i, j):
        assert(i.template == j.template)
        assert(i <= j)

        self.template = i.template

        self.start = i.left
        self.start_i = i.right
        self.end_i = j.left
        self.end = j.right

        self.fw = i.primer
        self.rv = j.primer

        self.head = self.fw.seq
        self.middle = self.template[self.start_i:self.end_i]
        self.tail = self.rv.seq.reverse_complement()

        self.seq = self.head + self.middle + self.tail

    def __repr__(self):
        return "PCRProduct(%s -> %s: %s)"%(self.fw.name, self.rv.name, self.seq)

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return str(self.seq)

    def detectable_cpg(self):
        return len(self.cpg_sites())

    def cpg_sites(self):
        s = self.start_i-1
        e = self.end_i+1
        return [i+s for i in cpg_sites(self.template[s:e])]

    def write_html(self, w):
        b = xmlwriter.builder(w)

        with b.div:
            s = self.start
            cm = ColorMap()
            for i in range(0, self.start_i-s):
                cm.add_color(i, 0, 200, 0)
            for i in range(self.end_i-s, self.end-s):
                cm.add_color(i, 0, 0, 200)
            for i in self.cpg_sites():
                cm.add_color(i-s, 255, 0, 0)
                cm.add_color(i+1-s, 255, 0, 0)

            pprint_sequence_html(w, self.seq, cm.get_color)
            with b.span(**{'class':'length'}):
                b.text('length=%d, CpG=%d, CpG between primers=%d'%(len(self.seq), count_cpg(self.seq), self.detectable_cpg()))

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

def load_primer_list_file(fileobj):
    primer = []
    for name, value, em in parser.parse_file(fileobj):
        primer.append(Primer(name, value))
    return primer
