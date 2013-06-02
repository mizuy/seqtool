import re

from . import to_seq, melt_temp, tm_gc
from ..util.memoize import memoize
from ..util import xmlwriter
from ..util.parser import TreekvParser
from .cpg import gc_ratio
from ..util.namedlist import NamedList
from . import iupac
from . import primer_cond as pcond
from . import melt_temp

__all__ = ['Primer', 'PrimerPair']

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
        for k in range(-(n-1),m-1 +1):
            if k<=0:
                # 5'- w[0]....w[-k]....w[n-1] -3'
                #         3'- v[0].....v[n+k-1]....v[m-1] -5'
                ss = [sv(w[-k+i],v[i]) for i in range(n+k)]
                av = max_(av, (sum(ss),(k,ss)))
                eav = max_(eav,(ea_lr(ss),(k,ss)))
            elif k<=m-n:
                #         w[0]....w[n-1]
                # v[0]....v[k]....v[k+n-1].....v[m-1]
                ss = [sv(w[0+i],v[k+i]) for i in range(n)]
                av = max_(av, (sum(ss),(k,ss)))
                eav = max_(eav,(ea_r(ss),(k,ss)))
            else:
                #        w[0]...w[m-k-1]....w[n-1]
                # v[0]...v[k]...v[m-1]
                ss = [sv(w[i],v[k+i]) for i in range(m-k)]
                av = max_(av, (sum(ss),(k,ss)))
    else:
        assert m-n <= 0
        for k in range(-(n-1),m-1 +1):
            if k<=m-n:
                # w[0]....w[-k]....w[n-1]
                #         v[0].....v[n+k-1]....v[m-1]
                ss = [sv(w[-k+i],v[i]) for i in range(n+k)]
                av = max_(av, (sum(ss),(k,ss)))
                eav = max_(eav,(ea_lr(ss),(k,ss)))
            elif k<=0:
                # w[0]....w[k]....w[m-k-1].....w[n-1]
                #         v[0]....v[m-1]
                ss = [sv(w[k+i],v[0+i]) for i in range(m)]
                av = max_(av, (sum(ss),(k,ss)))
                eav = max_(eav,(ea_l(ss),(k,ss)))
            else:
                #        w[0]...w[m-k-1]....w[n-1]
                # v[0]...v[k]...v[m-1]
                ss = [sv(w[i],v[k+i]) for i in range(m-k)]
                av = max_(av, (sum(ss),(k,ss)))

    if not end_annealing:
        return av[0], av[1]
    else:
        return eav[0], eav[1]


class PrimerAnnealing:
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
        w.insert('p','score=%s, index=%s'%(self.score,self.index))
        w.push('pre')
        w.text('\n'.join(self.get_bar()))
        w.pop()
        w.pop()

def count_while(iteration):
    count = 0

    for i in iteration:
        if i:
            count += 1
        else:
            break

    return count

class PrimerTemplateAnnealing:
    def __init__(self, primer, template, strand, loc_3p):
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

        self.loc_3p = loc_3p

        if self.strand:
            p = primer.seq
            l = min(len(p), loc_3p)
            self.length = count_while(iupac.base_match(template[loc_3p-i],p[len(p)-1-i]) for i in range(l))
            assert(self.length > 0)

            self.loc_5p = self.loc_3p - self.length + 1
            self.left =  self.loc_5p
            self.right = self.loc_3p + 1
        else:
            p = primer.seq.reverse_complement()
            l = min(len(p), len(template)-loc_3p)
            self.length = count_while(iupac.base_match(template[loc_3p+i],p[i]) for i in range(l))
            assert(self.length > 0)

            self.loc_5p = self.loc_3p + self.length - 1
            self.left =  self.loc_3p
            self.right = self.loc_5p + 1
        self.full = self.length == len(self.primer)
        self.match = self.template[self.left:self.right]

    def __le__(self, rhs):
        return (self.left <= rhs.left) and (self.right <= rhs.right)

    def n_percent(self):
        count = self.template[self.left:self.right].count('N')
        return 1.*count/(self.right-self.left)

    def tm(self):
        return melt_temp.melting_temperature_unambiguous(self.template[self.left:self.right])


class Primer:
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

    def search(self, template, template_ambiguous=False, min_length=10):
        """
        Return tuple of fowards anneal locations and reverse anneal locations.
        anneal locations are (5'location, 3'location)
        """
        if min_length and len(self.seq) > min_length > 0:
            l = len(self.seq) - min_length
            seq = self.seq[l:]
        else:
            seq = self.seq
        primer = seq
        cprimer = seq.reverse_complement()
        reg = re.compile('(%s)|(%s)'%(iupac.oligo_regex(primer),iupac.oligo_regex(cprimer)))

        template = str(template).upper()
        pp = []
        pc = []
        start = 0
        while True:
            m = reg.search(template, start)
            if not m:
                break
            start = m.start()+1
            if m.group(1):
                pta = PrimerTemplateAnnealing(self, template, True, m.end()-1)
                if pta.tm() < 40.:
                    continue
                if pta.n_percent() > 0.5:
                    continue
                pp.append(pta)
            if m.group(2):
                pta = PrimerTemplateAnnealing(self, template, False, m.start())
                if pta.tm() < 40.:
                    continue
                if pta.n_percent() > 0.5:
                    continue
                pc.append(pta)
        return pp,pc

    def debugprint(self):
        #print("%s: %s, len=%2d, Tm=%.2f, oTm=%.2f, MarmurTm: %s, GC=%.2f, sa=%2d, sea=%2d" % (**self.get_table_row()))
        print('sa=%s, index=%s'%(self.sa.score, self.sa.index))
        print('\n'.join(self.sa.get_bar()))
        print('sea=%s, index=%s'%(self.sea.score, self.sea.index))
        print('\n'.join(self.sea.get_bar()))


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

    @property
    @memoize
    def score(self):
        return pcond.NORMAL_PCR.score_primer(self)

    @property
    @memoize
    def score_bisulfite(self):
        return pcond.BISULFITE_PCR.score_primer(self)


    def sdss_length(self, anneal_temp=None, pcr_mix=melt_temp.DEFAULT_MIX, unmethyl=False):
        s = str(self.seq).upper()
        if unmethyl:
            s = s.replace('R','A').replace('Y','T')
        else:
            s = s.replace('R','G').replace('Y','C')

        if not anneal_temp:
            anneal_temp = self.melting_temperature()
            
        l = len(s)
        for i in range(2,l):
            f = melt_temp.complex_fraction(s[l-i:-1], anneal_temp, pcr_mix)
            if f>0.01:
                break
        return i



class PrimerPair:
    def __init__(self, fw, rv):
        self.fw = fw
        self.rv = rv

    @property
    @memoize
    def pair_annealing(self):
        return PrimerAnnealing(self.fw.seq, self.rv.seq)
    @property
    @memoize
    def pair_end_annealing(self):
        return PrimerAnnealing(self.fw.seq, self.rv.seq, True)

    pa = pair_annealing
    pea = pair_end_annealing

    @property
    @memoize
    def score(self):
        return pcond.NORMAL_PCR.score_primerpair(self)

    @property
    @memoize
    def score_bisulfite(self):
        return pcond.BISULFITE_PCR.score_primerpair(self)

    def debugprint(self):
        print('score=', self.score)
        print('bisulfite score=', self.score_bisulfite)
        for r in [self.fw, self.rv]:
            r.debugprint()
        pa = self.pair_annealing
        print('pa=%s, index=%s'%(pa.score,pa.index))
        print('\n'.join(pa.get_bar()))
        pea = self.pair_end_annealing
        print('pea=%s, index=%s'%(pea.score,pea.index))
        print('\n'.join(pea.get_bar()))


    def write_html(self, w, pair_values=True, annealings=False):
        head = Primer.get_table_head()
        if pair_values:
            head = head + ['pa.', 'pea.', 'score', 'bs score']

        b = xmlwriter.builder(w)
        # primer table
        with b.div(cls='primerpair'):
            with b.table(cls='primerpairtable', border=1):
                with b.tr:
                    for p in head:
                        b.th(p)

                with b.tr:
                    for v in self.fw.get_table_row():
                        b.td(str(v))
                    if pair_values:
                        for v in [self.pa.score, self.pea.score, '%.2f'%self.score, '%.2f'%self.score_bisulfite]:
                            b.td(str(v),rowspan='2')
                with b.tr:
                    for v in self.rv.get_table_row():
                        b.td(str(v))

            b.p('Tm melting temperature(SantaLucia), oTm melting temperature for bisulfite methyl-template, sa. self annealing, sea. self end annealing, pa. pair annealing, pea. pair end annealing', style='font-size:x-small')
            if annealings:
                self.pair_annealing.write_html(w)
                self.pair_end_annealing.write_html(w)


class Primers(NamedList):
    def __init__(self):
        super(Primers,self).__init__()


    def load_file(self, filename):
        with open(filename,'r') as f:
            parser = TreekvParser()
            tree = parser.readfp(f, filename)

            for kv in list(tree.items()):
                self.append(Primer(kv.key, kv.value))


    def write_csv(self):
        ret = ""
        ret += ", ".join(Primer.get_table_head())
        ret += "\n"

        for r in self:
            ret += (", ".join([str(v) for v in r.get_table_row()]) + "\n")

        return ret+"\n"

    def write_html(self, w):
        b = xmlwriter.builder(w)
        # primer table
        with b.div(cls='primerpair'):
            with b.table(cls='primerpairtable', border=1):
                with b.tr:
                    for p in Primer.get_table_head():
                        b.th(p)

                for r in self:
                    with b.tr:
                        for v in r.get_table_row():
                            b.td(str(v))

            b.p('Tm melting temperature(SantaLucia), oTm melting temperature for bisulfite methyl-template, sa. self annealing, sea. self end annealing, pa. pair annealing, pea. pair end annealing', style='font-size:x-small')

    def get_default(self, name, default_name):
        try:
            return self[name]
        except KeyError:
            r = Primer(default_name, name)
            self.append(r)
            return r

def test_pp(a,b):
    PrimerPair(Primer('fw',a),Primer('rv',b)).debugprint()

