import re

from . import to_seq, melt_temp, tm_gc, is_sequence_like
from ..util.memoize import memoize
from ..util import xmlwriter
from ..util.parser import TreekvParser
from .cpg import gc_ratio
from ..util.namedlist import NamedList
from . import iupac
from . import primer_cond as pcond

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
    '''
    >>> count_while([1,2,3,4,0])
    4
    >>> count_while([])
    0
    >>> count_while([0])
    0
    '''
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


        if strand

            a_l  a_r
             |    |
          5'-ATGCTCCGTATG-3'
                  CCGTATG
                  |      |
                left    right

        if not strand

                left    right
                  |       |
                  ATGCCATG
                  ATGCCATGAACAT-5'
                          |   |
                         a_l a_r


        """
        self.primer = primer
        self.template = template
        self.strand = strand

        self.loc_3p = loc_3p

        if self.strand:
            p = primer.seq
            l = min(len(p), loc_3p)
            self.length = count_while(iupac.base_match(template[loc_3p-i],p[len(p)-1-i]) for i in range(l))
            self.loc_5p = self.loc_3p - self.length + 1
            self.left =  self.loc_5p
            self.right = self.loc_3p + 1

            assert(self.length > 0)
            pp =  str(primer.seq)
            self.adapter_length = len(p)-self.length
            self.primer_match =  pp[self.adapter_length:]
            self.display_match = pp[self.adapter_length:]
            self.display_adapter = pp[:self.adapter_length]

            self.a_l = self.left -  self.adapter_length
            self.a_r = self.left
            self.leftmost =  self.a_l
        else:
            p = primer.seq.reverse_complement()
            l = min(len(p), len(template)-loc_3p)
            self.length = count_while(iupac.base_match(template[loc_3p+i],p[i]) for i in range(l))
            self.loc_5p = self.loc_3p + self.length - 1
            self.left =  self.loc_3p
            self.right = self.loc_5p + 1

            assert(self.length > 0)
            pp =  str(primer.seq)
            self.adapter_length = len(p)-self.length
            self.primer_match =  pp[self.adapter_length:]
            self.display_match = pp[self.adapter_length:][::-1]
            self.display_adapter = pp[:self.adapter_length][::-1]

            self.a_l = self.right
            self.a_r = self.right + self.adapter_length
            self.leftmost =  self.left


        self.full = (self.length == len(self.primer))
        self.match = self.template[self.left:self.right]


    def __le__(self, rhs):
        return (self.left <= rhs.left) and (self.right <= rhs.right)

    def n_percent(self):
        count = self.template[self.left:self.right].count('N')
        return 1.*count/(self.right-self.left)

    def tm(self, pcr_mix=melt_temp.STANDARD_MIX):
        return melt_temp.melting_temperature_unmethyl(self.primer_match, pcr_mix)


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

    def melting_temperature(self, pcr_mix=melt_temp.STANDARD_MIX, unmethyl=True):
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

    def search(self, template, template_ambiguous=False, min_length=16):
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
        if template_ambiguous:
            oregex = lambda x: iupac.oligo_regex(x, iupac.basematch_partial)
        else:
            oregex = lambda x: iupac.oligo_regex(x, iupac.basematch_unambiguous)
            
        reg = re.compile('(%s)|(%s)'%(oregex(primer),oregex(cprimer)))

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

    def write_text(self):
        print("-"*20)
        print("{}: 5'-{}-3'".format(self.name,self.seq))
        print("reverse_complement: 5'-{}-3'".format(self.seq.reverse_complement()))
        h = self.get_table_head()[2:8]
        r = self.get_table_row()[2:8]
        for n,v in zip(h,r):
            print('{}: {}'.format(n,v))
        print('sa=%s, index=%s'%(self.sa.score, self.sa.index))
        print('\n'.join(self.sa.get_bar()))
        print('sea=%s, index=%s'%(self.sea.score, self.sea.index))
        print('\n'.join(self.sea.get_bar()))
        print('check:')
        print('   longest continuous seq: TODO')
        print('   contains GGGG?: {}'.format(bool(self.seq.find('GGGG')>=0)))
        p35 = self.seq[-4:]
        p35gc = p35.count('C') + p35.count('G')
        print("   5bases at 3'end has GC <= 2 ?: {} ({}bp)".format(bool(p35gc<=2),p35gc))


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


class TaqmanProbe(Primer):
    def write_text(self):
        super().write_text()
        c = self.seq.count('C')
        g = self.seq.count('G')
        print('   Cs >Gs ?: {} ({} > {})'.format(c>g, c, g))
        print("   5' not G?: {}".format(bool(self.seq[0]!='G')))

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

    def write_text(self):
        print('PrimerPair: {} {}'.format(self.fw.name, self.rv.name))
        print('score=', self.score)
        print('bisulfite score=', self.score_bisulfite)
        print('pair-annealing:')
        pa = self.pair_annealing
        print('pa=%s, index=%s'%(pa.score,pa.index))
        print('\n'.join(pa.get_bar()))
        pea = self.pair_end_annealing
        print('pea=%s, index=%s'%(pea.score,pea.index))
        print('\n'.join(pea.get_bar()))
        for r in [self.fw, self.rv]:
            r.write_text()


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
            if is_sequence_like(name):
                r = Primer(default_name, name)
                self.append(r)
                return r
        raise KeyError('no such primer and this is not valid sequence: {}'.format(name))


def main():
    import sys
    ps = sys.argv[1:]

    if len(ps)==2:
        PrimerPair(Primer('fw',ps[0]),Primer('rv',ps[1])).write_text()
    else:
        for i,p in enumerate(ps):
            Primer('Primer{}'.format(i),p).write_text()


def probe():
    import sys
    ps = sys.argv[1:]

    for i,p in enumerate(ps):
        s = to_seq(p)
        TaqmanProbe('Probe_{}'.format(i),s).write_text()
        TaqmanProbe('ProbeR_{}'.format(i),s.reverse_complement()).write_text()
    
