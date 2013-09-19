import re

from . import to_seq, melt_temp, tm_gc, is_sequence_like
from ..util.memoize import memoize
from ..util import xmlwriter
from ..util.parser import TreekvParser
from .cpg import gc_ratio
from ..util.namedlist import NamedList
from . import iupac
from . import primer_cond as pcond
from .annealing import OligoAnnealing

__all__ = ['Primer', 'PrimerPair']

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

class CoordinateTransform:
    @classmethod
    def create_from_loc3p(cls, length, loc_3p, same_orientation=True):
        if same_orientation:
            loc_5p = loc_3p - (length - 1)
        else:
            loc_5p = loc_3p + (length - 1)
        return CoordinateTransform(loc_5p, same_orientation)
        
    def __init__(self, loc_5p, same_orientation=True):
        """
        ########### same_orientation = True #############
              012345678901234567         (seq_coordinate)
           5'-....ATGCCATG......-3'      (seq)
           0123456789012345678901234     (universal coordinate)
              |                |
             loc_5p=3        loc_3p=20

        universal_coordinate = 1 * seq_coordinate + loc_5p

        ########### same_orientation = False #############
                32109876543210           (seq_coordinate)
             3'-..ATGCCATG....-5'        (seq)
           0123456789012345678901234     (universal coordinate)
                |            |
            loc_3p=5       loc_5p=18

        universal_coordinate = -1 * seq_coordinate + loc_5p

        >>> c = CoordinateTransform(3, True)
        >>> c.get_u(0), c.get_u(10)
        (3, 13)
        >>> c.get_l(3), c.get_l(13)
        (0, 10)
        >>> c = CoordinateTransform(18, False)
        >>> c.get_u(0), c.get_u(10)
        (18, 8)
        >>> c.get_l(18), c.get_l(8)
        (0, 10)

        >>> c = CoordinateTransform.create_from_loc3p(18, 20, True)
        >>> c.get_u(0), c.get_u(10)
        (3, 13)
        >>> c.get_l(3), c.get_l(13)
        (0, 10)
        >>> c = CoordinateTransform.create_from_loc3p(14, 5, False)
        >>> c.get_u(0), c.get_u(10)
        (18, 8)
        >>> c.get_l(18), c.get_l(8)
        (0, 10)
        
        """
        self.a = 1 if same_orientation else -1
        self.b = loc_5p

    def get_u(self, local_location):
        return self.a * local_location + self.b
    def get_l(self, universal_location):
        return self.a * (universal_location - self.b)
        
    
class PrimerTemplateAnnealing:
    def __init__(self, primer, template, strand, loc_3p):
        """
           5'-ATGCATGCCATG-3'            (fw primer)
        5'-.......ATGCCATG..........-3'  (template)
               3'-TACGGTACATGC-5'        (rv primer)

        all location is template geometry

        ############### if strand = True (fw primer) ###############

           5'-ATGCATGCCATG-3'            (fw primer)
        5'-.......ATGCCATG..........-3'      (template)
           0123456789012345678901234
              |   |       |
            left middle  right
                 m_l     m_r
            a_l  a_r
                         |
                      loc_3p
        
        ############## if strand = False (rv primer) ###############

           0123456789012345678901234
        5'-.......ATGCCATG..........-3'  (template)
               3'-TACGGTACATGC-5'        (rv primer)
                  |       |   | 
                left   middle right
                m_l      m_r
                         a_l a_r
                  |
                loc_3p


        >>> fw = Primer('primer', 'ATGCATGCCATG')
        >>> rv = Primer('primer', 'CGTACATGGCAT')
        >>> template = to_seq('nnntacgATGCCATGatgcnnnnnn')
        >>> a = PrimerTemplateAnnealing(fw, template, True, 14)
        >>> a.match_left, a.match_right
        (7, 15)
        >>> a.adapter_left, a.adapter_right
        (3, 7)
        >>> a.match_length, a.adapter_length
        (8, 4)
        >>> a = PrimerTemplateAnnealing(rv, template, False, 7)
        >>> a.match_left, a.match_right
        (7, 15)
        >>> a.adapter_left, a.adapter_right
        (15, 19)
        >>> a.match_length, a.adapter_length
        (8, 4)
        """
        self.primer = primer
        self.template = template
        self.strand = strand

        #t_coord = CoordinateTransform(0, True)
        p_coord = CoordinateTransform.create_from_loc3p(len(self.primer), loc_3p, strand)
        def get(location_t):
            return (template[location_t], primer.seq[p_coord.get_l(location_t)])
        
        if self.strand:
            p = primer.seq
            l = min(len(p), loc_3p+1)

            self.match_length = count_while(
                iupac.base_match(template[i],p[p_coord.get_l(i)]) for i in range(loc_3p, loc_3p-l, -1) )
            assert(self.match_length > 0)
            self.adapter_length = len(p) - self.match_length

            self.leftmost =  loc_3p - len(self.primer) + 1
            self.middle = loc_3p - self.match_length + 1
            self.rightmost = loc_3p + 1

            self.adapter_left = self.leftmost
            self.adapter_right = self.middle
            self.match_left = self.middle
            self.match_right = self.rightmost
            
            #self.loc_3p = loc_3p
            #self.loc_5p = self.middle

            self.primer_match =  primer.seq[self.adapter_length:]
            self.primer_adapter = primer.seq[:self.adapter_length]
            self.primer_match_sense = self.primer_match
            self.primer_adapter_sense = self.primer_adapter
            self.display_match = str(self.primer_match)
            self.display_adapter = str(self.primer_adapter)
        else:
            p = primer.seq.complement()
            l = min(len(p), len(template)-loc_3p+1)
            self.match_length = count_while(
                iupac.base_match(template[i],p[p_coord.get_l(i)]) for i in range(loc_3p,loc_3p+l))
            assert(self.match_length > 0)
            self.adapter_length = len(p) - self.match_length

            self.leftmost = loc_3p
            self.middle = loc_3p + self.match_length
            self.rightmost = loc_3p + len(self.primer)

            self.match_left = self.leftmost
            self.match_right = self.middle
            self.adapter_left = self.middle
            self.adapter_right = self.rightmost

            assert(self.adapter_length == self.adapter_right - self.adapter_left)
            
            #self.loc_3p = loc_3p
            #self.loc_5p = self.middle - 1
            
            self.primer_match = primer.seq[self.adapter_length:]
            self.primer_adapter = primer.seq[:self.adapter_length]
            self.primer_match_sense = self.primer_match.reverse_complement()
            self.primer_adapter_sense = self.primer_adapter.reverse_complement()
            self.display_match = str(self.primer_match)[::-1]
            self.display_adapter = str(self.primer_adapter)[::-1]

        self.full = (self.match_length == len(self.primer))
        self.match = self.template[self.match_left:self.match_right]

    def __repr__(self):
        return 'PrimerTemplateAnnealing({}, match({}bp):{} -> {}, adapter({}bp):{} -> {})'.format(self.strand, self.match_length, self.match_left, self.match_right, self.adapter_length, self.adapter_left, self.adapter_right)

    def __le__(self, rhs):
        return (self.match_left <= rhs.match_left) and (self.match_right <= rhs.match_right)

    def n_percent(self):
        count = self.template[self.match_left:self.match_right].count('N')
        return 1.*count/self.match_length

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

    @memoize
    def annealings(self):
        return OligoAnnealing(self.seq, self.seq)
    @property
    def sa(self):
        return self.annealings().max_annealing
    @property
    def sea(self):
        return self.annealings().end_annealing
    self_annealing = sa
    self_end_annealing = sea

    def search(self, template, template_ambiguous=False, min_length=16):
        """
        Return tuple of fowards anneal locations and reverse anneal locations.
        anneal locations are (5'location, 3'location)
        """
        if min_length > 0 and len(self.seq) > min_length:
            seq = self.seq[-min_length:]
        else:
            seq = self.seq
        primer = seq
        cprimer = seq.reverse_complement()

        if template_ambiguous:
            oregex = lambda x: iupac.oligo_regex(x, iupac.basematch_subset)
        else:
            oregex = lambda x: iupac.oligo_regex(x, iupac.basematch_unambiguous)
            
        reg = re.compile('(%s)|(%s)'%(oregex(primer),oregex(cprimer)))

        template_str = str(template).upper()
        pp = []
        pc = []
        start = 0
        while True:
            m = reg.search(template_str, start)
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
        assert(isinstance(fw, Primer))
        assert(isinstance(rv, Primer))
        self.fw = fw
        self.rv = rv

    @memoize
    def annealings(self):
        return OligoAnnealing(self.fw.seq, self.rv.seq)
    @property
    def pa(self):
        return self.annealings().max_annealing
    @property
    def pea(self):
        return self.annealings().end_annealing
    pair_annealing = pa
    pair_end_annealing = pea
        
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
        pa = self.pa
        print('pa=%s, index=%s'%(pa.score,pa.index))
        print('\n'.join(pa.get_bar()))
        pea = self.pea
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
                self.pa.write_html(w)
                self.pea.write_html(w)


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
    
