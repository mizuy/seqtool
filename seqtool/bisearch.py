from __future__ import absolute_import

from .nucleotide.primer import Primer, PrimerPair, PrimerCondition

import bisect
from numpy import *
from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from threading import Lock, Thread

__all__ = ['PrimerSearchOptions', 'dynamic_primer_search', 'adhoc_primer_search']

def irange(from_,to):
    return xrange(from_,to+1)
def irange_r(from_,to):
    return xrange(to,from_-1,-1)

class MeltingTemperature(object):
    def __init__(self, na_conc=33.*10**-3, mg_conc=(2+4.5)*10**-3, primer_conc=0.5*10**-6):
        self.na_conc = na_conc
        self.mg_conc = mg_conc
        self.primer_conc = primer_conc
        c_salt = self.na_conc + 4*(self.mg_conc**0.5)
        self._cc_salt = 16.6*log10(c_salt / (1.+0.7*c_salt)) - 269.3
        self._cc_primer = 10.2-0.001*1.987*298.2*log(self.primer_conc)

        def ddict_to_table(ddict):
            ret = empty((4*4),float)
            for k in ddict.keys():
                i0 = {'A':0, 'T':1, 'G':2, 'C':3}[k[0]]
                i1 = {'A':0, 'T':1, 'G':2, 'C':3}[k[1]]
                i = i0 + i1 * 4
                ret[i] = ddict[k]
            return ret

        self._nnt_dh = ddict_to_table({
                'AA': 9.1, 'TT': 9.1,
                'AT': 8.6,
                'TA': 6.0,
                'CA': 5.8, 'TG': 5.8,
                'GT': 6.5, 'AC': 6.5,
                'CT': 7.8, 'AG': 7.8,
                'GA': 5.6, 'TC': 5.6,
                'CG': 11.9,
                'GC': 11.1,
                'GG': 11.0, 'CC': 11.0,
                })
        self._nnt_dg = ddict_to_table({
                'AA': 1.55, 'TT': 1.55,
                'AT': 1.25,
                'TA': 0.85,
                'CA': 1.15, 'TG': 1.15,
                'GT': 1.40, 'AC': 1.40,
                'CT': 1.45, 'AG': 1.45,
                'GA': 1.15, 'TC': 1.15,
                'CG': 3.05,
                'GC': 2.70,
                'GG': 2.3,  'CC': 2.3,
                })

    def tm(self, enthalpy, free_energy):
        return self._cc_salt + 298.2*(10. + enthalpy) / (enthalpy-free_energy + self._cc_primer)

    def nn_enthalpy(self, n0, n1):
        '''
        @param n0 first base in int. {'A':0, 'T':1, 'G':2, 'C':3}
        @param n1 second base in int. {'A':0, 'T':1, 'G':2, 'C':3}
        '''
        return self._nnt_dh[n0+n1*4]
    def nn_free_energy(self, n0, n1):
        '''
        @param n0 first base in int. {'A':0, 'T':1, 'G':2, 'C':3}
        @param n1 second base in int. {'A':0, 'T':1, 'G':2, 'C':3}
        '''
        return self._nnt_dg[n0+n1*4]

    def seq_tm(self, sequence):
        '''
        @param sequence int[]
        '''
        l = len(sequence)
        dh = sum(self.nn_enthalpy(sequence[i],seqence[i+1]) for i in irange(0,l-2))
        dg = sum(self.nn_free_energy(sequence[i],seqence[i+1]) for i in irange(0,l-2))
        return self.tm(dh, dg)

class TargetSequence(object):
    def __init__(self, target, tm_calculator):
        assert isinstance(target,str)
        self.seqstr = target
        self.length = len(target)
        self.seqint = empty((self.length),int)
        self.seq = Seq.Seq(target,IUPAC.unambiguous_dna)
        for i,n in enumerate(target):
            self.seqint[i] = {'A':0, 'T':1, 'G':2, 'C':3}[n.upper()]

        self.tm_calc = tm_calculator

        self.anneals = zeros((4,4))
        self.anneals[0,1] = 2
        self.anneals[1,0] = 2
        self.anneals[2,3] = 4
        self.anneals[3,2] = 4

        self.anneals_c = zeros((4,4))
        self.anneals_c[0,0] = 2
        self.anneals_c[1,1] = 2
        self.anneals_c[2,2] = 4
        self.anneals_c[3,3] = 4

        self.tm_dh = empty((len(self.seqint)),float)
        self.tm_dg = empty((len(self.seqint)),float)
        self.tm_dh_u = empty((len(self.seqint)),float)
        self.tm_dg_u = empty((len(self.seqint)),float)
        for i in xrange(len(self.seqint)-1):
            i0 = self.seqint[i]
            i1 = self.seqint[i+1]
            self.tm_dh[i] = self.tm_calc.nn_enthalpy(i0,i1)
            self.tm_dg[i] = self.tm_calc.nn_free_energy(i0,i1)
            i0u = 1 if i0==3 else i0
            i1u = 1 if i1==3 else i1
            self.tm_dh_u[i] = self.tm_calc.nn_enthalpy(i0u,i1u)
            self.tm_dg_u[i] = self.tm_calc.nn_free_energy(i0u,i1u)
        
    def s_c(self,i0,i1):
        assert 0 <= i0 < self.length
        assert 0 <= i1 < self.length
        return self.anneals[self.seqint[i0],self.seqint[i1]]
    def s_c_c(self,i0,i1):
        assert 0 <= i0 < self.length
        assert 0 <= i1 < self.length
        return self.anneals_c[self.seqint[i0],self.seqint[i1]]
    def __len__(self):
        return self.length
    
    def gc_ratio(self,start,end):
        gc = 0
        for i in range(start,end):
            if self.seqint[i] >= 2:
                gc += 1
        return 1.*gc/(end-start)

    def is_gc(self,start):
        return 1 if (self.seqint[start] >= 2) else 0
    def is_cpg(self, index):
        return 1 if (index+1<self.length and self.seqint[index]==3 and self.seqint[index+1]==2) else 0
    def is_cpg_c(self, index):
        return 1 if (0<=index-1 and self.seqint[index]==2 and self.seqint[index-1]==3) else 0


class PrimerSearchOptions(object):
    MODE_NORMAL = 0
    MODE_BISULFITE_SEQUENCE = 1
    
    def __init__(self):
        self.primer_cond = PrimerCondition()
        self.min_product = 100
        self.max_product = 600
        self.threshold = 100.
        self.na_conc = 33.*10**-3
        self.mg_conc = (2+4.5)*10**-3
        self.primer_conc = 0.5*10**-6
        self.max_met_tm_diff = 2.5
        self.mode = PrimerSearchOptions.MODE_NORMAL
        self.max_tm_diff = 8.

    def cond_gc(self):
        if self.mode == PrimerSearchOptions.MODE_BISULFITE_SEQUENCE:
            return self.primer_cond.gc_bsp
        else:
            return self.primer_cond.gc

def replace_seq(seq,f,t):
    Seq.Seq(str(seq).replace(f,t),IUPAC.ambiguous_dna)

class PrimerPairResult(object):
    def __init__(self, target, primerp, score, i, n, j, m, pa, pa_k, pea, pea_k, bsp=False):
        self.target = target
        self.primerp = primerp
        self.score = score
        self.i = i
        self.n = n
        self.j = j
        self.m = m
        self.jm = j+m
        self.pa = pa
        self.pa_k = pa_k
        self.pea = pea
        self.pea_k = pea_k
        if bsp:
            self.fw = Primer('FW',replace_seq(target.seq[i:i+n],'C','Y'))
            self.rv = Primer('RV',replace_seq(target.seq[j:j+m].reverse_complement(),'G','R'))
        else:
            self.fw = Primer('FW',target.seq[i:i+n])
            self.rv = Primer('RV',target.seq[j:j+m].reverse_complement())
        self.primers = PrimerPair(self.fw, self.rv)

    def debug(self):
        print 'score=%s, i,j:%s,%s j,m:%s,%s, pa=%s, pa_k=%s pea=%s pea_k=%s'%(self.score,self.i,self.n,self.j,self.m,self.pa,self.pa_k,self.pea,self.pea_k)
        print self.primerp.desc_pos(self.i,self.n)
        print self.primerp.desc_neg(self.j,self.m)
        if self.pa_k is not None:
            print 'pa=%s'%self.pa
            if self.pa_k > 0:
                print ' '*abs(self.pa_k)+"5'-%s-3'"%self.fw
                print "3'-%s-5'"%self.rv.reverse()
            else:
                print "5'-%s-3'"%self.fw
                print ' '*abs(self.pa_k)+"3'-%s-5'"%self.rv.reverse()
            print ''
        if self.pea_k is not None:
            print 'pea=%s'%self.pea
            if self.pea_k > 0:
                print ' '*abs(self.pea_k)+"5'-%s-3'"%self.fw
                print "3'-%s-5'"%self.rv.reverse()
            else:
                print "5'-%s-3'"%self.fw
                print ' '*abs(self.pea_k)+"3'-%s-5'"%self.rv.reverse()
            print ''
        self.primers.debugprint()
        print 'validation....'
        self. validate()

    def validate(self):
        # primer check
        for p,index in [(self.fw,(self.i,self.n)), (self.rv,(self.j,self.m))]:
            assert abs(p.melting_temperature() - self.primerp.tm[index]) < 0.01
            assert p.gc_ratio() == self.primerp.gc[index]
            assert p.self_annealing().score == self.primerp.sa[index]
            
            if p is self.fw:
                assert p.self_end_annealing().score == self.primerp.sea_pos[index]
            else:
                assert p.self_end_annealing().score == self.primerp.sea_neg[index]
        # primer pair check
        assert self.primers.pair_annealing().score == self.pa
        assert self.primers.pair_end_annealing().score == self.pea

    def write(self, output):
        pass

class PPRGroup(object):
    def __init__(self, ppr):
        assert isinstance(ppr, PrimerPairResult)
        self.best = ppr
        self.score = self.best.score
        self.min_i = ppr.i
        self.max_i = self.min_i
        self.min_jm = ppr.j+ppr.m
        self.max_jm = self.min_jm
        pass
    def add(self, other):
        assert isinstance(other, PrimerPairResult)
        if other.score < self.score:
            self.best = other
            self.score = self.best.score
        self.min_i = min(self.min_i,other.i)
        self.max_i = max(self.max_i,other.i)
        self.min_jm = min(self.min_jm, other.jm)
        self.max_jm = max(self.max_jm, other.jm)
        
    def similar(self, other, threshold=5):
        assert isinstance(other, PrimerPairResult)
        
        v = min(abs(self.min_i-other.i), abs(self.max_i-other.i))
        w = min(abs(self.min_jm-other.jm), abs(self.max_jm-other.jm))
        if v<threshold and w<threshold:
            return True
        else:
            return False

class Results(object):
    def __init__(self, size=None, mode_bsp=True):
        '''
        size: size of results buffer. if None, no limits for number of results to store
        mode_bsp: True if primer pairs are for bisulfite treated gDNA. just for write seqviewfile
        '''
        self.pl = []
        self.size = size
        self.mode_bsp = mode_bsp
        self.count = 0
    def add(self, ppr):
        self.count += 1
        assert isinstance(ppr, PrimerPairResult)
        
        # search primer pairs which results in the similar pcr product
        # only the best pairs within the similar groups is added.
        sim = []
        for i, (s, p) in enumerate(self.pl):
            if p.similar(ppr):
                sim.append(p)
                del self.pl[i]
        if sim:
            worstsim = sim[0]
            for s in sim[1:]:
                if worstsim.score < s.score:
                    worstsim = s
            worstsim.add(ppr)
            pg = worstsim
        else:
            pg = PPRGroup(ppr)
        bisect.insort(self.pl,(pg.score,pg))

        if self.size:
            self.pl = self.pl[0:self.size]
        
    def add_result(self, other):
        for (s,ppr) in other.pl:
            self.add(ppr)
    def debug(self):
        if len(self.pl)==0:
            print 'no result.'
        for s,ppr in self.pl[::-1]:
            print '-'*100
            ppr.debug()
    def write_seqviewfile(self, fileobj):
        c = 0
        primers_fw = {}
        for i,(s,pg) in enumerate(self.pl):
            k = (pg.best.i, pg.best.n)
            if not primers_fw.has_key(k):
                pg.best.fw.name = 'bi-%03d-FW'%c
                primers_fw[k] = pg.best.fw
                c += 1
        c = 0
        primers_rv = {}
        for i,(s,pg) in enumerate(self.pl):
            k = (pg.best.j, pg.best.m)
            if not primers_rv.has_key(k):
                pg.best.rv.name = 'bi-%03d-RV'%c
                primers_rv[k] = pg.best.rv
                c += 1

        fileobj.write('>primer\n')
        for p in primers_fw.values():
            fileobj.write('%s: %s\n'%(p.name,str(p.seq)))
        for p in primers_rv.values():
            fileobj.write('%s: %s\n'%(p.name,str(p.seq)))
        fileobj.write('>pcr\n')
        if self.mode_bsp:
            fileobj.write('@template=bisulfite_met\n')
        else:
            fileobj.write('@template=genome\n')
        for i,(s,pg) in enumerate(self.pl):
            fileobj.write('bi_%03d, %.1f: %s, %s\n'%(i, s, primers_fw[(pg.best.i, pg.best.n)].name, primers_rv[(pg.best.j, pg.best.m)].name))

def max_k(function, range_):
    ret = None
    ret_k = None
    for k in range_:
        v = function(k)
        if (not ret) or ret<v:
            ret = v
            ret_k = k
    return ret,ret_k

def dynamic_primer_search(target_seq, option=PrimerSearchOptions(), thread_num=1):
    '''
    primer search
    @param target_seq Seq.Seq object of target sequence. original sense or antisense sequence (before bisulfite treated)
    @param option PrimerSearchOptions object
    '''
    tm_calculator = MeltingTemperature(option.na_conc, option.mg_conc, option.primer_conc)
    
    target = TargetSequence(str(target_seq),tm_calculator)
    seq_length = len(target)
    primer_length = option.primer_cond.primer_length.maximum
    print 'sequence length=',seq_length

    mode_bsp = True if option.mode == PrimerSearchOptions.MODE_BISULFITE_SEQUENCE else False

    class PrimerProperties(object):
        def __init__(self, seq_len, primer_len):
            # +1 strand (sense strand)
            self.tm = empty((seq_len,primer_len+1), dtype=float)
            self.tm[:]=-1
            self.tm_u = empty((seq_len,primer_len+1), dtype=float)
            self.tm_u[:]=-1

            self.gc = empty((seq_len,primer_len+1), dtype=float)
            self.gc[:]=-1

            self.sa = empty((seq_len,primer_len+1), dtype=int)
            self.sa[:]=-1

            self.sa_k = empty((seq_len,primer_len+1), dtype=int)
            self.sa_k[:]=-1

            self.cpg = empty((seq_len,primer_len+1), dtype=int)
            self.cpg[:] = 0

            self.sea_pos = empty((seq_len,primer_len+1), dtype=int)
            self.sea_pos[:]=-1
            self.sea_pos_k = empty((seq_len,primer_len+1), dtype=int)
            self.sea_pos_k[:]=0
            self.sea_neg = empty((seq_len,primer_len+1), dtype=int)
            self.sea_neg[:]=-1
            self.sea_neg_k = empty((seq_len,primer_len+1), dtype=int)
            self.sea_neg_k[:]=0

            self.valid = empty((seq_len,primer_len+1), dtype=bool)
            self.valid[:]=True

            self.p_pos = empty((seq_len,primer_len+1), dtype=float)
            self.p_pos[:]=-1
            self.p_neg = empty((seq_len,primer_len+1), dtype=float)
            self.p_neg[:]=-1
            
        def set_invalid(self, i,n):
            self.valid[i,n]=False
        def is_valid(self,i,n):
            return self.valid[i,n]

        def desc(self,s,i,n):
            return '%s, len=%2d, Tm=%.2f, oTm=%.2f, GC=%.2f, sa=%2d, sa_k=%2d, CpG=%2d'%(("5'-%s-3'"%s).ljust(30), n, self.tm[i,n], self.tm_u[i,n], self.gc[i,n], self.sa[i,n], self.sa_k[i,n], self.cpg[i,n])
        
        def desc_pos(self,i,n):
            s = target.seq[i:i+n]
            return 'FW:'+self.desc(s,i,n)+', sea=%2d, sea_k=%2d'%(self.sea_pos[i,n], self.sea_pos_k[i,n])
        def desc_neg(self,i,n):
            s = target.seq[i:i+n].reverse_complement()
            return 'RV:'+self.desc(s,i,n)+', sea=%2d, sea_k=%2d'%(self.sea_neg[i,n], self.sea_neg_k[i,n])
        
    primerp = PrimerProperties(seq_length,primer_length)
    
    # Tm, length, GC calculation
    print 'calculating Tm, length, GC...'
    for i in xrange(seq_length):
        c_cpg = target.is_cpg(i)
        c_gc = target.is_gc(i)
        dh = 0.
        dg = 0.
        dh_u = 0.
        dg_u = 0.
        for n in irange(2,min(seq_length-i,primer_length)):
            assert i+n <= seq_length
            end = i+n-1
            c_cpg += target.is_cpg(end)
            c_gc += target.is_gc(end)
            dh += target.tm_dh[end-1]
            dg += target.tm_dg[end-1]
            dh_u += target.tm_dh_u[end-1]
            dg_u += target.tm_dg_u[end-1]

            # length
            if not option.primer_cond.primer_length.bound(n):
                primerp.set_invalid(i,n)
                continue
            
            # num of CpG in primer
            if mode_bsp and c_cpg >= 2:
                primerp.set_invalid(i,n)
                continue
            primerp.cpg[i,n] = c_cpg
            
            # Tm(methylated)
            tm = tm_calculator.tm(dh,dg)
            tm_u = tm_calculator.tm(dh_u,dg_u)
            if (not option.primer_cond.tm.bound(tm)) or (not option.primer_cond.tm.bound(tm_u)):
                primerp.set_invalid(i,n)
                continue
            if mode_bsp and abs(tm-tm_u)>option.max_met_tm_diff:
                primerp.set_invalid(i,n)
                continue
            primerp.tm[i,n] = tm
            primerp.tm_u[i,n] = tm_u

            # GC ratio
            gc = 100.*c_gc/n
            if not option.cond_gc().bound(gc):
                primerp.set_invalid(i,n)
                continue
            primerp.gc[i,n] = gc

    # x=table[startn(1,3)] : x is table value of 'v1,v2,v3' where v = sequence
    def startend(start,end):
        return (start, end-start+1)

    # self annealing calculation
    # sa(i...i+n-1) = max{k from -(n-1) to n-1} sa_k(i...i+n-1)
    # sa_k(i...i+n-1) = k<0 | sa_0(i.....i+n-1+k)
    #                   k=0 | sa_0(i.....i+n-1)
    #                   k>0 | sa_0(i+k...i+n-1)
    # sa first step, store all sa_0 value for sa
    # sa_0(i...i+n-1) = n=1 | score_c(i,i)
    #                   n=2 | 2*score_c(i,i+1)
    #                   n=n | sa_0(i+1,i+n-1-1) + 2*score_c(i,i+n-1)
    print 'calculating self annealing...'
    sa_0 = empty((seq_length,primer_length+1),dtype=int)
    sa_0[:] = -1
    for end in xrange(seq_length):
        for n in irange(1,min(end+1,primer_length)):
            i = end-n+1
            assert i+n <= seq_length
            # fill sa_0[i, n]
            if n==1:
                assert i==end, 'from %s to %s, n=%s'%(i,end,n)
                sa_0[i,n] = target.s_c(i,i)
            elif n==2:
                sa_0[i,n] = 2*target.s_c(i,end)
            else:
                assert sa_0[startend(i+1,end-1)]>=0, 'from %s to %s, n=%s'%(i,end,n)
                sa_0[i,n] = sa_0[startend(i+1,end-1)] + 2*target.s_c(i,end)
    # sa second step, 
    for i in xrange(seq_length):
        for n in irange(1,min(seq_length-i,primer_length)):
            assert i+n <= seq_length
            if not primerp.is_valid(i,n):
                continue
            sa, sa_k = max_k(lambda k:sa_0[(i,n+k) if k<0 else (i+k,n-k)], irange(-(n-1),n-1))
            primerp.sa[i,n] = sa
            primerp.sa_k[i,n] = sa_k
            if not option.primer_cond.sa.bound(sa):
                primerp.set_invalid(i,n)
        
    # Self End Annealing calculation
    # sea(i...i+n-1) = max{k from 0 to n-1} sea_k(i...i+n-1)
    # sea_k(i...i+n-1) = sea_0(i+k...i+n-1)
    # sea_0(i...i+n-1) = n=1 | score_c(i,i)
    #                    n=2 | 2*score_c(i,i+1)
    #                    n=n | 
    print 'calculating self end annealing...'
    sea_0 = empty((seq_length,primer_length+1),dtype=int)
    sea_0[:] = -1
    sea_0_full = empty((seq_length,primer_length+1),dtype=bool)
    sea_0_full[:] = None
    for end in xrange(seq_length):
        for n in irange(1,min(end+1,primer_length)):
            i = end-n+1
            assert 0 <= i <= end
            assert i+n <= seq_length
            # fill sea_0[i, n]
            if n==1:
                assert i==end, 'from %s to %s, n=%s'%(i,end,n)
                sea_0[i,n] = target.s_c(i,i)
                sea_0_full[i,n] = True
            elif n==2:
                s = target.s_c(i,end)
                sea_0[i,n] = 2*s
                sea_0_full[i,n] = True if s else False
            else:
                base = sea_0[startend(i+1,end-1)]
                assert base>=0, 'from %s to %s, n=%s'%(i,end,n)
                s = target.s_c(i,end)
                full = base>0 and s>0 and sea_0_full[startend(i+1,end-1)]
                sea_0_full[i,n] =full
                if full:
                    sea_0[i,n] = base+2*s
                else:
                    sea_0[i,n] = (base+s) if s else 0
                #print "1st: i=%d, %d=sea_0[%s,%s], %s, full=%s, base=%s, s=%s"%(i,v,i,n,target.seqstr[i:i+n],full,base,s)
                    
    # sea second step, 
    for i in xrange(seq_length):
        for n in irange(1,min(seq_length-i,primer_length)):
            assert i+n <= seq_length
            if not primerp.is_valid(i,n):
                continue
            pos, pos_k = max_k(lambda k:sea_0[i+k, n-k],irange(0,n-1))
            neg, neg_k = max_k(lambda k:sea_0[i, n+k],irange(-(n-1),0))
            primerp.sea_pos[i,n] = pos
            primerp.sea_pos_k[i,n] = pos_k
            primerp.sea_neg[i,n] = neg
            primerp.sea_neg_k[i,n] = neg_k
            if not (option.primer_cond.pea.bound(pos) or option.primer_cond.pea.bound(neg)):
                primerp.set_invalid(i,n)


    for i in xrange(seq_length):
        for n in irange(2,min(seq_length-i,primer_length)):
            assert i+n <= seq_length
            c = option.primer_cond.primer_length.score(n) + option.cond_gc().score(primerp.gc[i,n]) + option.primer_cond.sa.score(primerp.sa[i,n]) + option.primer_cond.tm.score(primerp.tm[i,n])
            pos = c + option.primer_cond.sea.score(primerp.sea_pos[i,n])
            neg = c + option.primer_cond.sea.score(primerp.sea_neg[i,n])
            if pos > option.threshold or neg > option.threshold:
                primerp.set_invalid(i,n)
            primerp.p_pos[i,n] = pos
            primerp.p_neg[i,n] = neg
    
    # Pair Annealing
    pa_0 = empty((seq_length, seq_length, primer_length+1),dtype=int)
    pa_0[:] = -1
    pea_0_l = empty((seq_length, seq_length, primer_length+1),dtype=int)
    pea_0_l[:] = -1
    pea_0_l_full = empty((seq_length, seq_length, primer_length+1),dtype=bool)
    pea_0_l_full[:] = True
    pea_0_r = empty((seq_length, seq_length, primer_length+1),dtype=int)
    pea_0_r[:] = -1

    # pair annealing values for each pair w=target[i:i+n], v=target[j:j+m]
    # pa(w,v) = max{k from -(n-1) to m-1} 
    #                 pa_0[1,1,n+k] if k<0, n+k<m
    # lazy calculation of pa0
    print 'calculating pair annealing...'

    lock_pa = Lock()
    lock_pea = Lock()

    result = [Results(None, mode_bsp) for i in xrange(thread_num)]

    def multithread_division_problem(problem, threadno, threadsnum):
        block = problem/threadsnum
        if threadno == threadsnum-1:
            return threadno*block, problem
        else:
            return threadno*block, (threadno+1)*block4

    def pair_calc(threadno, threadsnum):
        #for i in xrange(seq_length):
        def get_pa0(i,j,l):
            if pa_0[i,j,l] < 0:
                if l==1:
                    s = target.s_c_c(i,j)
                    with lock_pa:
                        pa_0[i,j,l] = s
                else:
                    assert get_pa0(i,j,l-1)>=0
                    s = get_pa0(i,j,l-1) + target.s_c_c(i+l-1,j+l-1)
                    with lock_pa:
                        pa_0[i,j,l] = s
            return pa_0[i,j,l]

        def get_pea0(i,j,l):
            if pea_0_l[i,j,l] < 0:
                if l==1:
                    s = target.s_c_c(i,j)
                    with lock_pea:
                        pea_0_l_full[i,j,l] = (s!=0)
                        pea_0_l[i,j,l] = s
                        pea_0_r[i,j,l] = s
                else:
                    left,right = get_pea0(i,j,l-1)
                    assert left>=0 and right>=0
                    s = target.s_c_c(i+l-1,j+l-1)
                    full = pea_0_l_full[i,j,l-1] and (s!=0)
                    newleft = left+s if full else left
                    newright = right+s if s else 0
                    with lock_pea:
                        pea_0_l_full[i,j,l] = full
                        pea_0_l[i,j,l]=newleft
                        pea_0_r[i,j,l]=newright
            return (pea_0_l[i,j,l],pea_0_r[i,j,l])

        starti, endi = multithread_division_problem(seq_length, threadno, threadsnum)

        for i in xrange(starti, endi):
            for n in irange(2,min(seq_length-i,primer_length)):
                assert i+n <= seq_length
                if not primerp.is_valid(i,n):
                    continue
                for j in irange(i+option.min_product-primer_length,min(i+option.max_product-primer_length,seq_length)):
                    for m in irange(2,min(seq_length-j,primer_length)):
                        if not (option.min_product <= j+m-i+1 <= option.max_product):
                            continue
                        if not primerp.is_valid(j,m):
                            continue
                        if abs(primerp.tm[i,n]-primerp.tm[j,m]) > option.max_tm_diff:
                            continue
                        score = primerp.p_pos[i,n]+primerp.p_neg[j,m]
                        if score > option.threshold:
                            continue
                        def calcpa(k):
                            if n<=m:
                                if k<=0:
                                    return get_pa0(i+max(0,-k), j+max(0,k), n+k)
                                elif k<=m-n:
                                    return get_pa0(i+max(0,-k), j+max(0,k), n)
                                else:
                                    return get_pa0(i+max(0,-k), j+max(0,k), m-k)
                            else:
                                if k<=m-n:
                                    return get_pa0(i+max(0,-k), j+max(0,k), n+k)
                                elif k<=0:
                                    return get_pa0(i+max(0,-k), j+max(0,k), m)
                                else:
                                    return get_pa0(i+max(0,-k), j+max(0,k), m-k)
                        def calcpea(k):
                            if n<=m:
                                if k<=0:
                                    l,r = get_pea0(i+max(0,-k), j+max(0,k), n+k)
                                    return max(l,r)
                                elif k<=m-n:
                                    l,r = get_pea0(i+max(0,-k), j+max(0,k), n)
                                    return r
                            else:
                                if k<=m-n:
                                    l,r = get_pea0(i+max(0,-k), j+max(0,k), n+k)
                                    return max(l,r)
                                elif k<=0:
                                    l,r = get_pea0(i+max(0,-k), j+max(0,k), m)
                                    return l
                                return 0

                        pa,pa_k = max_k(calcpa, irange(-(n-1),m-1))
                        if not option.primer_cond.pa.bound(pa):
                            continue
                        pea,pea_k = max_k(calcpea, irange(-(n-1),m-1))
                        if not option.primer_cond.pea.bound(pea):
                            continue
                        score += option.primer_cond.pa.score(pa)
                        score += option.primer_cond.pea.score(pea)
                        ppr = PrimerPairResult(target, primerp, score,i,n,j,m,pa,pa_k,pea,pea_k,mode_bsp)
                        result[threadno].add(ppr)

    if thread_num>1:
        threads = []
        for threadno in xrange(thread_num):
            th = Thread(target=pair_calc, args=(threadno,thread_num))
            threads.append(th)
            th.start()
        for t in threads:
            print 'join...',t
            t.join()
            print 'end',t
        for i in result[1:]:
            result[0].add_result(i)
    else:
        pair_calc(0,1)

    '''
    total = 0
    used = 0
    for i in xrange(seq_length):
        for j in xrange(seq_length):
            for n in xrange(primer_length):
                total += 1
                if pa_0[i,j,n] >= 0:
                    used += 1
    print 100.*used/total
    '''
                                  
    return result[0]

def gen_primer_pairs(target,length_min=10,length_max=30):
    '''
    [target[i,j] for all i,j where length_min <= len(target[i,j]) <= length_max]
    '''
    l = len(target)
    for p in xrange(l):
        for q in xrange(p+length_min,min(p+length_max,l)):
            for r in xrange(l):
                for s in xrange(r+length_min,min(r+length_max,l)):
                    yield target[p:q],target[r:s]


def adhoc_primer_search(sequence):
    '''
    ad hoc primer search.
    VERY slow.
    just for experiment
    '''
    pl = []
    for a,b in gen_primer_pairs(sequence,10,30):
        s = PrimerPair(Primer(None,a),Primer(None,b)).score()
        bisect.insort(pl,(s,a,b))
        pl = pl[0:50]


def main():
    import sys, os
    from optparse import OptionParser

    parser = OptionParser('usage: %prog [options] input.fasta -o outputfile')

    parser.add_option('-o', '--output', dest='output', help='output filename')
    parser.add_option('--min_product', type='int', dest='min_product', default=100)
    parser.add_option('--max_product', type='int', dest='max_product', default=600)
    parser.add_option('--min_primer_length', type='int', dest='min_primer_length', default=10)
    parser.add_option('--max_primer_length', type='int', dest='max_primer_length', default=30)
    parser.add_option('-m', '--mode', dest='mode', default='normal', help='bisearch mode: normal or bisulfite [default: %default]')

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error('no input file')
        exit()

    option = PrimerSearchOptions()
    option.min_product = options.min_product
    option.max_product = options.max_product
    option.primer_cond.primer_length.minimum = options.min_primer_length
    option.primer_cond.primer_length.maximum = options.max_primer_length
    if options.mode == 'normal':
        option.mode = PrimerSearchOptions.MODE_NORMAL
    elif options.mode == 'bisulfite':
        option.mode = PrimerSearchOptions.MODE_BISULFITE_SEQUENCE
    else:
        parser.error('unkown mode: %s'%options.mode)

    if options.output:
        output = open(options.output,'w')
    else:
        output = sys.stdout

    for i in args:
        for record in SeqIO.parse(open(i,'r'), 'fasta'):
            name = '%s: %s'%(i, record.id)
            target = record.seq

            output.write('#%s\n'%name)
            result = dynamic_primer_search(target,option)
            result.write_seqviewfile(output)

if __name__=='__main__':
    main()
