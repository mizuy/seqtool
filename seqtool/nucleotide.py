from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from math import log, log10

import re, operator

def cpg_obs_per_exp(seq):
    """
    CpG islands in vertebrate genomes {Gardiner-Garden, 1987, p02206}
    'Obs/Exp CpG' = N * 'Number of CpG' / 'Number of C' * 'Number of G'
    where, N is the total number of nucleotide in the sequence being analyzed.
    """
    n = len(seq)
    c = 0
    g = 0
    cpg = 0
    for i,b in enumerate(seq):
        b = b.upper()
        if b=='C': c+=1
        if b=='G': g+=1
        if i+1<n and b=='C' and seq[i+1]=='G':
            cpg += 1

    if c*g == 0:
        return 0
    return 1.*n*cpg / (c*g)


def cpg_sites(seq):
    return [i for i in range(0,len(seq)-1) if seq[i]=='C' and seq[i+1]=='G']

def base_color(n):
    return {'A':'#00FF00',
            'T':'#FF0000',
            'G':'#000000',
            'C':'#0000FF'}[n.upper()]

def tm_gc(seq):
    gc = seq.count('G')+seq.count('C')
    at = seq.count('A')+seq.count('T')
    return '4x%s+2x%s=%s'%(gc,at,4*gc+2*at)

def is_cpg(seq,i):
    if i+1>=len(seq):
        return False
    if seq[i]=='C' and seq[i+1]=='G':
        return True
    return False

def is_repeat(seq,i,repeatno=8):
    v = seq[i]
    tail = 0
    for t in seq[i+1:]:
        if t!=v:
            break
        tail += 1
    head = 0
    for t in seq[:i][::-1]:
        if t!=v:
            break
        head += 1
    return tail+head+1 >= repeatno

def count_cpg(seq):
    count = 0
    for i,n in enumerate(seq[:-1]):
        if n=='C' and seq[i+1]=='G':
            count += 1
    return count

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

def search_primer(primer, seq):
    seqs = str(seq)
    cprimer = primer.reverse_complement()
    reg = re.compile('(%s)|(%s)'%(_gen_re_seq(primer),_gen_re_seq(cprimer)))
    pp = []
    pc = []
    start = 0
    while True:
        m = reg.search(seqs, start)
        if not m:
            break
        start = m.start()+1
        if m.group(1):
            pp.append(m.start())
        if m.group(2):
            pc.append(m.start())
    return pp,pc

def bisulfite(seq, methyl):
    muta = seq.tomutable()
    old = 'X'
    for i,c in enumerate(muta[:-2]):
        if c=='C':
            if not(muta[i+1]=='G' and methyl):
                muta[i] = 'T'
    return muta.toseq()

class Annealing:
    def __init__(self, p, q, end_annealing=False):
        s, (i, ss) = Annealing.annealing_score(p,q,end_annealing)
        self.p = p
        self.q = q
        self.score = s
        self.scores = ss
        self.index = i

    @classmethod
    def annealing_score_n(cls,x,y):
        if (x=='A' and y=='T') or (x=='T' and y=='A'):
            return 2
        elif (x=='G' and y=='C') or (x=='C' and y=='G'):
            return 4
        else:
            return 0

    @classmethod
    def annealing_score(cls, p,q,end_annealing=False,getindex=False):
        def sv(x,y):
            return Annealing.annealing_score_n(str(x),str(y))
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

    def get_bar(self):
        sv = Annealing.annealing_score_n

        i = self.index
        p = self.p
        q = self.q
        lp = len(p)
        lq = len(q)
        spc = ' '*abs(i)
        ss = self.scores
        bar = ''.join(['|' if s>0 else ' ' for s in ss])

        if i>0:
            return [spc+"5'-%s-3'"%p, spc+"  <"+bar+">", "3'-%s-5'"%q[::-1] ]
        else:
            return ["5'-%s-3'"%p, spc+"  <"+bar+">", spc+"3'-%s-5'"%q[::-1] ]

    def write_html(self, w):
        w.push('div',style='annealing')
        w.push('p','pea=%s, index=%s'%(pea.score,self.index))
        w.push('pre')
        w.text('\n'.join(self.get_bar()))
        w.pop()
        w.pop()

_nnt_dh = {
    'AA': 9.1,
    'TT': 9.1,
    'AT': 8.6,
    'TA': 6.0,
    'CA': 5.8,
    'TG': 5.8,
    'GT': 6.5,
    'AC': 6.5,
    'CT': 7.8,
    'AG': 7.8,
    'GA': 5.6,
    'TC': 5.6,
    'CG': 11.9,
    'GC': 11.1,
    'GG': 11.0,
    'CC': 11.0,
}
_nnt_dg = {
    'AA': 1.55,
    'TT': 1.55,
    'AT': 1.25,
    'TA': 0.85,
    'CA': 1.15,
    'TG': 1.15,
    'GT': 1.40,
    'AC': 1.40,
    'CT': 1.45,
    'AG': 1.45,
    'GA': 1.15,
    'TC': 1.15,
    'CG': 3.05,
    'GC': 2.70,
    'GG': 2.3,
    'CC': 2.3,
}

class Primer:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
    def __repr__(self):
        return "Primer(%s: %s)"%(self.name,self.seq)
    def __len__(self):
        return len(self.seq)
    def __str__(self):
        return str(self.seq)
    def gc_ratio(self):
        return GC(self.seq)
    def reverse(self):
        return self.seq[::-1]
    def melting_temperature(self, c_na=33.*10**-3, c_mg=(2+4.5)*10**-3, c_primer=0.5*10**-6, unmethyl=True):
        '''
        c_na, c_mg: final concentration in molar of Na and Mg in PCR reaction mix
                    these are used for calculation of salt concentration
        c_primer:   final concentration in molar of each primer in PCR reaction mix
        unmethyl:   if 'unmethyl' is true, all CpGs of template are assumed to be unmethyled
                    then unmethyl version of primer are used for calculation
        
        '''
        if unmethyl:
            seq = Seq.Seq(str(self.seq).replace('R','A').replace('Y','T'),IUPAC.unambiguous_dna)
        else:
            seq = Seq.Seq(str(self.seq).replace('R','G').replace('Y','C'),IUPAC.unambiguous_dna)

        '''
        c_salt: salt concentration in molar, = 
        c_primer: primer concentration in molar
        '''
        c_salt = c_na + 4*(c_mg**0.5)

        l = len(self)
        t0 = 298.2
        r = 1.987
        d_h_e = 5.
        d_h_p = -1000. * ( 2*d_h_e + sum([_nnt_dh[str(seq[i:i+2]).upper()] for i in range(l-1)]) )
        d_g_e = 1.
        d_g_i = -2.2
        d_g_p = -1000. * ( 2*d_g_e + d_g_i + sum([_nnt_dg[str(seq[i:i+2]).upper()] for i in range(l-1)]) )
        t_p = t0*d_h_p / (d_h_p-d_g_p + r*t0*log(c_primer) ) + 16.6*log10(c_salt / (1.+0.7*c_salt)) - 269.3
        return t_p

    def self_annealing(self,getindex=False):
        return Annealing(self.seq, self.seq)
    def self_end_annealing(self,getindex=False):
        return Annealing(self.seq, self.seq, True)

    def debugprint(self):
        print "%s: %s, len=%2d, Tm=%.2f, oTm=%.2f, GC=%.2f, sa=%2d, sea=%2d, MarmurTm: %s" % (self.name, ("5'-%s-3'"%self.seq).ljust(30), len(self), self.melting_temperature(), self.melting_temperature(unmethyl=False), self.gc_ratio(), self.self_annealing().score, self.self_end_annealing().score, tm_gc(self.seq))

        sa = self.self_annealing()
        print 'sa=%s, index=%s'%(sa.score,sa.index)
        print '\n'.join(sa.get_bar())
        sea = self.self_end_annealing()
        print 'sea=%s, index=%s'%(sea.score,sea.index)
        print '\n'.join(sea.get_bar())

class PCRProduct:
    def __init__(self, template, start, end, primer_fw, primer_rv):
        self.template = template
        self.seq = template[start:end]
        self.start = start
        self.start_i = start+len(primer_fw)
        self.end = end
        self.end_i = end-len(primer_rv)
        self.primer_fw = primer_fw
        self.primer_rv = primer_rv
        self.head = template[self.start:self.start_i]
        self.middle = template[self.start_i:self.end_i]
        self.tail = template[self.end_i:self.end]
    def __repr__(self):
        return "PCRProduct(%s -> %s: %s)"%(self.primer_fw.name, self.primer_rv.name, self.seq)
    def __len__(self):
        return len(self.seq)
    def __str__(self):
        return str(self.seq)
    def detectable_cpg(self):
        return len(self.cpg_sites())
    def cpg_sites(self):
        return [i for i in range(self.start_i, self.end_i) if self.template[i]=='C' and self.template[i+1]=='G']

    def write_html(self, w):
        w.push('div')
        cpg = count_cpg(self.seq)
        w.text('length=%d, CpG=%d, detectable CpG=%d'%(len(self), cpg, self.detectable_cpg()))
        w.insertc('br')

        seqstr = self.seq

        s = self.start
        style = [[0,0,0] for i in range(len(seqstr))]
        for i in range(0,self.start_i-s):
            style[i][0]=1
        for i in range(self.end_i-s, self.end-s):
            style[i][1]=1
        j = 0
        while 1:
            j = seqstr.find('CG', j)
            if j<0:
                break
            style[j][2]=1
            style[j+1][2]=1
            j += 1

        colors = [[0,200,0], [0,0,200], [255,0,0]]
        product = lambda m:reduce(operator.mul,m,1)

        w.push('pre')
        i = 0
        length = len(seqstr)
        for i in range(0,length,100):
            w.write('<span style="color:#000">%4d: </span>'%i)
            oldcol = None
            for ii in range(i,min(i+100,length),10):
                for iii in range(ii,min(ii+10,length),1):
                    color = '#%02x%02x%02x'%tuple(sum(map(product,zip(style[iii],col))) for col in zip(*colors))
                    if oldcol!=color:
                        if oldcol:
                            w.write('</span>')
                        w.write('<span style="color:%s">'%color)
                        oldcol = color
                    w.write(seqstr[iii])
                w.write(' ')
            w.write('\n')
            if oldcol:
                w.write('</span>')
        w.pop()

        
        w.push('textarea', cols='10', rows='1', cls='copybox')
        w.write(str(seqstr))
        w.pop()

        w.pop()

class PrimerPair(object):
    def __init__(self, fw, rv):
        self.fw = fw
        self.rv = rv

    def pair_annealing(self):
        return Annealing(self.fw.seq, self.rv.seq)
    def pair_end_annealing(self=False):
        return Annealing(self.fw.seq, self.rv.seq, True)

    def score(self):
        sc_primer_length = [0.5, 23.]
        #sc_gc_content = [1.0, 30.]
        sc_gc_content = [1.0, 50.]
        sc_melting_temperature = [1.0, 60.]
        sc_self_annealing = [0.1, 0.]
        sc_self_end_annealing = [0.2, 0.]
        sc_pair_annealing = [0.1, 0.]
        sc_pair_end_annealing = [0.2, 0.]

        def s(value, parameter):
            return parameter[0]*abs(value-parameter[1])

        score = 0.
        score += s(len(self.fw),sc_primer_length)
        score += s(len(self.rv),sc_primer_length)
        score += s(self.fw.gc_ratio(), sc_gc_content)
        score += s(self.rv.gc_ratio(), sc_gc_content)
        score += s(self.fw.melting_temperature(), sc_melting_temperature)
        score += s(self.rv.melting_temperature(), sc_melting_temperature)
        score += s(self.fw.self_annealing().score, sc_self_annealing)
        score += s(self.fw.self_end_annealing().score, sc_self_end_annealing)
        score += s(self.rv.self_annealing().score, sc_self_annealing)
        score += s(self.rv.self_end_annealing().score, sc_self_end_annealing)
        score += s(self.pair_annealing().score, sc_pair_annealing)
        score += s(self.pair_end_annealing().score, sc_pair_end_annealing)

        return score

    def debugprint(self):
        print 'score=', self.score()
        for r in [self.fw, self.rv]:
            r.debugprint()
        pa = self.pair_annealing()
        print 'pa=%s, index=%s'%(pa.score,pa.index)
        print '\n'.join(pa.get_bar())
        pea = self.pair_end_annealing()
        print 'pea=%s, index=%s'%(pea.score,pea.index)
        print '\n'.join(pea.get_bar())

    def write_html(self, w):
        pa = self.pair_annealing()
        pea = self.pair_end_annealing()
        score = self.score()

        # primer table
        w.push('div',cls='primerpair')

        w.push('table',cls='primerpairtable',border=1)
        pt = ['name', 'sequence', 'length[bp]', 'Tm[C]', 'oTm[C]', 'old Tm[C]','GC[%]', 'sa.', 'sea.', 'pa.', 'pea.', 'pair score']
        w.push('tr')
        for p in pt:
            w.insert('th',p)
        w.pop()
        for r in [self.fw, self.rv]:
            sa = r.self_annealing()
            sea = r.self_end_annealing()
            vt = [r.name, "5'-%s-3'"%r.seq, len(r),
                  '%.2f'%r.melting_temperature(),
                  '%.2f'%r.melting_temperature(unmethyl=False),
                  tm_gc(r.seq),
                  '%.2f'%r.gc_ratio(), 
                  sa.score, sea.score, pa.score, pea.score, 
                  '%.2f'%score]
            w.push('tr')
            for v in vt:
                w.insert('td',str(v))
            w.pop()
        w.pop()

        w.insert('p','Tm melting temperature(SantaLucia), oTm melting temperature for bisulfite methyl-template, sa. self annealing, sea. self end annealing, pa. pair annealing, pea. pair end annealing', style='font-size:x-small')
        w.pop()

class PCR:
    def __init__(self, name, template, primer_fw, primer_rv):
        self.name = name
        self.template = template
        self.primers = PrimerPair(primer_fw,primer_rv)
        self.products = None

    @property
    def fw(self):
        return self.primers.fw
    @property
    def rv(self):
        return self.primers.rv

    def get_products(self):
        if not self.products:
            self.products = self._calc_products()
        return self.products
            
    def _calc_products(self):
        ret = [] 
        fpp,fpc = search_primer(self.primers.fw.seq,self.template)
        rpp,rpc = search_primer(self.primers.rv.seq,self.template)
        f = fpp+rpp
        r = fpc+rpc
        def g(fw,rv,i,j):
            return PCRProduct(self.template,i,j+len(rv),fw, rv)
        for i in fpp:
            for j in fpc:
                ret.append(g(self.primers.fw,self.primers.fw,i,j))
            for j in rpc:
                ret.append(g(self.primers.fw,self.primers.rv,i,j))
        for i in rpp:
            for j in fpc:
                ret.append(g(self.primers.rv,self.primers.fw,i,j))
            for j in rpc:
                ret.append(g(self.primers.rv,self.primers.rv,i,j))
        return ret

    def pair_annealing(self):
        return self.primers.pair_annealing()
    def pair_end_annealing(self=False):
        return self.primers.pair_end_annealing()

    def primer_score(self):
        return self.primers.score()

    def debugprint(self):
        print '%s: score=%.2f'%(self.name, self.primer_score())
        self.primers.debugprint()
        for c in self.get_products():
            print 'product: len=%d, detectable CpG=%d'%(len(c),c.detectable_cpg())
            print c.seq

    def write_html(self, w):
        w.push('div')
        w.insert('h2', self.name)

        self.primers.write_html(w)
        #self.pair_annealing().write_html(w)
        #self.pair_end_annealing().write_html(w)
        for c in self.get_products():
            c.write_html(w)
        w.pop()

def print_html_sequence(w,seq):
    colors = [[0,200,0], [0,0,200], [255,0,0]]
    product = lambda m:reduce(operator.mul,m,1)
    seqstr=str(seq)

    style = [[0,0,0] for i in range(len(seqstr))]
    j = 0
    while 1:
        j = seqstr.find('CG', j)
        if j<0:
            break
        style[j][2]=1
        style[j+1][2]=1
        j += 1

    w.push('pre')
    i = 0
    length = len(seqstr)
    for i in range(0,length,100):
        w.write('<span style="color:#000">%4d: </span>'%i)
        oldcol = None
        for ii in range(i,min(i+100,length),10):
            for iii in range(ii,min(ii+10,length),1):
                color = '#%02x%02x%02x'%tuple(sum(map(product,zip(style[iii],col))) for col in zip(*colors))
                if oldcol!=color:
                    if oldcol:
                        w.write('</span>')
                    w.write('<span style="color:%s">'%color)
                    oldcol = color
                w.write(seqstr[iii])
            w.write(' ')
        w.write('\n')
        if oldcol:
            w.write('</span>')
    w.pop()
    
