from . import to_unambiguous_seq, to_ambiguous_seq

'''
sense
CG -> YG
C* -> T*

anti-sense
CG -> CR
*G -> *A

met-sense
CG -> CG
C* -> T*
# this==C and next==G -> C(this)
  this==C and next!=G -> T
  this!=C -> this

unmet-sense
CG -> TG
C* -> T*
# just convert every C to T

met-asense
CG -> CG
*G -> *A
# this=G and last==C -> G(this)
  this=G and last!=C -> A
  this!=G -> this

unmet-sense
CG -> CA
*G -> *A
# just convert every G to A

'''

import re
re_CpG = re.compile('(CG)|(YG)|(CR)')

def _cpg_range(range_, length):
    p,q = range_
    p = p or 0
    if q:
        q = min(q+1, length)
    else:
        q = length
    return p,q

def cpg_sites(seq, range_=(None,None)):
    """
    return all cpg location of seq

    >>> seq = to_unambiguous_seq('ATGCCGCGATCG')
    >>> list(cpg_sites(seq))
    [4, 6, 10]
    >>> list(cpg_sites(seq,(2,7)))
    [4, 6]
    """
    seqstr = str(seq)

    p,q = _cpg_range(range_, len(seq))

    seqstr = seqstr

    yield from [p+m.start() for m in re_CpG.finditer(seqstr[p:q])]
    '''
    j = p
    while 1:
        j = seqstr.find('CG', j)
        if j<0 or q<=j:
            break
        yield j
        j += 1
    '''

def is_cpg(seq,i):
    """
    >>> seq = to_unambiguous_seq('ATGCGC')
    >>> is_cpg(seq, 3)
    True
    >>> is_cpg(seq, 4)
    False
    """
    if not (i+1<len(seq)):
        raise IndexError('out of range')
    if str(seq)[i:i+2] in ['CG','YG','CR']:
        return True
    return False

def count_cpg(seq, range_=(None,None)):
    """
    >>> seq = to_ambiguous_seq('ATGCYGCGATCG')
    >>> count_cpg(seq)
    3
    >>> count_cpg(seq[2:7])
    1
    >>> count_cpg(seq,(2,7))
    2
    """
    p,q = _cpg_range(range_, len(seq))

    return len(re_CpG.findall(str(seq)[p:q]))

def bisulfite_conversion(seq, sense=True):
    """
    >>> seq = to_unambiguous_seq('ATGCGC')
    >>> bisulfite_conversion(seq, sense=True)
    Seq('ATGYGT', IUPACAmbiguousDNA())
    >>> bisulfite_conversion(seq, sense=False)
    Seq('ATACRC', IUPACAmbiguousDNA())
    """
    muta = str(seq)
    if sense:
        muta = muta.replace('CG','YG').replace('C','T')
    else:
        muta = muta.replace('CG','CR').replace('G','A')
    return to_ambiguous_seq(muta)
    
        
def to_unambiguous(bsseq, methyl=True):
    """
    >>> to_unambiguous('ATGTYG', True)
    Seq('ATGTCG', IUPACUnambiguousDNA())
    >>> to_unambiguous('ATGTYG', False)
    Seq('ATGTTG', IUPACUnambiguousDNA())
    """
    trans = str.maketrans('YR','CG') if methyl else str.maketrans('YR','TA')
    return to_unambiguous_seq(str(bsseq).translate(trans))

def bisulfite_conversion_unambiguous(seq, sense, methyl):
    return to_unambiguous(bisulfite_conversion(seq, sense=sense), methyl=methyl)

def bisulfite(seq, methyl, sense=True):
    """
    >>> import Bio.Seq as Seq
    >>> bisulfite(Seq.Seq('ATGCGC'), methyl=True)
    Seq('ATGCGT', IUPACUnambiguousDNA())
    >>> bisulfite(Seq.Seq('ATGCGC'), methyl=False)
    Seq('ATGTGT', IUPACUnambiguousDNA())
    >>> bisulfite(Seq.Seq('ATGCGC'), methyl=True, sense=False)
    Seq('ATACGC', IUPACUnambiguousDNA())
    >>> bisulfite(Seq.Seq('ATGCGC'), methyl=False, sense=False)
    Seq('ATACAC', IUPACUnambiguousDNA())
    """
    key = '_b' + ('m' if methyl else 'u') + ('p' if sense else 'n')

    if not hasattr(seq, key):
        setattr(seq, key, bisulfite_conversion_unambiguous(seq, sense=sense, methyl=methyl))
    return getattr(seq, key)

def gc_ratio(seq):
    """
    >>> gc_ratio('ATGCCGGATT')
    50.0
    >>> gc_ratio('ATGYYGGATT')
    40.0
    """
    c = seq.count('C')
    c += seq.count('G')
    c += .5 * seq.count('Y')
    c += .5 * seq.count('R')
    return 100. * c/len(seq)

class BisulfiteTemplate:
    def __init__(self, origin_seq):
        self.origin = origin_seq
        self.sense = bisulfite_conversion(True)
        self.asense = bisulfite_conversion(False)


# C->T conversion

def _c2t_conversion(seq):
    muta = seq.tomutable()
    for i,c in enumerate(muta):
        if c=='C':
            muta[i] = 'Y'
    return muta.toseq()

def c2t_conversion(seq, sense=True):
    return asymmetric_conversion(seq, lambda x: _c2t_conversion(x), sense)
    
# reference

def _bisulfite_conversion(seq):
    """
    >>> seq = to_unambiguous_seq('ATGCCGATGC')
    >>> _bisulfite_conversion(seq)
    Seq('ATGTYGATGT', IUPACAmbiguousDNA())
    """
    seqstr = str(seq)
    l = len(seqstr)
    muta = seq.tomutable()

    j = 0
    while j<l:
        j = seqstr.find('C', j)

        if not (0<=j and j<l):
            break
        if j + 1<l and seqstr[j+1]=='G':
            muta[j] = 'Y'
        else:
            muta[j] = 'T'

        j += 1

    return to_ambiguous_seq(str(muta))
    
def bisulfite_conversion_slow(seq, sense=True):
    """
    >>> seq = to_unambiguous_seq('ATGCGC')
    >>> bisulfite_conversion_slow(seq, sense=True)
    Seq('ATGYGT', IUPACAmbiguousDNA())
    >>> bisulfite_conversion_slow(seq, sense=False)
    Seq('ATACRC', IUPACAmbiguousDNA())
    """
    return asymmetric_conversion(seq, _bisulfite_conversion, sense=sense)

def asymmetric_conversion(seq, conv, sense):
    if sense:
        return conv(seq)
    else:
        return conv(seq.reverse_complement()).reverse_complement()


