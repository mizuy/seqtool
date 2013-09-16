from . import to_unambiguous_seq

'''
sense
CG -> YG
anti-sense
CG -> CR


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

    >>> list(cpg_sites('ATGCCGCGATCG'))
    [4, 6, 10]
    >>> list(cpg_sites('ATGCCGCGATCG',(2,7)))
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
    >>> is_cpg('ATGCGC', 3)
    True
    >>> is_cpg('ATGCGC', 4)
    False
    """
    if i+1>=len(seq):
        return False
    if seq[i:i+2] in ['CG','YG','CR']:
        return True
    return False

def count_cpg(seq, range_=(None,None)):
    """
    >>> count_cpg('ATGCYGCGATCG')
    3
    >>> count_cpg('ATGCYGCGATCG'[2:7])
    1
    >>> count_cpg('ATGCYGCGATCG',(2,7))
    2
    """
    p,q = _cpg_range(range_, len(seq))

    return len(re_CpG.findall(str(seq)[p:q]))

def _bisulfite_conversion(seq):
    """
    >>> import Bio.Seq as Seq
    >>> _bisulfite_conversion(Seq.Seq('ATGCCGATGC'))
    Seq('ATGTYGATGT', Alphabet())
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

    return muta.toseq()

def bisulfite_conversion(seq, sense=True):
    """
    >>> import Bio.Seq as Seq
    >>> bisulfite_conversion(Seq.Seq('ATGCGC'), sense=True)
    Seq('ATGYGT', Alphabet())
    >>> bisulfite_conversion(Seq.Seq('ATGCGC'), sense=False)
    Seq('ATACRC', Alphabet())
    """
    return asymmetric_conversion(seq, _bisulfite_conversion, sense=sense)

def asymmetric_conversion(seq, conv, sense):
    if sense:
        return conv(seq)
    else:
        return conv(seq.reverse_complement()).reverse_complement()

def to_unambiguous(bsseq, methyl=True):
    """
    >>> str(to_unambiguous('ATGTYG',True))
    'ATGTCG'
    >>> str(to_unambiguous('ATGTYG',False))
    'ATGTTG'
    """
    trans = str.maketrans('YR','CG') if methyl else str.maketrans('YR','TA')
    return to_unambiguous_seq(str(bsseq).translate(trans))

def bisulfite_conversion_unambiguous(seq, sense, methyl):
    return to_unambiguous(bisulfite_conversion(seq, sense=sense), methyl=methyl)

def bisulfite(seq, methyl, sense=True):
    """
    >>> import Bio.Seq as Seq
    >>> str(bisulfite(Seq.Seq('ATGCGC'), methyl=True))
    'ATGCGT'
    >>> str(bisulfite(Seq.Seq('ATGCGC'), methyl=False))
    'ATGTGT'
    """
    key = '_bisulfite_' + ('met' if methyl else 'unmet') + '_' + ('sense' if sense else 'asense')

    if not hasattr(seq, key):
        ret = bisulfite_conversion_unambiguous(seq, sense=sense, methyl=methyl)
        setattr(seq, key, ret)

    return getattr(seq, key)

def gc_ratio(seq):
    """
    >>> gc_ratio('ATGCCGGATT')
    50.0
    """
    c = seq.count('C')
    c += seq.count('G')
    c += seq.count('Y')
    c += seq.count('R')
    return 100. * c/len(seq)

def _c2t_conversion(seq):
    muta = seq.tomutable()
    for i,c in enumerate(muta):
        if c=='C':
            muta[i] = 'Y'
    return muta.toseq()

def c2t_conversion(seq, sense=True):
    return asymmetric_conversion(seq, lambda x: _c2t_conversion(x), sense)

class BisulfiteTemplate:
    def __init__(self, origin_seq):
        self.origin = origin_seq
        self.sense = bisulfite_conversion(True)
        self.asense = bisulfite_conversion(False)
