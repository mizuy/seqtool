import numpy

# http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
# http://www.ibm.com/developerworks/jp/java/library/j-seqalign/

class Score(object):
    def __init__(self, match, mismatch, gap):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

SCORE = Score(2,-2,-1)

def match_nucleotide(a, b):
    if a==b:
        return True
    else:
        return False

def match_bar(s0, s1):
    assert(len(s0)==len(s1))
    return ''.join('|' if match_nucleotide(s0[i],s1[i]) else ' ' for i in range(len(s0)))

class AlingnedSeq(object):
    def __init__(self, seq, first, mid, last):
        self.seq = seq
        self.first = first
        self.mid = mid
        self.last = last

        self.adjust = len(self.first)

    def combined(self):
        return self.first + self.mid + self.last

    def local(self, upstream, downstream):
        return self.first[-upstream:].rjust(upstream) + self.mid + self.last[:downstream].ljust(downstream)



def smith_waterman(s0, s1):
    m = len(s0)
    n = len(s1)

    matrix = numpy.zeros((m+1,n+1),dtype=numpy.int32)
    direction = numpy.zeros((m+1,n+1),dtype=numpy.int32)

    def match(i,j):
        if match_nucleotide(s0[i-1],s1[j-1]):
            return SCORE.match
        else:
            return SCORE.mismatch

    for j in range(1,n+1):
        for i in range(1,m+1):
            v = numpy.array([0, matrix[i-1,j-1] + match(i,j), 
                               matrix[i-1, j] + SCORE.gap,
                               matrix[i, j-1] + SCORE.gap])
            # TODO: list all maximum
            w = v.argmax()
            matrix[i,j] = v[w]
            direction[i,j] = w

    def backward(i,j):
        d = direction[i,j]
        return [(0,0), (i-1,j-1), (i-1,j), (i,j-1)][d]

    # TODO: list all maximum
    p,q = numpy.unravel_index(matrix.argmax(), matrix.shape)
    tracer = [(p,q)]
    while True:
        i,j = tracer[-1]
        if matrix[i,j]==0:
            break
        ii,jj = backward(i,j)
        tracer.append((ii,jj))

    p,q = tracer[0]

    mid0 = ''
    mid1 = ''
    for p,q in tracer:
        d = direction[p,q]
        if d==0:
            break
        elif d==1:
            mid0 = s0[p-1] + mid0
            mid1 = s1[q-1] + mid1
        elif d==2:
            mid0 = s0[p-1] + mid0
            mid1 = '-' + mid1
        elif d==3:
            mid0 = '-' + mid0
            mid1 = s1[q-1] + mid1

    first = tracer[-1]
    last = tracer[0]

    aseq0 = AlingnedSeq(s0, s0[:first[0]], mid0, s0[last[0]:])
    aseq1 = AlingnedSeq(s1, s1[:first[1]], mid1, s1[last[1]:])

    return aseq0, aseq1

class Alignment(object):
    def __init__(self, seq_s0, seq_s1):
        self.aseq0, self.aseq1 = smith_waterman(seq_s0,seq_s1)

    def text_all(self):
        p,q = self.aseq0.adjust, self.aseq1.adjust
        r = max(p,q)
        return '{}\n{}'.format(' '*(r-p) + self.aseq0.combined(), 
                               ' '*(r-q) + self.aseq1.combined() )

    def text_local(self, upstream=30, downstream=30):
        p,q = self.aseq0.adjust, self.aseq1.adjust
        return '{}\n{}\n{}'.format(self.aseq0.local(upstream,downstream), 
                               ' '*upstream+match_bar(self.aseq0.mid, self.aseq1.mid),
                               self.aseq1.local(upstream,downstream) )

def print_sw(s0, s1):
    a = Alignment(s0, s1)

    print 'alignment of :'
    print a.aseq0.seq
    print a.aseq1.seq
    print 'is:'
    print a.text()

if __name__=='__main__':
    print_sw('ACACACTA','AGCACACA')
    print_sw('GCCCTAGCG', 'GCGCAATG')
