# distutils: extra_compile_args = -w

import numpy as np
cimport numpy as np
#np.import_array()

# http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
# http://www.ibm.com/developerworks/jp/java/library/j-seqalign/


'''
no optomization
bin/py script/pcmv.py GGAGCTCCTGATTTAGAGCTTGACGGGGAAAG  14.25s user 0.15s system 94% cpu 15.215 total

little optimization
bin/py script/pcmv.py GGAGCTCCTGATTTAGAGCTTGACGGGGAAAG  0.26s user 0.11s system 98% cpu 0.373 total
'''
cdef int SCORE_MATCH = 2
cdef int SCORE_MISMATCH = -2
cdef int SCORE_GAP = -1

cdef inline int match(int i,int j):
    if i == j:
        return SCORE_MATCH
    else:
        return SCORE_MISMATCH

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

    def location(self):
        return (len(self.first), len(self.first)+len(self.mid))

cdef _sw(bytes s0_, bytes s1_):
    cdef char* s0 = s0_
    cdef char* s1 = s1_
    cdef int m = len(s0_)
    cdef int n = len(s1_)

    cdef np.ndarray[np.int_t, ndim=2] matrix = np.zeros((m+1,n+1),dtype=np.int)
    cdef np.ndarray[np.int8_t, ndim=2] direction = np.zeros((m+1,n+1),dtype=np.int8)

    cdef int j, i, w, p, q, ii, jj, a0, a1
    cdef int diag, up, down, mv, md, v

    for j in range(1,n+1):
        for i in range(1,m+1):
            a0 = <int>s0[i-1]
            a1 = <int>s1[j-1]
            diag = matrix[i-1,j-1] + match(a0,a1)
            down = matrix[i-1, j] + SCORE_GAP
            up = matrix[i, j-1] + SCORE_GAP

            mv = 0
            md = 0

            if diag > mv:
                mv = diag
                md = 1
            if down > mv:
                mv = down
                md = 2
            if up > mv:
                mv = up
                md = 3

            matrix[i,j] = mv
            direction[i,j] = md

    # TODO: list all maximum
    #p,q = np.unravel_index(matrix.argmax(), (matrix.shape[0],matrix.shape[1]))
    mv = 0
    p = 0
    q = 0
    for j in range(1,n+1):
        for i in range(1,m+1):
            v = matrix[i,j]
            if v>mv:
                mv = v
                p = i
                q = j
    tracer = [(p,q,direction[p,q])]

    while True:
        i,j,d = tracer[-1]
        if matrix[i,j]==0:
            break
        ii,jj = [(0,0), (i-1,j-1), (i-1,j), (i,j-1)][d]
        tracer.append((ii,jj,direction[ii,jj]))

    return tracer, mv

def smith_waterman(str s0, str s1):
    cdef bytes s0__ = s0.encode()
    cdef bytes s1__ = s1.encode()
    tracer,mv = _sw(s0__, s1__)

    p,q,d = tracer[0]

    mid0 = ''
    mid1 = ''
    for p,q,d in tracer:
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

    first = tracer[-1][:2]
    last = tracer[0][:2]

    aseq0 = AlingnedSeq(s0, s0[:first[0]], mid0, s0[last[0]:])
    aseq1 = AlingnedSeq(s1, s1[:first[1]], mid1, s1[last[1]:])

    return aseq0, aseq1, mv
