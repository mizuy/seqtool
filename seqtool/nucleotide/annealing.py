__all__ = ['OligoAnnealing']

def annealing_score_n(x,y):
    if (x=='A' and y=='T') or (x=='T' and y=='A'):
        return 2
    elif (x=='G' and y=='C') or (x=='C' and y=='G'):
        return 4
    else:
        return 0

class AnnealingState:
    def __init__(self, p, q, score, index, scores):
        self.p = p
        self.q = q
        self.score = score
        self.scores = scores
        self.index = index

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

def annealing_score(p,q):
    """
    >>> fw = 'GAAGGAGACCCAAATTCAAAGTT'
    >>> rv = 'CCTTTCTCCCTTCGTAGGT'
    >>> annealing_score(fw, rv)
    ((18, -7, [2, 4, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 4, 0, 0]), (10, -7, [2, 4, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 4, 0, 0]))
    >>> annealing_score(fw, fw)
    ((24, 6, [4, 2, 2, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 2, 2, 4]), (4, -18, [2, 2, 0, 2, 2]))
    """
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

    ea_v = (-1, None, None)
    def u_ea(score, index, scores):
        nonlocal ea_v
        os,oi,oss = ea_v
        if os < score:
            ea_v = (score, index, scores)

    a_v = (-1, None, None)
    def u_a(score, index, scores):
        nonlocal a_v
        os,oi,oss = a_v
        if os < score:
            a_v = (score, index, scores)
            
    if n<=m:
        assert m-n >= 0
        for k in range(-(n-1),m-1 +1):
            if k<=0:
                # 5'- w[0]....w[-k]....w[n-1] -3'
                #         3'- v[0].....v[n+k-1]....v[m-1] -5'
                ss = [sv(w[-k+i],v[i]) for i in range(n+k)]
                u_a(sum(ss),k,ss)
                u_ea(ea_lr(ss),k,ss)
            elif k<=m-n:
                #         w[0]....w[n-1]
                # v[0]....v[k]....v[k+n-1].....v[m-1]
                ss = [sv(w[0+i],v[k+i]) for i in range(n)]
                u_a(sum(ss),k,ss)
                u_ea(ea_r(ss),k,ss)
            else:
                #        w[0]...w[m-k-1]....w[n-1]
                # v[0]...v[k]...v[m-1]
                ss = [sv(w[i],v[k+i]) for i in range(m-k)]
                u_a(sum(ss),k,ss)
    else:
        assert m-n <= 0
        for k in range(-(n-1),m-1 +1):
            if k<=m-n:
                # w[0]....w[-k]....w[n-1]
                #         v[0].....v[n+k-1]....v[m-1]
                ss = [sv(w[-k+i],v[i]) for i in range(n+k)]
                u_a(sum(ss),k,ss)
                u_ea(ea_lr(ss),k,ss)
            elif k<=0:
                # w[0]....w[k]....w[m-k-1].....w[n-1]
                #         v[0]....v[m-1]
                ss = [sv(w[k+i],v[0+i]) for i in range(m)]
                u_a(sum(ss),k,ss)
                u_ea(ea_l(ss),k,ss)
            else:
                #        w[0]...w[m-k-1]....w[n-1]
                # v[0]...v[k]...v[m-1]
                ss = [sv(w[i],v[k+i]) for i in range(m-k)]
                u_a(sum(ss),k,ss)

    return a_v, ea_v

class OligoAnnealing:
    def __init__(self, p, q):
        a_v, ea_v = annealing_score(p,q)
        self.max_annealing = AnnealingState(p, q, a_v[0], a_v[1], a_v[2])
        self.end_annealing = AnnealingState(p, q, ea_v[0], ea_v[1], ea_v[2])
        
