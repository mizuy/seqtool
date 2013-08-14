
import sys

from ..nucleotide.primer import Primer
from ..nucleotide import to_seq, primer_cond

def main():
    pc = primer_cond.PrimerCondition()
    pc.tm = primer_cond.Condition(1.0, 60., 59., 61.)
    pc.primer_length = primer_cond.Condition(.5, 23., 17, 30)
    pc.sa = primer_cond.Condition(0.1, 0., 0, 100)
    pc.sea = primer_cond.Condition(0.1, 0., 0, 50)
    pc.gc = primer_cond.Condition(1., 50., 30, 90)

    a,b = pc.primer_length.minimum, pc.primer_length.maximum

    if len(sys.argv)!=2:
        print('usage: pdesign [sequence]')
        return
    s = to_seq(sys.argv[1])
    l = len(s)

    for i in range(l):
        for j in range(i+a, i+b):
            name = "p{}-{}".format(i,j)
            pp = Primer(name, s[i:j])
            if not pc.bound_primer(pp):
                continue
            print("{}: {}".format(pp.name, pp.seq))

if __name__=='__main__':
    main()
