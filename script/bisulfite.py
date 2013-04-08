# encoding: utf-8
from __future__ import absolute_import

from Bio import Seq
from Bio.Alphabet import IUPAC

from seqtool import nucleotide as nuc


def reverse(seq):
    return str(seq)[::-1]

def negative(seq):
    return "3'-{0}-5'".format(reverse(seq))

def positive(seq):
    return "5'-{0}-3'".format(seq)

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(prog='bisulfite', description='Print bisulfite-modified sequence')
    parser.add_argument('sequence', nargs=1, help="sequence from 5' to 3'")
    args = parser.parse_args()

    s = Seq.Seq(args.sequence[0],IUPAC.unambiguous_dna)

    sense = nuc.bisulfite(s, True)
    asense = nuc.bisulfite(s.reverse_complement(), True)

    print 'Fw = ',positive(s)
    print 'Rv = ',positive(s.reverse_complement())
    print ""

    print negative(sense.reverse_complement())
    print positive(sense)
    print "↑"
    print positive(s)
    print negative(s.reverse_complement())
    print "↓"
    print negative(asense)
    print positive(asense.reverse_complement())


if __name__=='__main__':
    main()
