# encoding: utf-8


from Bio import Seq
from Bio.Alphabet import IUPAC

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(prog='translate', description='Print translation of nucleotide')
    parser.add_argument('sequence', nargs=1, help="sequence from 5' to 3'")
    args = parser.parse_args()

    s = Seq.Seq(args.sequence[0],IUPAC.unambiguous_dna)
    print('+0',s.translate())
    print('+1',s[1:].translate())
    print('+2',s[2:].translate())

    s = s.reverse_complement()
    print('+0',s.translate())
    print('+1',s[1:].translate())
    print('+2',s[2:].translate())


if __name__=='__main__':
    main()
