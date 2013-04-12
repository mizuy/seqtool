#!bin/py

import sys, os
import readline

REGEX = False

from seqtool.nucleotide import to_seq
from seqtool.nucleotide.cpg import bisulfite
from seqtool.nucleotide import sw
from Bio import SeqIO

def load_fasta(filename):
    seqs = []

    fff = open(filename,'r')
    for record in SeqIO.parse(fff, "fasta"):
        s = record.seq
        seqs.append((record.id + ' origin(MH)', s))
        seqs.append((record.id + ' BS+ met', bisulfite(s, True, True)))
        seqs.append((record.id + ' BS+ unm(N)', bisulfite(s, False, True)))
        seqs.append((record.id + ' BS- met', bisulfite(s, True, False)))
        seqs.append((record.id + ' BS- unm(N)', bisulfite(s, False, False)))

    return seqs

def print_result(target, seqs):
    print 'target length = ',len(target)
    for (name,s) in seqs:
        print name

        sense = sw.Alignment(str(target), str(s))
        asense = sw.Alignment(str(to_seq(target).reverse_complement()), str(s))
        print '  sense: target {}, template {}'.format(sense.aseq0.location(), sense.aseq1.location())
        print sense.text_local()
        print '  asense: target {}, template {}'.format(asense.aseq1.location(), asense.aseq1.location())
        print asense.text_local()

def clear_screen():
    #os.system('clear')
    #sys.stderr.write("\x1b[2J\x1b[H")
    print chr(27) + "[2J"

def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(prog='align', description='align sequence result to sequences in fasta')
    parser.add_argument('fasta', nargs=1, help="fasta")
    args = parser.parse_args()

    seqs = load_fasta(args.fasta[0])

    while True:
        s = raw_input('> ')
        if not s:
            continue

        try:
            target = s.upper()
            print_result(target, seqs)
        except:
            print 'error. try again.'
            continue


if __name__=='__main__':
    main()

def c2t(seq):
    muta = seq.tomutable()
    for i,c in enumerate(muta):
        if c=='C':
            muta[i] = 'T'
    return muta.toseq()
