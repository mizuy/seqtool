#!bin/py

import sys, os
import readline
import traceback

from seqtool.nucleotide import to_seq
from seqtool.nucleotide.cpg import bisulfite_conversion_ambiguous, c2t_conversion
from seqtool.nucleotide import sw
from Bio import SeqIO

conversions = [
    ('origin(MH)', lambda x: x),
    ('BS+ (N)', lambda x: bisulfite_conversion_ambiguous(x, True)),
    ('BS- (N)', lambda x: bisulfite_conversion_ambiguous(x, False)),
    ('C2T+ (NMH)', lambda x: c2t_conversion(x, True)),
    ('C2T- (NMH)', lambda x: c2t_conversion(x, False)) ]

def foreach_conversions(seq):
    for cname, func in conversions:
        for sense in [True,False]:
            name = '{} {}'.format(cname,'sense' if sense else 'asens')
            yield name, str(func(seq)), sense

def load_fasta(filename):
    seqs = []

    fff = open(filename,'r')
    for record in SeqIO.parse(fff, "fasta"):
        s = record.seq
        seqs.append((record.id, s))

    return seqs


def print_result(target, seqs):
    print 'target length = ',len(target)

    for (name,s) in seqs:
        print name

        results = []

        for cname, cseq, sense in foreach_conversions(s):
            if sense:
                mtarget = str(target)
            else:
                mtarget = str(to_seq(target).reverse_complement())

            al = sw.Alignment(mtarget, cseq)

            pc =  '     {}: score {}'.format(cname, al.score_text())

            loc = '     target {}, template {}'.format(
                                        al.aseq0.location(),al.aseq1.location())
            pp = '{}\n{}\n{}\n\n'.format(pc, loc, al.text_local())

            print pc
            results.append((al.score, pp))

        mscore = 0
        mpp = []
        for score, pp in results:
            if score > mscore:
                mscore = score
                mpp = [pp]
            elif score == mscore:
                mpp.append(pp)
        print ' maximum score: ',mscore
        for m in mpp:
            print m

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
            print traceback.format_exc()
#            print sys.exc_info()[0]
#            print 'error. try again.'
            continue


if __name__=='__main__':
    main()

def c2t(seq):
    muta = seq.tomutable()
    for i,c in enumerate(muta):
        if c=='C':
            muta[i] = 'T'
    return muta.toseq()
