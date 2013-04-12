import sys, os

REGEX = False

from seqtool.nucleotide import to_seq
from seqtool.nucleotide.pcr import Primer
from seqtool.nucleotide.cpg import bisulfite
from seqtool.nucleotide import sw
from Bio import SeqIO

seqs = []


fff = open(os.path.join(os.path.dirname(__file__), 'test_sequence.txt') ,'r')
for record in SeqIO.parse(fff, "fasta"):
    s = record.seq
    seqs.append((record.id + ' origin(MH)', s))
    seqs.append((record.id + ' BS+ met', bisulfite(s, True, True)))
    seqs.append((record.id + ' BS+ unm(N)', bisulfite(s, False, True)))
    seqs.append((record.id + ' BS- met', bisulfite(s, True, False)))
    seqs.append((record.id + ' BS- unm(N)', bisulfite(s, False, False)))


def preceeding(seq, loc, width):
    return seq[max(0, loc-width):loc]
def subsequent(seq, loc, width):
    return seq[loc:min(len(seq),loc+width)]

def main():
    if len(sys.argv) < 1:
        print 'no sequence'
        return
    #os.system('clear')
    #sys.stderr.write("\x1b[2J\x1b[H")
    print chr(27) + "[2J"
    target = sys.argv[1].upper()
    tp = Primer('target',target)
    for (name,s) in seqs:
        print name

        if REGEX:
            def fo(left, right):
                return '{}-{}-{}'.format(preceeding(s, left, 20),s[left:right], subsequent(s, right, 20))

            pp, pc = tp.search(s, None)
            for p in pp:
                print ' ','sense', p.loc_3p, fo(p.left, p.right)
            for p in pc:
                print ' ','asense', p.loc_3p, fo(p.left, p.right)
        else:
            sense = sw.Alignment(str(target), str(s))
            asense = sw.Alignment(str(to_seq(target).reverse_complement()), str(s))
            print '  sense:\n',sense.text_local()
            print '  asense:\n',asense.text_local()


if __name__=='__main__':
    main()

def c2t(seq):
    muta = seq.tomutable()
    for i,c in enumerate(muta):
        if c=='C':
            muta[i] = 'T'
    return muta.toseq()
