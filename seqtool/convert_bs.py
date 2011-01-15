#!/usr/bin/env python 
import commands
from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

def bisulfite(seq, methyl):
    muta = seq.tomutable()
    for i in xrange(len(muta)-1):
        if muta[i].upper()=='C' and not (methyl and muta[i+1].upper()=='G'):
            muta[i] = 'T'
    if muta[-1].upper()=='C':
        muta[-1] = 'T'
    return muta

def reverse(seq):
    return seq.reverse_complement()

def convert_fasta(records, output_file, conv, name_addition, description_addition):
    for record in records:
        rid = record.id+name_addition
        desc = record.description + description_addition
        SeqIO.write([SeqRecord(id=rid, seq=conv(record.seq), description=desc)], output_file, 'fasta')

def main():
    import sys, os
    from optparse import OptionParser

    parser = OptionParser('usage: %prog [options] fasta.fa fasta2.fa ...')

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.error("no input file")
    inputfiles = args
    
    for i in inputfiles:
        print i
        records = list(rec.upper() for rec in SeqIO.parse(open(i,'r'), "fasta"))

        b = os.path.splitext(i)[0]

        print 'met+'
        out = open(b+'_bs_sense_met.fa', 'w')
        convert_fasta(records, out, lambda x:bisulfite(x, True), 'bs_sense_met', '; with bisulfite converseion (sense methyl)')

        print 'unmet+'
        out = open(b+'_bs_sense_unmet.fa', 'w')
        convert_fasta(records, out, lambda x:bisulfite(x, False), 'bs_sense_unmet', '; with bisulfite converseion (sense unmethyl)')

        print 'met-'
        out = open(b+'_bs_antisense_met.fa', 'w')
        convert_fasta(records, out, lambda x:bisulfite(reverse(x), True), 'bs_antisense_met', '; with bisulfite converseion (antisense methyl)')

        print 'unmet-'
        out = open(b+'_bs_antisense_unmet.fa', 'w')
        convert_fasta(records, out, lambda x:bisulfite(reverse(x), False), 'bs_antisense_unmet', '; with bisulfite converseion (antisense unmethyl)')
