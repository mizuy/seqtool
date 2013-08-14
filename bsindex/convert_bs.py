#!/usr/bin/env python 
import subprocess

from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

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
        print(i)
        records = list(rec.upper() for rec in SeqIO.parse(open(i,'r'), "fasta"))

        b = os.path.splitext(i)[0]

        print('met+')
        out = open(b+'_bs_pos_met.fa', 'w')
        convert_fasta(records, out, lambda x:bisulfite(x, True), 'bs_pos_met', '; with bisulfite converseion (pos methyl)')

        print('unm+')
        out = open(b+'_bs_pos_unm.fa', 'w')
        convert_fasta(records, out, lambda x:bisulfite(x, False), 'bs_pos_unm', '; with bisulfite converseion (pos unmethyl)')

        print('met-')
        out = open(b+'_bs_neg_met.fa', 'w')
        convert_fasta(records, out, lambda x:bisulfite(reverse(x), True), 'bs_neg_met', '; with bisulfite converseion (neg methyl)')

        print('unm-')
        out = open(b+'_bs_neg_unm.fa', 'w')
        convert_fasta(records, out, lambda x:bisulfite(reverse(x), False), 'bs_neg_unm', '; with bisulfite converseion (neg unmethyl)')
