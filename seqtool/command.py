from __future__ import absolute_import

import sys, os
from argparse import ArgumentParser

from .dirutils import Filepath
from . import seqview as sv
from . import tssview as tv
from . import primers
from .db import entrez

def seqview():
    parser = ArgumentParser(prog='seqview', description='pretty HTML+SVG report of sequence and adittional data')
    parser.add_argument("seqvfile", help=".seqv file")
    parser.add_argument("-o", "--output", dest="output", help="output filename")
    
    args = parser.parse_args()

    inputp = Filepath(args.seqvfile)
    if args.output:
        outputp = Filepath(args.output)
    else:
        outputp = inputp.change_ext('.html')

    p = sv.SeqvFile()
    p.load_seqvfileentry(inputp.path)
    p.write_html(outputp)

def geneview():
    parser = ArgumentParser(prog='geneview', description='pretty HTML+SVG report of gene')
    parser.add_argument("gene_symbol", help="Gene Symbol")
    parser.add_argument("-o", "--output", dest="output", help="output filename")
    
    args = parser.parse_args()

    if not args.output:
        parser.error('no output')
        return

    outputp = Filepath(args.output)
    gene_id, gene_symbol = entrez.get_gene_from_text(args.gene_symbol)

    e = sv.GeneBankEntry(gene_symbol)
    e.load_genbank(entrez.get_genomic_context_genbank(gene_id))

    p = sv.SeqvFile()
    p.load_genbankentry(e)
    p.write_html(outputp)


def tssview():
    parser = ArgumentParser(prog='tssview', description='pretty HTML+SVG report of dbtss of multiple genes')
    parser.add_argument("tssvfile", help=".tssv file")
    parser.add_argument("-o", "--output", dest="output", help="output filename")
    
    args = parser.parse_args()

    inputp = Filepath(args.tssvfile)
    if args.output:
        outputp = Filepath(args.output)
    else:
        outputp = inputp.change_ext('.html')

    with open(inputp.path,'r') as f:
        tssv = tv.TssvFile(f)

    tssv.write_csv(outputp.change_ext('.csv').path)
    tssv.write_html(outputp)

def get_genbank():
    parser = ArgumentParser(prog='get_genbank', description='Retrieve Genbank from Gene ID or Gene Symbol')
    parser.add_argument("gene", help="Gene ID or Gene Symbol")
    parser.add_argument("-o", "--output", dest="output", help="output filename")

    args = parser.parse_args()

    if args.output:
        output = open(outputfile,'w')
    else:
        output = sys.stdout

    try:
        genbank = entrez.get_genomic_context_genbank(args.gene)
    except Exception,e:
        print str(e)
        return

    if not genbank:
        print "GenBank retrieve error: ", gene_id
        return
    print >>output, genbank

def primers():
    from optparse import OptionParser

    parser = OptionParser('usage: %prog inputfile -o outputfile')
    parser.add_option('-o', '--output', dest='output', help='output filename')
    parser.add_option('-c', '--csv', action='store_true', dest='csv', default='False', help='write a csv file(default is html file)')

    if len(sys.argv) == 0:
        parser.error('no input file')

    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        return

    inputfiles = args
    
    primers = []
    with prompt('loading primers: ') as pr:
        for filename in inputfiles:
            with open(filename,'r') as f:
                for i in load_primer_list_file(f):
                    primers.append(i)
                    pr.progress()
    
    if options.output:
        output = open(options.output, 'w')
    else:
        output = sys.stdout

    if options.csv:
        output.write(primers_write_csv(primers))
    else:
        print 'writing html...'
        html = XmlWriter(output)
        b = builder(html)
        with b.html:
            with b.head:
                pass
            with b.body(style='font-family:monospace;font-size:small'):
                b.h1('Primers')
                primers.primers_write_html(html, primers)

