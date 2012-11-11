from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

from .xmlwriter import XmlWriter, builder
from .pcr import primers_write_html, primers_write_csv, Primer
from .parser import parse_file
from .prompt import prompt

import sys

__all__ = ['load_primer_list_file']

def load_primer_list_file(fileobj):
    primer = []
    for name, value, em in parse_file(fileobj):
        primer.append(Primer(name, Seq.Seq(value.upper(),IUPAC.ambiguous_dna)))
    return primer

def main():
    import sys, os
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
                primers_write_html(html, primers)

if __name__=='__main__':
    main()
