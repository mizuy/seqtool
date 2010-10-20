from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

from .xmlwriter import XmlWriter, builder
from .nucleotide import tm_gc, Primer
from .parser import parse_file

__all__ = ['load_primer_list_file','primers_write_html']

def load_primer_list_file(fileobj):
    primer = []
    for name, value, em in parse_file(fileobj):
        primer.append(Primer(name, Seq.Seq(value.upper(),IUPAC.ambiguous_dna)))
    return primer

def primers_write_html(html, primers):
    pt = ['name', 'sequence', 'len[bp]', 'Tm[C]', 'oTm[C]', 'old Tm[C]','GC[%]', 'sa', 'sea']
    width = [16, 26, None, None, None, None, None, None, None, None]

    b = builder(html)
    with b.table(cls='primer',border=1):
        with b.tr(style='max-height:1em'):
            for p in pt:
                b.th(p)

        for p in primers:
            sa = p.self_annealing()
            sea = p.self_end_annealing()
            vt = [p.name,
                  "5'-%s-3'"%p.seq,
                  len(p),
                  '%.2f'%p.melting_temperature(),
                  '%.2f'%p.melting_temperature(unmethyl=False),
                  tm_gc(p.seq),
                  '%.2f'%p.gc_ratio(),
                  sa.score,
                  sea.score]
            b.tr(style='max-height:1em')
            for i,v in enumerate(vt):
                if width[i]:
                    b.td(str(v),style='width:%sem'%width[i])
                else:
                    b.td(str(v))

def main():
    import sys, os
    from optparse import OptionParser

    parser = OptionParser('usage: %prog inputfile -o outputfile')
    parser.add_option('-o', '--output', dest='output', help='output filename')

    if len(sys.argv) == 0:
        parser.error('no input file')

    (options, args) = parser.parse_args()
    
    inputfiles = args
    
    if options.output:
        sys.stdout = open(options.output, 'w')

    primers = []
    for filename in inputfiles:
        with open(filename,'r') as f:
            print 'loading primers: '+filename,
            for p in load_primer_list_file(f):
                primers.append(p)
                print '.',
            print '...done.'
    
    print 'writing html...'
    html = XmlWriter(sys.stdout)
    b = builder(html)
    with b.html:
        with b.head:
            pass
        with b.body(style='font-family:monospace;font-size:small'):
            b.h1('Primers')
            primers_write_html(html, primers)

if __name__=='__main__':
    main()
