from __future__ import absolute_import

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC

from .seqview import load_primer_list_file
from .htmlwriter import HtmlWriter
from .nucleotide import tm_gc

def primer_table_html(html, primers):
    pt = ['name', 'sequence', 'len[bp]', 'Tm[C]', 'oTm[C]', 'old Tm[C]','GC[%]', 'sa', 'sea']
    width = [16, 26, None, None, None, None, None, None, None, None]

    html.push('table',cls='primer',border=1)
    html.push('tr',style='max-height:1em')
    for p in pt:
        html.insert('th',p)
    html.pop()

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
        html.push('tr',style='max-height:1em')
        for i,v in enumerate(vt):
            if width[i]:
                html.insert('td',str(v),style='width:%sem'%width[i])
            else:
                html.insert('td',str(v))
        html.pop()
    html.pop()

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
    html = HtmlWriter(sys.stdout)
    html.push('html')
    html.push('head')
    html.pop()
    html.push('body',style='font-family:monospace;font-size:small')

    html.insert('h1','Primers')
    primer_table_html(html, primers)

    html.pop()
    html.pop()

if __name__=='__main__':
    main()
