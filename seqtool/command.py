from __future__ import absolute_import

import sys, os, glob
import traceback
from argparse import ArgumentParser, FileType

from .util.dirutils import Filepath
from .view import seqview as seqv
from . import db
from .util import xmlwriter
from .nucleotide.pcr import primers_write_html, load_primer_list_file, primers_write_csv
from contextlib import contextmanager
from .script.sequencing import SequencingAnalysis

def log(message):
    sys.stderr.write('{0}\n'.format(message))
    sys.stderr.flush()

@contextmanager
def report_exceptions():
    try:
        yield
    except Exception:
        info = sys.exc_info()
        tbinfo = traceback.format_tb( info[2] )
        print 'exception traceback:'.ljust( 80, '=' )
        for tbi in tbinfo:
            print tbi
        print '  %s' % str( info[1] )
        print ''.rjust( 80, '=' )

def get_genomic_context_genbank(gene_text, upstream=1000, downstream=1000):
    gene_id, gene_symbol = db.get_gene_from_text(gene_text)
    locus = db.get_gene_locus(gene_id).expand(upstream, downstream)
    return db.get_locus_genbank(locus)

def get_output_path(input, output, multiple_inputs, ext='.html'):
    """
    input and output option policy

    - IF no output option specified, use input file directory as output directory.
    - IF single input file, output option is considered as output filename.
    - IF multiple input files, output option is considered as directory name.
    - for geneview and get_genebank, which has no inputfile, use currentdirectory as default output directory.

    """
    prefix, suffix = os.path.splitext(os.path.basename(input))
    default_output = prefix + ext

    if not output:
        return Filepath(input).change_ext(ext)
    elif not multiple_inputs:
        if os.path.isdir(output):
            return Filepath(os.path.join(output, default_output) )
        else:
            return Filepath(output)
    else:
        if os.path.isdir(output):
            return Filepath(os.path.join(output, default_output) )
        else:
            print "If you set multiple input files, output option must be existing directroy OR ignored."
            return Filepath(input).change_ext(ext)

def seqview():
    parser = ArgumentParser(prog='seqview', description='pretty HTML+SVG report of sequence and adittional data')
    parser.add_argument("inputs", nargs='+', metavar="seqvfiles", help=".seqv files")
    parser.add_argument("-o", "--output", dest="output", help="output html filename or directory")
    parser.add_argument("--open", dest="open", help="open output file using Mac command 'open'", action='store_true')

    args = parser.parse_args()

    outputs = []

    for input in args.inputs:
        with report_exceptions():
            inputp = Filepath(input)
            outputp = get_output_path(input, args.output, len(args.inputs)>1)
            print "processing input: {0}, output:{1}".format(inputp.path, outputp.path)

            p = seqv.Seqviews()
            p.load_seqv(inputp.path)
            p.write_html(outputp)

            outputs.append(outputp.path)

    if outputs and args.open:
        os.system('open {0}'.format(' '.join(outputs)))

def tssview():
    parser = ArgumentParser(prog='tssview', description='pretty HTML+SVG report of dbtss of multiple genes')
    parser.add_argument("inputs", nargs='+', metavar='tssviewfiles', help=".tssv files")
    parser.add_argument("-o", "--output", dest="output", help="output filename or directory")
    parser.add_argument("--open", dest="open", help="open output file using Mac command 'open'", action='store_true')

    args = parser.parse_args()

    outputs = []

    for input in args.inputs:
        with report_exceptions():
            inputp = Filepath(input)
            outputp = get_output_path(input, args.output, len(args.inputs)>1)
            print "processing input: {0}, output:{1}".format(inputp.path, outputp.path)

            sv = seqv.TssvFile()
            sv.load_tssv(inputp.path)
            sv.write_csv(outputp.change_ext('.csv').path)
            sv.write_html(outputp)

            outputs.append(outputp.path)

    if outputs and args.open:
        os.system('open {0}'.format(' '.join(outputs)))

def geneview():
    parser = ArgumentParser(prog='geneview', description='pretty HTML+SVG report of gene')
    parser.add_argument("gene_symbols", nargs='+', help="Gene Symbols")
    parser.add_argument("-o", "--output", dest="output", help="output filename")
    parser.add_argument("--open", dest="open", help="open output file using Mac command 'open'", action='store_true')

    args = parser.parse_args()

    outputs = []

    for gene_symbol in args.gene_symbols:
        with report_exceptions():
            try:
                gene_id, gene_symbol = db.get_gene_from_text(gene_symbol)
            except db.NoSuchGene:
                log('gene entry: No such Gene %s'%gene_symbol)
                continue

            sv = seqv.Seqviews()
            e = sv.load_gene(gene_symbol, gene_id)
            e.dbtss.set_tissueset(db.get_dbtss_tissues())

            outputp = get_output_path(gene_symbol, args.output, len(args.gene_symbols)>1)

            sv.write_html(outputp)

            outputs.append(outputp.path)

    if outputs and args.open:
        os.system('open {0}'.format(' '.join(outputs)))


def sequencing():
    parser = ArgumentParser(prog='sequencing', description='sequencing result aligner')
    parser.add_argument('input', help="directory contains .seq files.")
    parser.add_argument("-t", "--template", dest="template", help="template fasta file.")
    parser.add_argument("-o", "--output", dest="output", help="output html filename or directory")
    parser.add_argument("--open", dest="open", help="open output file using Mac command 'open'", action='store_true')

    args = parser.parse_args()
    assert(os.path.isdir(args.input))
    p = SequencingAnalysis()
    p.load_fasta(args.template)

    if args.output:
        if os.path.isdir(args.output):
            outputp = Filepath(os.path.join(args.output, 'sequencing_analysis.html') )
        else:
            outputp = Filepath(args.output)
    else:
        outputp = Filepath(os.path.join(args.input, 'sequencing_analysis.html') )

    for sf in glob.glob(os.path.join(args.input,'*.seq')):
        print ' processing:',sf
        p.load_seqfile(sf)
    
    p.write_html(outputp)

    if args.open:
        os.system('open {0}'.format(outputp.path))


def get_genbank():
    parser = ArgumentParser(prog='get_genbank', description='Retrieve Genbank from Gene ID or Gene Symbol')
    parser.add_argument("gene", help="Gene ID or Gene Symbol")
    parser.add_argument("-o", "--output", dest="output", nargs='?',
                         type=FileType('w'), default=sys.stdout, help="output filename")

    args = parser.parse_args()
    genbank = get_genomic_context_genbank(args.gene)

    if not genbank:
        print "GenBank retrieve error: ", args.gene
        return
    print >>args.output, genbank

def primers():
    parser = ArgumentParser(prog='primers', description='Print primer properties')
    parser.add_argument('primers', nargs='+', help='primer text file')
    parser.add_argument('-c', '--csv', action='store_true', dest='csv', help='write a csv file(default is html file)')
    parser.add_argument("-o", "--output", dest="output", nargs='?',
                         type=FileType('w'), default=sys.stdout, help="output filename")

    args = parser.parse_args()

    inputfiles = args.primers
    
    ps = []
    for filename in inputfiles:
        with open(filename,'r') as f:
            for i in load_primer_list_file(f):
                ps.append(i)

    if args.csv:
        args.output.write(primers_write_csv(ps))
    else:
        html = xmlwriter.XmlWriter(args.output)
        b = xmlwriter.builder(html)
        with b.html:
            with b.head:
                pass
            with b.body(style='font-family:monospace;font-size:small'):
                b.h1('Primers')
                primers_write_html(html, ps)


