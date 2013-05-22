import sys, os, glob
import traceback
from argparse import ArgumentParser, FileType

from ..util.dirutils import Filepath
from ..view import seqview as seqv
from ..util import xmlwriter
from ..nucleotide.primer import Primers
from contextlib import contextmanager
from ..script.sequencing import SequencingAnalysis
import pdb

from .configuration import GeneralConfiguration

from .. import db

def log(message):
    sys.stderr.write('{0}\n'.format(message))
    sys.stderr.flush()

@contextmanager
def report_exceptions():
    try:
        yield
    except Exception:
        e, m, tb = sys.exc_info()
        print('exception traceback:'.ljust( 80, '=' ))
        for tbi in traceback.format_tb( tb ):
            print(tbi)
        print('  %s' % str( m ))
        print(''.rjust( 80, '=' ))
        pdb.post_mortem(tb)

def init_db():
    gc = GeneralConfiguration()
    db.initialize(gc.get_cache_dir(), gc.get_email())

def get_output_path(input, output, multiple_inputs, ext='.html'):
    """
    input and output option policy

    - IF no output option specified, use input file directory as output directory.
    - IF output options specified, output options are directory or filename.
    - IF multiple input files, output option must be directory or ignored.
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
            print("If you set multiple input files, output option must be existing directroy OR ignored.")
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
            print("processing input: {0}, output:{1}".format(inputp.path, outputp.path))

            init_db()
            p = seqv.Seqviews()
            p.load_seqv(inputp.path)
            p.write_html(outputp)

            outputs.append(outputp.path)

    if outputs and args.open:
        os.system('open {0}'.format(' '.join(outputs)))

def tssview():
    print('not implemented.')
    return 
    
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
            print("processing input: {0}, output:{1}".format(inputp.path, outputp.path))

            init_db()
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
                gene_id, gene_symbol = db.genetable.get_gene_from_text(gene_symbol)
            except db.NoSuchGene:
                log('gene entry: No such Gene %s'%gene_symbol)
                continue

            init_db()
            sv = seqv.Seqviews()
            e = sv.load_gene(gene_symbol, gene_id)
            e.dbtss.set_tissueset(db.dbtss.get_dbtss_tissues())

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
        print(' processing:',sf)
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

    init_db()
    genbank = db.genetable.get_genomic_context_genbank(args.gene)

    if not genbank:
        print("GenBank retrieve error: ", args.gene)
        return
    print(genbank, file=args.output)

def primers():
    parser = ArgumentParser(prog='primers', description='Print primer properties')
    parser.add_argument('primers', nargs='+', help='primer text file')
    parser.add_argument('-c', '--csv', action='store_true', dest='csv', help='write a csv file(default is html file)')
    parser.add_argument("-o", "--output", dest="output", nargs='?',
                         type=FileType('w'), default=sys.stdout, help="output filename")

    args = parser.parse_args()

    inputfiles = args.primers
    
    ps = Primers()

    for filename in inputfiles:
        ps.load_file(filename)

    if args.csv:
        args.output.write(ps.write_csv())
    else:
        html = xmlwriter.XmlWriter(args.output)
        b = xmlwriter.builder(html)
        with b.html:
            with b.head:
                pass
            with b.body(style='font-family:monospace;font-size:small'):
                b.h1('Primers')
                ps.write_html(html)

def seqdb_command():
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='seqdb', description='seqtool database administration.')
    parser.add_argument("--cache_dir", help='database directory')

    subparsers = parser.add_subparsers(help='sub-command help')
    parser_load = subparsers.add_parser('load', help='load database')
    parser_load.add_argument("--ucsc_tab_file", default='ucsc.tab', help='see Makefile')
    parser_load.add_argument("--hgnc_tab_file", default='hgnc.tab', help='see Makefile')
    parser_load.add_argument("--chrom_tab_file", default='chrom.tab', help='see Makefile')
    parser_load.set_defaults(func=seqdb_load)

    parser_dbtss = subparsers.add_parser('dbtss', help='dbtss database')
    parser_dbtss.add_argument("--bed_dir", help='bed directory')
    parser_dbtss.set_defaults(func=seqdb_dbtss)

    parser_clear = subparsers.add_parser('clear', help='clear database')
    parser_clear.set_defaults(func=seqdb_clear)

    if (len(sys.argv) < 2):
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args(sys.argv[1:])

    if not args.cache_dir:
        gc = GeneralConfiguration()
        args.cache_dir = gc.get_cache_dir()

    args.func(args)


def seqdb_clear(args):
    db.clear(args.cache_dir)
    print("genetable cleared.")

def seqdb_dbtss(args):
    print("loading dbtss...")
    db.load_dbtss(bed_dir=args.bed_dir, cache_dir=args.cache_dir)

def seqdb_load(args):
    print("loading gene tables...")
    db.load_table(args.cache_dir, args.chrom_tab_file, args.hgnc_tab_file, args.ucsc_tab_file)

