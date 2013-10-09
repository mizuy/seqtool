import sys, os, glob
from argparse import ArgumentParser
from collections import defaultdict
from ..util.dirutils import Filepath
from ..util.debug import report_exceptions
from ..util import DefaultOrderedDict

def log(message):
    sys.stderr.write('{0}\n'.format(message))
    sys.stderr.flush()

def init_db():
    from .. import db
    from .configuration import GeneralConfiguration
    gc = GeneralConfiguration()
    db.initialize(gc.get_cache_dir(), gc.get_email())

'''
import pkgutil
def module_mtime(module):
    pp = pkgutil.walk_packages(module.__path__, prefix=module.__name__+'.')
    for i in pp:
        print(i[0])
    return max(os.path.getmtime(i[0].path) for i in pp)
'''

def check_update(dest_path, source_paths, force_update):
    if force_update:
        return True

    def _update(dest_path, source_paths):
        if not os.path.exists(dest_path):
            return True

        dst_mtime = os.path.getmtime(dest_path)

        srcs = [os.path.getmtime(s) for s in source_paths if os.path.exists(s)]
        if not srcs:
            return False

        if max(srcs) > dst_mtime:
            return True
        else:
            return False

    if _update(dest_path, source_paths):
        print("Need update. output:{}".format(dest_path))
        return True
    else:
        print("Nothing to be done for {}".format(dest_path))
        #print("set -f if you want a force update.")
        return False

def get_output_path(input, output, ext='.html'):
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
    else:
        if os.path.isdir(output):
            return Filepath(os.path.join(output, default_output) )
        else:
            return Filepath(output)

def mac_open(outputp):
    if outputp:
        command = 'open {}'.format(outputp.path)
        print(command)
        os.system(command)

def seqview():
    parser = ArgumentParser(prog='seqview', description='pretty HTML+SVG report of sequence and adittional data')
    parser.add_argument("input", help=".seqv file")
    parser.add_argument("-o", "--output", dest="output", help="output html filename or directory")
    parser.add_argument("--open", dest="open", help="open output file using Mac command 'open'", action='store_true')
    parser.add_argument("-f", "--force", dest="force", help="force update", action='store_true')

    args = parser.parse_args()
    inputp = Filepath(args.input)
    outputp = get_output_path(args.input, args.output, ext='.html')


    if check_update(outputp.path, [inputp.path], args.force):
        init_db()
        print("processing input:{}, output:{}".format(inputp.path,outputp.path))
        from ..view import seqview as seqv
        p = seqv.Seqview.load_seqv(inputp.path)
        p.write_html(outputp)

        if args.open:
            mac_open(outputp)


def geneview():
    parser = ArgumentParser(prog='geneview', description='pretty HTML+SVG report of gene')
    parser.add_argument("gene_symbol", help="Gene Symbols")
    parser.add_argument("-o", "--output", dest="output", help="output filename or directory")
    parser.add_argument("--open", dest="open", help="open output file using Mac command 'open'", action='store_true')

    args = parser.parse_args()
    gene_symbol = args.gene_symbol
    outputp = get_output_path(gene_symbol, args.output, False)

    init_db()
    try:
        from .. import db
        genbank = db.instance.get_genbank(gene_symbol)
    except db.NoSuchGene:
        log('gene entry: No such Gene %s'%gene_symbol)
        return

    from ..view import seqview as seqv
    sv = seqv.Seqview.create_gene(gene_symbol)
    sv.add_dbtss(db.instance.beddb.get_dbtss_tissues())
    sv.write_html(outputp)

    if args.open:
        mac_open(outputp)

def sequencing():
    parser = ArgumentParser(prog='sequencing', description='sequencing result aligner')
    parser.add_argument('inputs', nargs='*', help=".seq files.")
    parser.add_argument("-t", "--template", dest="template", help="template fasta file.")
    parser.add_argument("-o", "--output", dest="output", help="output html filename")
    parser.add_argument("--open", dest="open", help="open output file using Mac command 'open'", action='store_true')
    parser.add_argument("-f", "--force", dest="force", help="force update", action='store_true')
    parser.add_argument("-d", "--default", dest="default", help='set default arguments', action='store_true')
    args = parser.parse_args()

    if args.default:
        args.template = 'reference.fasta'
        args.inputs = '.'

    if not args.template:
        log('No template file.')
        return

    if not args.output:
        args.output = 'aligned.html'

    if not os.path.exists(args.template):
        print('No such file: {}'.format(args.template))
        return

    outputp = Filepath(args.output)
        
    srcs = DefaultOrderedDict(list)
    inputs = []

    def sf(f):
        n,ext = os.path.splitext(f)
        return ext in ['.ab1']
        #return ext in ['.seq','.ab1']
        
    for input in args.inputs:
        if os.path.isdir(input):
            for dpath, dnames, fnames in os.walk(input):
                dnames.sort()
                fnames.sort()
                #if os.path.basename(dpath) == 'trash':
                #    continue
                for f in fnames:
                    if sf(f):
                        srcs[dpath].append(Filepath(os.path.join(dpath,f)).path)
        else:
            if sf(input):
                srcs['top'].append(Filepath(input).path)

    inputps = [Filepath(i).path for i in sum(srcs.values(), [])]
            
    if check_update(outputp.path, inputps, args.force):
        with report_exceptions():
            from ..script.sequencing import sequencing_alignment_dir
            sequencing_alignment_dir(args.template, srcs, outputp)

            if args.open:
                mac_open(outputp)
                
                
def abiview():
    parser = ArgumentParser(prog='abiview', description='ABI sequencing output viewer')
    parser.add_argument('ab1_filename', help=".ab1 file.")
    parser.add_argument("-o", "--output", dest="output", help="output html filename")
    parser.add_argument("--open", dest="open", help="open output file using Mac command 'open'", action='store_true')
    parser.add_argument("-f", "--force", dest="force", help="force update", action='store_true')

    args = parser.parse_args()

    inputp = Filepath(args.ab1_filename)
    outputp = get_output_path(inputp.path, args.output, ext='.svg')
    
    if check_update(outputp.path, [inputp.path], args.force):
        with report_exceptions():
            from ..format.abi import write_svg
            write_svg(inputp.path, outputp.path)

            if args.open:
                mac_open(outputp)
                
def get_genbank():
    parser = ArgumentParser(prog='get_genbank', description='Retrieve Genbank from Gene ID or Gene Symbol')
    parser.add_argument("gene_symbol", help="Gene Symbol")

    args = parser.parse_args()

    init_db()
    from .. import db
    genbank = db.instance.get_genbank(args.gene_symbol)

    if not genbank:
        print("Failed to retrieve GenBank: ", args.gene_symbol)
        return
    sys.stdout.write(genbank)


def primers():
    parser = ArgumentParser(prog='primers', description='Print primer properties')
    parser.add_argument('primers', nargs='+', help='primer text file')
    parser.add_argument('-c', '--csv', action='store_true', dest='csv', help='write a csv file(default is html file)')

    args = parser.parse_args()

    inputfiles = args.primers
    
    from ..nucleotide.primer import Primers
    from ..util import xmlwriter
    ps = Primers()

    for filename in inputfiles:
        ps.load_file(filename)

    if args.csv:
        sys.stdout.write(ps.write_csv())
    else:
        html = xmlwriter.XmlWriter(sys.stdout)
        b = xmlwriter.builder(html)
        with b.html:
            with b.head:
                pass
            with b.body(style='font-family:monospace;font-size:small'):
                b.h1('Primers')
                ps.write_html(html)

