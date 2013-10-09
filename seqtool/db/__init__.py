from .locus import Locus
from .bed import BedDB
from .ucsc import GenomeDB
from .entrez import CachedEntrez

__all__ = ['DB', 'instance', 'initialize', 'NoSuchGene', 'seqdb_command']

class NoSuchGene(Exception):
    def __init__(self, m):
        self.gene_text = str(m)
    def __str__(self):
        return 'No Such Gene "%s"'%self.gene_text


instance = None
def initialize(cache_dir, email):
    global instance
    instance = DB(cache_dir, email)
    
class DB:
    def __init__(self, cached_dir, email):
        self.cache_dir = cached_dir
        self.email = email

        self.beddb = BedDB()
        self.beddb.load_cache(self.cache_dir)

        self.genomedb = GenomeDB.load_sqlite(self.cache_dir)
        self.centrez = CachedEntrez(self.cache_dir, email)

    def clear(self):
        self.genomedb = GenomeDB.load_sqlite(self.cache_dir)
        self.genomedb.drop_all()
        """
        beddb = BedDB()
        beddb.load_file(args.bed_dir, args.cache_dir)
        beddb.clear()
        """

    def load_bed(self, bed_dir):
        print('loading bed files from {}'.format(bed_dir))
        self.beddb.load_file(bed_dir, self.cache_dir)

    def load_table(self, host='localhost', password=None, force=False):
        print('loading sql tables from {}'.format(host))
        mysql = GenomeDB.load_mysql('genome', host=host, password=password)
        GenomeDB.copy_db(mysql, self.genomedb, force=force)

    def get_genbank(self, gene_symbol):
        return self.centrez.get_genbank(self.genomedb.get_gene_locus(gene_symbol))
        
def seqdb_command():
    from argparse import ArgumentParser
    import sys
    from ..frontend.configuration import GeneralConfiguration
    parser = ArgumentParser(prog='seqdb', description='seqtool database administration.')
    parser.add_argument("--cache_dir", help='database directory')

    subparsers = parser.add_subparsers(help='sub-command help')
    parser_load = subparsers.add_parser('load', help='load database from mysql server')
    parser_load.add_argument("--host", default='localhost', help='see Makefile')
    parser_load.add_argument("--password", help='see Makefile')
    parser_load.add_argument("-f", "--force", dest="force", help="force update", action='store_true')
    parser_load.set_defaults(func=seqdb_load)

    parser_bed = subparsers.add_parser('bed', help='load bed files.')
    parser_bed.add_argument("--bed_dir", help='bed directory')
    parser_bed.set_defaults(func=seqdb_bed)

    parser_clear = subparsers.add_parser('clear', help='clear database')
    parser_clear.set_defaults(func=seqdb_clear)

    if (len(sys.argv) < 2):
        args = parser.parse_args(['-h'])
    else:
        args = parser.parse_args(sys.argv[1:])

    if not args.cache_dir:
        gc = GeneralConfiguration()
        args.cache_dir = gc.get_cache_dir()

    args.instance = DB(args.cache_dir, None)

    args.func(args)

def seqdb_clear(args):
    args.instance.clear()
    print("genetable cleared.")

def seqdb_bed(args):
    print("loading bed...")
    args.instance.load_bed(bed_dir=args.bed_dir)

def seqdb_load(args):
    print("loading gene tables...")
    args.instance.load_table(host=args.host, password=args.password, force=args.force)

