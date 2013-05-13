import os

directory = os.path.dirname(os.path.abspath(__file__))
#hgnc_tab_file = os.path.join(directory,'../../_db/hgnc.tab')
#ucsc_tab_file = os.path.join(directory,'../../_db/ucsc_gene_table.tab')

#from sqlalchemy import ForeignKey
from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
#from sqlalchemy.orm import relationship, backref
#from sqlalchemy.schema import UniqueConstraint
#from sqlalchemy.sql import select, or_, and_
#from ..util.prompt import prompt
#from collections import defaultdict

db_file = os.path.join(os.path.dirname(__file__),'../../_db/seqtool.db')

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
# echo=True for debug
engine = create_engine('sqlite:///'+db_file, echo=False)
Session = sessionmaker(bind=engine)

Base = declarative_base()

class Chromosome(Base):
    __tablename__ = 'chromosome'
    id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)
    accession = Column(String, unique=True)
    gid = Column(Integer, unique=True)

    def __init__(self, ch, ac, gid):
        self.name = ch
        self.accession = ac
        self.gid = gid

    @classmethod
    def load(cls, filename):
        chrom = open(filename, 'r').read()
        t_chromosome = [[x.strip() for x in c.split(',')] for c in chrom.strip().split('\n')]
        session = Session()
        for ch, ac, gid in t_chromosome:
            session.add(Chromosome(ch, ac, gid))
        session.commit()

class GeneTable(Base):
    __tablename__ = 'gene_table'
    
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, unique=True) # Gene ID
    symbol = Column(String, unique=True)

    def __init__(self, gene_id, symbol):
        self.gene_id = gene_id
        self.symbol = symbol

    @classmethod
    def load(cls, filename):
        session = Session()
        with open(filename, 'r') as f:
            lines = f.readlines()
            title = lines[0].strip().split('\t')
            assert(title[0]=="Approved Symbol")
            assert(title[1].startswith("Entrez Gene ID"))
            for line in lines[1:]:
                l = line.strip().split('\t')
                if len(l) < 2:
                    continue
                try:
                    symbol = l[0]
                    gene_id = int(l[1])
                except:
                    continue
                session.add(GeneTable(gene_id, symbol))
        session.commit()


class UcscTable(Base):
    __tablename__ = 'ucsc_table'

    id = Column(Integer, primary_key=True)
    n_bin = Column(Integer)
    accession = Column(String)
    chrom = Column(String)
    strand = Column(String)
    txStart = Column(Integer)
    txEnd = Column(Integer)
    cdsStart = Column(Integer)
    cdsEnd = Column(Integer)
    exonCount = Column(Integer)
    exonStarts  = Column(String)
    exonEnds = Column(String)
    score = Column(Integer)
    symbol = Column(String)
    cdsStartStat = Column(String)
    cdsEndStat = Column(String)
    exonFrames  = Column(String)

    def __init__(self):
        pass

    @classmethod
    def load(cls, filename):
        session = Session()
        with open(filename, 'r') as f:
            lines = f.readlines()
            #title = lines[0].strip().split('\t')
            for line in lines[1:]:
                l = line.strip().split('\t')
                assert len(l) == 16
                
                u = UcscTable()
                u.n_bin = l[0]
                u.accession = l[1]
                u.chrom = l[2]
                u.strand = l[3]
                u.txStart = int(l[4])
                u.txEnd = int(l[5])
                u.cdsStart = int(l[6])
                u.cdsEnd = int(l[7])
                u.exonCount = int(l[8])
                u.exonStarts = l[9]
                u.exonEnds = l[10]
                u.score = int(l[11])
                u.symbol = l[12]
                u.cdsStartStat = l[13]
                u.cdsEndStat = l[14]
                u.exonFrames = l[15]
                
                session.add(u)
        session.commit()

Base.metadata.create_all(engine)

def clear_all():
    con = engine.connect() 

    trans = con.begin() 

    for name, table in Base.metadata.tables.items(): 
        print table.delete() 
        con.execute(table.delete()) 

    trans.commit() 

def database_load():
    from argparse import ArgumentParser
    parser = ArgumentParser(prog='database_load', description='load database')
    parser.add_argument("--ucsc_tab_file", default='ucsc.tab', help='see Makefile')
    parser.add_argument("--hgnc_tab_file", default='hgnc.tab', help='see Makefile')
    parser.add_argument("--chrom_tab_file", default='chrom.tab', help='see Makefile')
    args = parser.parse_args()

    clear_all()
    session = Session()
    print "loading Chromosome..."
    Chromosome.load(args.chrom_tab_file)
    print "loading GeneTable..."
    GeneTable.load(args.hgnc_tab_file)
    print "loading UcscTable..."
    UcscTable.load(args.ucsc_tab_file)
    print '...done.'
    session.commit()

