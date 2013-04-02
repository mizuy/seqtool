import os

directory = os.path.dirname(os.path.abspath(__file__))
hgnc_tab_file = os.path.join(directory,'../../_db/hgnc.tab')
ucsc_tab_file = os.path.join(directory,'../../_db/ucsc_gene_table.tab')
                

from sqlalchemy import ForeignKey
from sqlalchemy import Column, Date, Integer, String, Boolean
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.schema import UniqueConstraint
from sqlalchemy.sql import select, or_, and_
from ..util.prompt import prompt
from collections import defaultdict

db_file = os.path.join(os.path.dirname(__file__),'../../_db/seqtool.db')

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
# echo=True for debug
engine = create_engine('sqlite:///'+db_file, echo=False)
Session = sessionmaker(bind=engine)

Base = declarative_base()

t_chromosome = [[x.strip() for x in c.split(',')] for c in """
chr1,  NC_000001, 224589800
chr2,  NC_000002, 224589811
chr3,  NC_000003, 224589815
chr4,  NC_000004, 224589816
chr5,  NC_000005, 224589817
chr6,  NC_000006, 224589818
chr7,  NC_000007, 224589819
chr8,  NC_000008, 224589820
chr9,  NC_000009, 224589821
chr10, NC_000010, 224589801
chr11, NC_000011, 224589802
chr12, NC_000012, 224589803
chr13, NC_000013, 224589804
chr14, NC_000014, 224589805
chr15, NC_000015, 224589806
chr16, NC_000016, 224589807
chr17, NC_000017, 224589808
chr18, NC_000018, 224589809
chr19, NC_000019, 224589810
chr20, NC_000020, 224589812
chr21, NC_000021, 224589813
chr22, NC_000022, 224589814
chrX,  NC_000023, 224589822
chrY,  NC_000024, 224589823
""".strip().split('\n')]

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
    def load(cls):
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
    def load(cls, hgnc_tab_file=hgnc_tab_file):
        session = Session()
        with open(hgnc_tab_file, 'r') as f:
            lines = f.readlines()
            title = lines[0].strip().split('\t')
            assert(title[1]=="Approved Symbol")
            assert(title[2]=="Entrez Gene ID (mapped data supplied by NCBI)")
            for line in lines[1:]:
                l = line.strip().split('\t')
                if len(l) < 4:
                    continue
                try:
                    symbol = l[1]
                    gene_id = int(l[2])
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
    def load(cls, ucsc_tab_file=ucsc_tab_file):
        session = Session()
        with open(ucsc_tab_file, 'r') as f:
            lines = f.readlines()
            title = lines[0].strip().split('\t')
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

def load_all():
    clear_all()
    session = Session()
    print "loading Chromosome..."
    Chromosome.load()
    print "loading GeneTable..."
    GeneTable.load()
    print "loading UcscTable..."
    UcscTable.load()
    print '...done.'
    session.commit()

