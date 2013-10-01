import os
from .locus import Locus
#from sqlalchemy import ForeignKey
from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
#from sqlalchemy.orm import relationship, backref
#from sqlalchemy.schema import UniqueConstraint
#from sqlalchemy.sql import select, or_, and_

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
# echo=True for debug

from . import entrez

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
    def load(cls, filename, session):
        chrom = open(filename, 'r').read()
        t_chromosome = [[x.strip() for x in c.split(',')] for c in chrom.strip().split('\n')]
        for ch, ac, gid in t_chromosome:
            session.add(Chromosome(ch, ac, gid))

class HgncTable(Base):
    __tablename__ = 'hgnc_table'
    
    id = Column(Integer, primary_key=True)
    gene_id = Column(Integer, unique=True) # Gene ID
    symbol = Column(String, unique=True)

    def __init__(self, gene_id, symbol):
        self.gene_id = gene_id
        self.symbol = symbol

    @classmethod
    def load(cls, filename, session):
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
                session.add(HgncTable(gene_id, symbol))


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
    def load(cls, filename, session):
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

def make_engine(cache_dir):
    database_file = os.path.join(cache_dir,'seqtool.db')
    engine = create_engine('sqlite:///'+database_file, echo=False)
    Base.metadata.create_all(engine)
    return engine


class GeneTable(object):
    def __init__(self, cache_dir, email=None):
        self.engine = make_engine(cache_dir)
        self.Session = sessionmaker(bind=self.engine)

        if email:
            self.entrez = entrez.CachedEntrez(cache_dir, email)

    def clear(self):
        Base.metadata.create_all(self.engine)
        con = self.engine.connect()

        trans = con.begin()

        for name, table in list(Base.metadata.tables.items()): 
            print(table.delete())
            con.execute(table.delete())

        trans.commit() 

        print("all database cleared.")

    def load(self, chrom_tab_file, hgnc_tab_file, ucsc_tab_file):
        self.clear()

        session = self.Session()
        print("loading Chromosome...")
        Chromosome.load(chrom_tab_file, session)
        session.commit()

        print("loading HgncTable...")
        HgncTable.load(hgnc_tab_file, session)
        session.commit()

        print("loading UcscTable...")
        UcscTable.load(ucsc_tab_file, session)
        session.commit()

        print('...done.')


    def get_gene_symbol(self, gene_id):
        session = self.Session()
        try:
            return session.query(HgncTable.symbol).filter_by(gene_id=gene_id).one()[0]
        except:
            return None

    def get_gene_id(self, gene_symbol):
        session = self.Session()
        try:
            return session.query(HgncTable.gene_id).filter_by(symbol=gene_symbol).one()[0]
        except:
            return None

    def get_gene_from_text(self, text):
        try:
            gene_id = int(text)
            gene_symbol = self.get_gene_symbol(gene_id)
        except ValueError:
            gene_symbol = text
            gene_id = self.get_gene_id(gene_symbol)
        return gene_id, gene_symbol

    def get_gene_locus(self, gene_id):
        symbol = self.get_gene_symbol(gene_id)
        if not symbol:
            return None
        session = self.Session()
        # TODO: gene has multiple mRNA.
        # TODO: gene might have even multiple loci !!
        try:
            f = session.query(UcscTable).filter_by(symbol=symbol).first()
            return Locus(f.chrom, f.strand=='+', f.txStart, f.txEnd)
        except:
            return None

    def get_chrom_gid(self, locus):
        session = self.Session()

        try:
            return session.query(Chromosome.gid).filter_by(name=locus.chrom).one()[0]
        except:
            return None

    def get_locus_genbank(self, locus):
        if not self.entrez:
            print('entrez. not loaded.')
            return None
        session = self.Session()
        try:
            chrom_gid = session.query(Chromosome.gid).filter_by(name=locus.chrom).one()[0]
            return self.entrez.get_genbank(chrom_gid, locus)
        except:
            return None

    def get_gene_genbank(self, gene_id):
        locus = self.get_gene_locus(gene_id)
        return self.get_locus_genbank(locus)

    def get_genomic_context_genbank(self, gene_text, upstream=1000, downstream=1000):
        gene_id, gene_symbol = self.get_gene_from_text(gene_text)
        locus = self.get_gene_locus(gene_id).expand(upstream, downstream)
        return self.get_locus_genbank(locus)

if __name__=='__main__':
    db = GeneTable('/Users/mizuy/tmp/')
    print(db.get_genomic_context_genbank('MAL'))
