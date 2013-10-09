import os
from .locus import Locus
from sqlalchemy import Column, Integer, String
from sqlalchemy.ext.declarative import declarative_base
#from sqlalchemy import ForeignKey
#from sqlalchemy.orm import relationship, backref
#from sqlalchemy.schema import UniqueConstraint
#from sqlalchemy.sql import select, or_, and_

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from ..util import prompt

def windowed_query(q, session, limit=8000):
    offset = 0
    while True:
        elem = None
        for elem in session.execute(q.offset(offset).limit(limit)):
            yield elem
        offset += limit
        if elem is None:
            break

class batch_insert:
    def __init__(self, insert, session, batchsize=20000, prompt=None):
        self._records = []
        self._insert = insert
        self._session = session
        self._batchsize = batchsize
        self._progress = prompt

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._records:
            self.execute()
        
    def insert(self, record):
        self._records.append(record)
        if len(self._records) >= self._batchsize:
            self.execute()

    def execute(self):
        self._session.execute(self._insert, self._records)
        self._session.commit()
        self._records = []
        if self._progress:
            self._progress.progress()
        
            
Base = declarative_base()

class knownGene(Base):
    __tablename__ = 'knownGene'
    name = Column(String, unique=True, primary_key=True)
    chrom = Column(String)
    strand = Column(String)
    txStart = Column(Integer)
    txEnd = Column(Integer)
    cdsStart = Column(Integer)
    cdsEnd = Column(Integer)
    exonCount = Column(Integer)
    exonStarts = Column(String)
    exonEnds = Column(String)
    proteinID = Column(String)
    alignID = Column(String)

    def __repr__(self):
        return '<knownGene({}, {}, {}, {}, {})>'.format(self.name, self.chrom, self.strand, self.txStart, self.txEnd)

    @property
    def locus(self):
        return Locus(self.chrom, self.strand=='+', self.txStart, self.txEnd)

class kgXref(Base):
    __tablename__ = 'kgXref'
    kgID = Column(String, primary_key=True)
    mRNA = Column(String)
    spID = Column(String)
    spDisplayID = Column(String)
    geneSymbol = Column(String)
    refseq = Column(String)
    protAcc = Column(String)
    description = Column(String)
    rfamAcc = Column(String)
    tRnaName = Column(String)


# http://code.activestate.com/recipes/576644-diff-two-dictionaries/
def diff_dict(d1, d2, NO_KEY='<KEYNOTFOUND>'):
    both = d1.keys() & d2.keys()
    diff = {k:(d1[k], d2[k]) for k in both if d1[k] != d2[k]}
    diff.update({k:(d1[k], NO_KEY) for k in d1.keys() - both})
    diff.update({k:(NO_KEY, d2[k]) for k in d2.keys() - both})
    return diff

    
class GenomeDB:
    def __init__(self, connection_string, echo=False):
        print('connecting: {}'.format(connection_string))
        self.connection_string = connection_string
        self.engine = create_engine(connection_string, echo=echo)
        Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)

    def get_knownGene(self, gene_symbol):
        session = self.Session()
        for xr in session.query(kgXref).filter_by(geneSymbol=gene_symbol).all():
            print('knowngene id: ', xr.kgID, xr.mRNA, xr.spID, xr.geneSymbol, xr.refseq)
            print(session.query(knownGene).filter_by(name=xr.kgID).one())
            pass
        
    def get_gene_locus(self, gene_symbol):
        session = self.Session()
        xr = session.query(kgXref).filter_by(geneSymbol=gene_symbol).first()
        if not xr:
            return None
        kg = session.query(knownGene).filter_by(name=xr.kgID).one()
        if not kg:
            return None
        return Locus(kg.chrom, kg.strand=='+', kg.txStart, kg.txEnd)

    def drop_all(self):
        for table in Base.metadata.sorted_tables:
            table.drop(self.engine, checkfirst=True)
        Base.metadata.create_all(self.engine)

    @classmethod
    def load_mysql(self, user='genome', host='localhost', db='hg19', password=None):
        if password:
            user = '{}:{}'.format(user, password)
        return GenomeDB('mysql+mysqlconnector://{}@{}/{}'.format(user,host,db))

    @classmethod
    def load_sqlite(self, cache_dir):
        cache_dir = os.path.abspath(os.path.expanduser(cache_dir))
        db = GenomeDB('sqlite:///'+os.path.join(cache_dir,'seqtool.ucsc.db'))
        return db

    @classmethod
    def copy_db(cls, src, dst, force=False):
        ssession = src.Session()
        dsession = dst.Session()

        if not force:
            if src.rows() == dst.rows():
                return

        dst.drop_all()

        for table in Base.metadata.sorted_tables:
            columns = table.columns.keys()

            with prompt.prompt('Copying {}'.format(table.name)) as progress:
                with batch_insert(table.insert(), dsession, prompt=progress) as ins:
                    ts = table.select()
                    for record in windowed_query(ts, ssession):
                        a = {str(col) : getattr(record, col) for col in columns }
                        ins.insert(a)

        dd = diff_dict(src.rows(), dst.rows())
        if dd:
            raise Exception("Mirroed Database difference error: {}'".format(dd))

    def rows(self):
        session = self.Session()
        return {table.name : session.query(table).count() for table in Base.metadata.sorted_tables}

if __name__=='__main__':
    import sys
    password = sys.argv[1]
    
    for table in Base.metadata.sorted_tables:
        print (table.name)
        
    src = GenomeDB.load_mysql('genome', password=password, host='localhost')
    print('mysql')
    print(src.get_gene_locus('B2M'))
    dst = GenomeDB.load_sqlite('~/tmp/')
    GenomeDB.copy_db(src, dst)
    print('mirror')
    print(src.get_gene_locus('B2M'))


    
